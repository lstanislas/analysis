#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <TMath.h>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"

using o2::mch::Digit;

//_________________________________________________________________________________________________
std::function<double(double, double, double)> MathiesonIntegrate(bool fixK3, double k3 = 0.3)
{
  auto init = [](double k3) -> std::tuple<double, double, double> {
    double sqrtK3 = std::sqrt(k3);
    double k2 = (TMath::Pi() / 2.) * (1. - 0.5 * sqrtK3);
    double k4 = 1. / (4 * TMath::ATan(sqrtK3));
    return std::make_tuple(sqrtK3, k2, k4);
  };

  auto integrate = [](double min, double max, double sqrtK3, double k2, double k4) -> double {
    double umin = sqrtK3 * TMath::TanH(k2 * min);
    double umax = sqrtK3 * TMath::TanH(k2 * max);
    return 2. * k4 * (TMath::ATan(umax) - TMath::ATan(umin));
  };

  if (fixK3) {
    auto [sqrtK3, k2, k4] = init(k3);
    return [sqrtK3, k2, k4, integrate](double min, double max, double) -> double {
      return integrate(min, max, sqrtK3, k2, k4);
    };
  } else {
    return [init, integrate](double min, double max, double k3) -> double {
      auto [sqrtK3, k2, k4] = init(k3);
      return integrate(min, max, sqrtK3, k2, k4);
    };
  }
}

//_________________________________________________________________________________________________
std::function<double(double, double, double)> Error(std::string mode, double alpha)
{
  if (mode == "LS") {
    return [alpha](double, double, double adcFit) -> double {
      return alpha * std::sqrt(adcFit);
    };
  } else if (mode == "MLS") {
    return [alpha](double adc, double, double) -> double {
      return alpha * std::sqrt(adc);
    };
  } else if (mode == "MC") {
    return [alpha](double, double sample, double) -> double {
      return 0.5 * (std::sqrt(sample) + alpha);
    };
  } else if (mode.starts_with("MULT_")) {
    auto params = mode.substr(5);
    std::replace(params.begin(), params.end(), 'p', '.');
    double a = 1.0, b = 0., g = 0.;
    size_t pos1 = params.find('_');
    size_t pos2 = (pos1 == std::string::npos) ? std::string::npos : params.find('_', pos1 + 1);
    if (pos1 != std::string::npos && pos2 != std::string::npos) {
      try {
        a = std::stod(params.substr(0, pos1));
        b = std::stod(params.substr(pos1 + 1, pos2 - pos1 - 1));
        g = std::stod(params.substr(pos2 + 1));
      } catch (...) {
      }
    }
    return [a, b, g](double adc, double, double) -> double {
      const double sqrtQ = std::sqrt(adc);
      return a * sqrtQ + b * adc + g * adc * sqrtQ;
    };
  } else {
    return [alpha](double, double, double) -> double {
      return alpha;
    };
  }
}

//_________________________________________________________________________________________________
ROOT::Fit::FitResult Fit(const std::vector<Digit>& digits, const std::array<double, 6> param,
                         const std::array<int, 6> fix, bool fitAsymm, std::string errorMode, double alpha = 0.)
{
  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());
  double pitch = (digits[0].getDetID() < 300) ? 0.21 : 0.25; // set pitch value

  // store digit data (x, y, dx, dy, ADC, nSamples) per cathode
  std::vector<std::array<double, 6>> digitData[2] = {{}, {}};
  for (const auto& digit : digits) {
    auto id = digit.getPadID();
    int iCath = segmentation.isBendingPad(id) ? 0 : 1;
    digitData[iCath].push_back({segmentation.padPositionX(id), segmentation.padPositionY(id),
                                segmentation.padSizeX(id) / 2., segmentation.padSizeY(id) / 2.,
                                static_cast<double>(digit.getADC()), static_cast<double>(digit.getNofSamples())});
  }

  int iCharge[2] = {4, 5}; // charge parameter index for bending and non-bending
  int nParam = 6;          // number of parameters to fit
  double param0[6];        // starting parameters
  for (int i = 0; i < 6; ++i) {
    param0[i] = param[i];
  }
  if (!fitAsymm) {
    iCharge[1] = 4;
    nParam = 5;
    param0[4] = std::sqrt(param[4] * param[5]);
  }

  // error function for chi2 calculation
  auto error = Error(errorMode, alpha);

  // mathieson functions for chi2 calculation
  auto mathiesonIntegrateX = MathiesonIntegrate(fix[2] == 1, param0[2]);
  auto mathiesonIntegrateY = MathiesonIntegrate(fix[3] == 1, param0[3]);

  // chi2 function to be minimized
  auto chi2Function = [&](const double* par) {
    double chi2 = 0.;
    for (int iCath = 0; iCath < 2; ++iCath) {
      for (const auto& data : digitData[iCath]) {
        double xmin = (data[0] - par[0] - data[2]) / pitch;
        double xmax = (data[0] - par[0] + data[2]) / pitch;
        double ymin = (data[1] - par[1] - data[3]) / pitch;
        double ymax = (data[1] - par[1] + data[3]) / pitch;
        double mgIntegrateX = mathiesonIntegrateX(xmin, xmax, par[2]);
        double mgIntegrateY = mathiesonIntegrateY(ymin, ymax, par[3]);
        double charge = par[iCharge[iCath]] * mgIntegrateX * mgIntegrateY;
        double delta = (data[4] - charge) / error(data[4], data[5], charge);
        chi2 += delta * delta;
      }
    }
    return chi2;
  };

  ROOT::Fit::Fitter fitter;
  ROOT::Math::Functor fcn(chi2Function, nParam);
  fitter.SetFCN(fcn, param0, digits.size(), 1);

  // Parameters config
  fitter.Config().ParSettings(0).SetName("x");
  fitter.Config().ParSettings(0).SetStepSize(0.1);
  fitter.Config().ParSettings(1).SetName("y");
  fitter.Config().ParSettings(1).SetStepSize(0.05);
  fitter.Config().ParSettings(2).SetName("kx");
  fitter.Config().ParSettings(2).SetLimits(0.005, 1.);
  fitter.Config().ParSettings(2).SetStepSize(0.01);
  fitter.Config().ParSettings(3).SetName("ky");
  fitter.Config().ParSettings(3).SetLimits(0.005, 1.);
  fitter.Config().ParSettings(3).SetStepSize(0.01);
  fitter.Config().ParSettings(4).SetName(fitAsymm ? "Qt_b" : "Qt");
  fitter.Config().ParSettings(4).SetStepSize(0.1);
  if (fitAsymm) {
    fitter.Config().ParSettings(5).SetName("Qt_nb");
    fitter.Config().ParSettings(5).SetStepSize(0.1);
  }

  // fix parameters if required
  for (int i = 0; i < nParam; i++) {
    if (fix[i] == 1) {
      fitter.Config().ParSettings(i).Fix();
    }
  }

  fitter.FitFCN();
  return fitter.Result();
}
