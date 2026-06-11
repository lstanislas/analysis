#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHSimulation/Response.h"
#include "Framework/Logger.h"

#include "DigitUtils.h"
#include "PreClusterUtils.h"

using o2::mch::Digit;
using o2::mch::Response;

// generate random numbers
std::mt19937 mRandom{std::random_device{}()};

//_________________________________________________________________________________________________
void ConfiguredNoise(double& charge, int iSt, bool isBending)
{

  // tuning per-station and cathode
  static constexpr std::array<std::array<double, 2>, 3> alpha = {{{1.3, 1.3}, {0.79, 0.82}, {0.77, 0.93}}};
  static constexpr std::array<std::array<double, 2>, 3> beta = {{{-0.085, -0.095}, {0.0058, 0.012}, {-0.034, -0.068}}};
  static constexpr std::array<std::array<double, 2>, 3> gamma = {{{0.0048, 0.005}, {0.0067, 0.0071}, {0.004, 0.006}}};

  static std::normal_distribution mNoise{0., 1.};
  const int cathodeIndex = isBending ? 0 : 1;

  const double a = alpha[iSt][cathodeIndex];
  const double b = beta[iSt][cathodeIndex];
  const double g = gamma[iSt][cathodeIndex];

  const double sqrtQ = std::sqrt(charge);
  const double sigma = a * sqrtQ + b * charge + g * charge * sqrtQ;
  charge += mNoise(mRandom) * sigma;
}
//_________________________________________________________________________________________________
void AddNoise(double& charge, uint32_t nSamples, std::string mode, int deId, bool isBending)
{
  // define station from deId
  int iSt = (deId < 300) ? 0 : ((deId < 500) ? 1 : 2);

  if (mode == "none") {
    return;
  }

  if (mode.rfind("predefined", 0) == 0) { // starts_with
    ConfiguredNoise(charge, iSt, isBending);
    return;
  }

  if (mode.rfind("MC_", 0) == 0) { // gaussian MC noise: 0.5*(sqrt(nSamples)+alpha)
    auto salpha = mode.substr(3);
    std::replace(salpha.begin(), salpha.end(), 'p', '.');
    double alpha = 1.0;
    try {
      alpha = std::stod(salpha);
    } catch (...) {
      LOGP(warn, "MC_ noise: invalid alpha '{}', defaulting to 1.0", salpha);
    }
    static std::normal_distribution mNoiseMC{0., 0.5};
    charge += mNoiseMC(mRandom) * (std::sqrt(static_cast<double>(nSamples)) + alpha);
    return;
  }

  if (mode.rfind("MULT_", 0) == 0) {
    auto params = mode.substr(5);
    std::replace(params.begin(), params.end(), 'p', '.');
    size_t pos1 = params.find('_');
    size_t pos2 = (pos1 == std::string::npos) ? std::string::npos : params.find('_', pos1 + 1);
    if (pos1 != std::string::npos && pos2 != std::string::npos) {
      double a = 1.0, b = 0.0, g = 0.0;
      try {
        a = std::stod(params.substr(0, pos1));
        b = std::stod(params.substr(pos1 + 1, pos2 - pos1 - 1));
        g = std::stod(params.substr(pos2 + 1));
      } catch (...) {
        LOGP(warn, "MULT_ noise: failed to parse parameters '{}', using defaults", params);
      }
      static std::normal_distribution mNoise{0., 1.};
      const double sqrtQ = std::sqrt(charge);
      const double sigma = a * sqrtQ + b * charge + g * charge * sqrtQ;
      charge += mNoise(mRandom) * sigma;
    }
    return;
  }

  LOGP(error, "unknown noise mode '{}'", mode);
}

//_________________________________________________________________________________________________
void GenerateAsymm(double& chargeB, double& chargeNB, std::string mode)
{
  double charge = std::sqrt(chargeB * chargeNB);

  if (mode == "none") {

    chargeB = charge;
    chargeNB = charge;

  } else if (mode.starts_with("gaus_")) {

    mode.erase(0, 5);
    std::replace(mode.begin(), mode.end(), 'p', '.');
    auto gamma = std::stod(mode);

    static std::normal_distribution mAsym{0., 0.055};
    double Y = gamma * mAsym(mRandom);
    chargeB = charge * exp(Y);
    chargeNB = charge * exp(-Y);

  } else if (mode == "tripleGaus") {
    // to be done
  }
}

//_________________________________________________________________________________________________
bool IsAboveThreshold(double charge, std::string mode)
{
  if (mode == "gaus") {

    static std::normal_distribution mThreshold{22.2, 2.8};
    return charge > mThreshold(mRandom);

  } else if (mode == "uniform") {

    return charge > 22.2;
  }

  return charge > 0.;
}

//_________________________________________________________________________________________________
void TMC(std::vector<Digit>& digits, int32_t time, int deId, std::array<double, 6> param,
         std::string asymm, std::string noise, std::string threshold)
{
  static const Response response[] = {{o2::mch::Station::Type1}, {o2::mch::Station::Type2345}};
  const auto& mathieson = GetMathieson((deId / 100 - 1) / 2, sqrt(param[2]), sqrt(param[3]));
  const auto& mSegmentation = o2::mch::mapping::segmentation(deId);
  int iSt = (deId < 300) ? 0 : 1;

  // generate charge asymmetry if needed
  if (asymm != "copy") {
    GenerateAsymm(param[4], param[5], asymm);
  }

  // borders of charge integration area
  auto dxy = response[iSt].getSigmaIntegration() * response[iSt].getChargeSpread();
  auto xMin = param[0] - dxy;
  auto xMax = param[0] + dxy;
  auto yMin = param[1] - dxy;
  auto yMax = param[1] + dxy;

  mSegmentation.forEachPadInArea(xMin, yMin, xMax, yMax, [&](int padid) {
    auto dx = mSegmentation.padSizeX(padid) * 0.5;
    auto dy = mSegmentation.padSizeY(padid) * 0.5;
    auto xPad = mSegmentation.padPositionX(padid) - param[0];
    auto yPad = mSegmentation.padPositionY(padid) - param[1];
    double q = mathieson.integrate(xPad - dx, yPad - dy, xPad + dx, yPad + dy);
    if (response[iSt].isAboveThreshold(q)) {
      q *= mSegmentation.isBendingPad(padid) ? param[4] : param[5];
      auto nSamples = response[iSt].nSamples(q);
      if (noise != "none") {
        AddNoise(q, nSamples, noise, deId, mSegmentation.isBendingPad(padid));
      }
      if (IsAboveThreshold(q, threshold)) {
        digits.emplace_back(deId, padid, std::round(q), time - 2, std::min(nSamples, 0x3FFU), false);
      }
    }
  });
}
