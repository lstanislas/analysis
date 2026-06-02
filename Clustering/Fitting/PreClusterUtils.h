#ifndef PRECLUSTERUTILS_H_
#define PRECLUSTERUTILS_H_

#include <cmath>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <gsl/span>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include "DataFormatsMCH/Digit.h"
#include "Framework/Logger.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "MCHPreClustering/PreClusterFinderParam.h"
#include "MCHSimulation/Response.h"

#include "DigitUtils.h"

using o2::mch::Digit;

/*
 * This file contains utility functions to determine the
 * characteristics of a precluster and produce control plots.
 * To be used, they require the MCH mapping to be loaded:
 * gSystem->Load("libO2MCHMappingImpl4")
 */

//_________________________________________________________________________________________________
bool IsMonoCathode(const gsl::span<const Digit> digits)
{
  /// return true if all the digits are on the same cathode

  if (digits.empty()) {
    LOGP(warning, "IsMonoCathode: list of digits is empty");
    return true;
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  bool hasBendingPad(false);
  bool hasNonBendingPad(false);

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      if (hasNonBendingPad) {
        return false;
      }
      hasBendingPad = true;
    } else {
      if (hasBendingPad) {
        return false;
      }
      hasNonBendingPad = true;
    }
  }

  return true;
}

//_________________________________________________________________________________________________
bool IsComposite(const gsl::span<const Digit> digits, bool excludeCorners = false)
{
  /// return true if the precluster is made of disjoint groups of pads on either cathode
  /// if excludeCorners, pads touching only at one corner are considered to be disjoint

  if (digits.empty()) {
    LOGP(warning, "IsComposite: list of digits is empty");
    return false;
  }

  static o2::mch::PreClusterFinder preClusterFinder{};
  static bool init = false;
  if (o2::mch::PreClusterFinderParam::Instance().excludeCorners != excludeCorners) {
    o2::conf::ConfigurableParam::setValue("MCHPreClustering.excludeCorners", excludeCorners ? "1" : "0");
    if (init) {
      preClusterFinder.deinit();
      init = false;
    }
  }
  if (!init) {
    preClusterFinder.init();
    init = true;
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());
  std::vector<Digit> cathodeDigits[2] = {{}, {}};
  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      cathodeDigits[1].emplace_back(digit);
    } else {
      cathodeDigits[0].emplace_back(digit);
    }
  }

  for (int iCath = 0; iCath < 2; ++iCath) {
    preClusterFinder.reset();
    preClusterFinder.loadDigits(cathodeDigits[iCath]);
    int nPreClusters = preClusterFinder.run();
    if (nPreClusters > 1) {
      return true;
    }
  }

  return false;
}

//_________________________________________________________________________________________________
std::pair<int, int> GetNPads(const gsl::span<const Digit> digits)
{
  /// return the number of pads on both cathodes

  if (digits.empty()) {
    LOGP(warning, "GetNPads: list of digits is empty");
    return std::make_pair(0, 0);
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  int nPadsNB = 0;
  int nPadsB = 0;

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      ++nPadsB;
    } else {
      ++nPadsNB;
    }
  }

  return std::make_pair(nPadsNB, nPadsB);
}

//_________________________________________________________________________________________________
std::pair<int, int> GetSize(const gsl::span<const Digit> digits)
{
  /// return the size of the precluster in pad unit
  /// the size in x (y) direction is given by the number of lines of pads
  /// on the non-bending (bending) plane in that direction
  /// note: the pad size is constant in the x (y) direction on the non-bending (bending) plane

  if (digits.empty()) {
    LOGP(warning, "GetSize: list of digits is empty");
    return std::make_pair(0, 0);
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double dimX[2] = {1.e6, -1.e6};
  double padSizeX = -1.;
  double dimY[2] = {1.e6, -1.e6};
  double padSizeY = -1.;

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      double padY = segmentation.padPositionY(digit.getPadID());
      dimY[0] = std::min(dimY[0], padY);
      dimY[1] = std::max(dimY[1], padY);
      padSizeY = segmentation.padSizeY(digit.getPadID());
    } else {
      double padX = segmentation.padPositionX(digit.getPadID());
      dimX[0] = std::min(dimX[0], padX);
      dimX[1] = std::max(dimX[1], padX);
      padSizeX = segmentation.padSizeX(digit.getPadID());
    }
  }

  int sizeX = (padSizeX > 0.) ? std::lround((dimX[1] - dimX[0]) / padSizeX) + 1 : 0;
  int sizeY = (padSizeY > 0.) ? std::lround((dimY[1] - dimY[0]) / padSizeY) + 1 : 0;

  return std::make_pair(sizeX, sizeY);
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCharge(const gsl::span<const Digit> digits, bool run2 = false)
{
  /// return the total charge of the digits on both cathodes

  if (digits.empty()) {
    LOGP(warning, "GetCharge: list of digits is empty");
    return std::make_pair(0., 0.);
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  std::pair<double, double> charge{0., 0.};

  for (const auto& digit : digits) {
    if (segmentation.isBendingPad(digit.getPadID())) {
      charge.second += run2 ? adcToCharge(digit.getADC()) : digit.getADC();
    } else {
      charge.first += run2 ? adcToCharge(digit.getADC()) : digit.getADC();
    }
  }

  return charge;
}

//_________________________________________________________________________________________________
std::pair<double, double> GetChargeFraction(const gsl::span<const Digit> digits, double localX, double localY)
{
  /// return the total charge fraction seen by digits on both cathodes given the cluster position

  if (digits.empty()) {
    LOGP(warning, "GetChargeFraction: list of digits is empty");
    return std::make_pair(0., 0.);
  }

  static const o2::mch::Response response[] = {{o2::mch::Station::Type1}, {o2::mch::Station::Type2345}};

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());
  int iSt = (digits[0].getDetID() < 300) ? 0 : 1;

  std::pair<double, double> qCl{0., 0.};

  for (const auto& digit : digits) {
    auto padid = digit.getPadID();
    auto dx = segmentation.padSizeX(padid) * 0.5;
    auto dy = segmentation.padSizeY(padid) * 0.5;
    auto xPad = segmentation.padPositionX(padid) - localX;
    auto yPad = segmentation.padPositionY(padid) - localY;
    auto qPad = response[iSt].chargePadfraction(xPad - dx, xPad + dx, yPad - dy, yPad + dy);
    segmentation.isBendingPad(padid) ? qCl.second += qPad : qCl.first += qPad;
  }

  return qCl;
}

//_________________________________________________________________________________________________
bool IsFittable(const gsl::span<const Digit> digits)
{
  /// return true if the precluster has enough constraints to be fitted
  /// with 4 free parameters (2 in both X and Y directions)

  if (digits.empty()) {
    LOGP(warning, "IsFittable: list of digits is empty");
    std::cout << std::flush;
    return false;
  }

  auto [sizeX, sizeY] = GetSize(digits);
  auto [nPadsNB, nPadsB] = GetNPads(digits);

  // >= 3 pads along X on non-bending plane
  bool conditionX1 = sizeX > 2;

  // 2 pads along X on non-bending plane and at least 1 extra pad along X on bending plane
  bool conditionX2 = sizeX == 2 && sizeY != nPadsB;

  // >= 3 pads along Y on bending plane
  bool conditionY1 = sizeY > 2;

  // 2 pads along Y on bending plane and at least 1 extra pad along Y on non-bending plane
  bool conditionY2 = sizeY == 2 && sizeX != nPadsNB;

  return (conditionX1 || conditionX2) && (conditionY1 || conditionY2);
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCOG(const gsl::span<const Digit> digits)
{
  /// return the center of gravity of the digits
  /// the weight of each digit is given by its ADC charge
  /// for bi-cathode clusters, x (y) position is given by digits in the non-bending (bending) plane
  /// note: the pad size is constant in the x (y) direction on the non-bending (bending) plane

  if (digits.empty()) {
    LOGP(error, "GetCOG: list of digits is empty");
    exit(-1);
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double x[2] = {0., 0.};
  double y[2] = {0., 0.};
  double q[2] = {0., 0.};

  for (const auto& digit : digits) {
    int i = segmentation.isBendingPad(digit.getPadID()) ? 1 : 0;
    x[i] += digit.getADC() * segmentation.padPositionX(digit.getPadID());
    y[i] += digit.getADC() * segmentation.padPositionY(digit.getPadID());
    q[i] += digit.getADC();
  }

  double xCl = (q[0] > 0.) ? x[0] / q[0] : x[1] / q[1];
  double yCl = (q[1] > 0.) ? y[1] / q[1] : y[0] / q[0];

  return std::make_pair(xCl, yCl);
}

//_________________________________________________________________________________________________
std::pair<double, double> GetCOG2(const gsl::span<const Digit> digits)
{
  /// return the center of gravity of the digits
  /// the weight of each digit is given by its ADC charge / (its size / sqrt(12))^2
  /// note: the pad size is constant in the x (y) direction on the non-bending (bending) plane
  /// the purpose of weighting by the pad size is to combine the digits from both planes

  if (digits.empty()) {
    LOGP(error, "GetCOG2: list of digits is empty");
    exit(-1);
  }

  const auto& segmentation = o2::mch::mapping::segmentation(digits[0].getDetID());

  double x(0.);
  double y(0.);
  double wx(0.);
  double wy(0.);

  for (const auto& digit : digits) {
    double wxi = digit.getADC() / std::pow(segmentation.padSizeX(digit.getPadID()), 2.);
    double wyi = digit.getADC() / std::pow(segmentation.padSizeY(digit.getPadID()), 2.);
    x += wxi * segmentation.padPositionX(digit.getPadID());
    y += wyi * segmentation.padPositionY(digit.getPadID());
    wx += wxi;
    wy += wyi;
  }

  return std::make_pair(x / wx, y / wy);
}

//_________________________________________________________________________________________________
void CreatePreClusterInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create histograms of precluster characteristics

  histos.emplace_back(new TH1F(Form("hCharge%s", extension),
                               "cluster charge;charge (ADC)", 5000, -0.25, 49999.75));
  histos.emplace_back(new TH1F(Form("hChargeB%s", extension),
                               "cluster charge bending;charge (ADC)", 5000, -0.25, 49999.75));
  histos.emplace_back(new TH1F(Form("hChargeNB%s", extension),
                               "cluster charge non bending;charge (ADC)", 5000, -0.25, 49999.75));
  histos.emplace_back(new TH1F(Form("hChargeAsymm%s", extension),
                               "cluster charge asymmetry;(NB-B)/(NB+B)", 201, -1.005, 1.005));
  histos.emplace_back(new TH2F(Form("hChargeAsymm2D%s", extension),
                               "cluster charge asymmetry vs charge;(NB+B)/2 (ADC);(NB-B)/(NB+B)",
                               200, -0.25, 19999.75, 201, -1.005, 1.005));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);
  histos.emplace_back(new TH2F(Form("hChargeCorr%s", extension),
                               "cluster charge bending vs non-bending;charge NB (ADC);charge B (ADC)",
                               200, -0.25, 19999.75, 200, -0.25, 19999.75));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);
  histos.emplace_back(new TH1F(Form("hChargeAsymm2%s", extension),
                               "cluster charge asymmetry2;1/2 * ln(NB/B)", 201, -1.005, 1.005));
  histos.emplace_back(new TH2F(Form("hChargeAsymm22D%s", extension),
                               "cluster charge asymmetry2 vs charge2;#sqrt{NB*B} (ADC);1/2 * ln(NB/B)",
                               200, -0.25, 19999.75, 201, -1.005, 1.005));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);
  histos.emplace_back(new TH1F(Form("hDimX%s", extension),
                               "cluster size X;size X (#pads)", 21, -0.5, 20.5));
  histos.emplace_back(new TH1F(Form("hDimY%s", extension),
                               "cluster size Y;size Y (#pads)", 31, -0.5, 30.5));
  histos.emplace_back(new TH2F(Form("hDimXY%s", extension),
                               "cluster size Y vs X;size X (#pads);size Y (#pads)",
                               21, -0.5, 20.5, 31, -0.5, 30.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillPreClusterInfo(double chargeNB, double chargeB, int sizeX, int sizeY, std::vector<TH1*>& histos)
{
  /// fill histograms of precluster characteristics

  double charge = 0.5 * (chargeNB + chargeB);
  double charge2 = std::sqrt(chargeB * chargeNB);
  double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
  double chargeRatio = 0.5 * std::log(chargeNB / chargeB);

  histos[0]->Fill(charge);
  histos[1]->Fill(chargeB);
  histos[2]->Fill(chargeNB);
  histos[3]->Fill(chargeAsymm);
  histos[4]->Fill(charge, chargeAsymm);
  histos[5]->Fill(chargeNB, chargeB);
  histos[6]->Fill(chargeRatio);
  histos[7]->Fill(charge2, chargeRatio);
  histos[8]->Fill(sizeX);
  histos[9]->Fill(sizeY);
  histos[10]->Fill(sizeX, sizeY);
}

//_________________________________________________________________________________________________
TCanvas* DrawPreClusterInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// draw histograms of precluster characteristics

  static int logy[] = {1, 1, 0, 0, 1, 0, 1, 1, 0};
  static int logz[] = {0, 0, 1, 1, 0, 1, 0, 0, 1};
  static const char* opt[] = {"", "", "colz", "colz", "", "colz", "", "", "colz"};

  TCanvas* c = new TCanvas(Form("preClusterInfo%s", extension),
                           Form("precluster characteristics %s", extension), 10, 10, 900, 900);
  c->Divide(3, 3);

  TLegend* l = new TLegend(0.45, 0.7, 0.85, 0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  c->cd(1);
  gPad->SetLogy();
  histos[0]->SetLineColor(1);
  histos[0]->Draw();
  l->AddEntry(histos[0], "(NB+B)/2", "l");
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(2);
  histos[1]->Draw("sames");
  l->AddEntry(histos[1], "Bending", "l");
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(4);
  histos[2]->Draw("sames");
  l->AddEntry(histos[2], "Non Bending", "l");
  l->Draw("same");

  int i(0);
  for (int i = 1; i < 9; ++i) {
    c->cd(i + 1);
    gPad->SetLogy(logy[i]);
    gPad->SetLogz(logz[i]);
    histos[i + 2]->Draw(opt[i]);
  }

  return c;
}

//_________________________________________________________________________________________________
static const std::vector<double> chargeLimits{200., 400., 700., 1200., 2200., 4000.};

//_________________________________________________________________________________________________
void CreatePreClusterInfoVsWire(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create histograms of precluster characteristics versus distance to closest wire

  histos.emplace_back(new TH2F(Form("hChargeAsymmVsWire%s", extension),
                               "cluster charge asymmetry vs distance to wire;#Delta_x (cm);(NB-B)/(NB+B)",
                               25, -0.125, 0.125, 101, -1.01, 1.01));
  histos.back()->GetYaxis()->SetTitleOffset(1.5);

  for (size_t i = 0; i <= chargeLimits.size(); ++i) {
    auto chargeRange = (i == 0)                    ? fmt::format("charge < {}", chargeLimits[i])
                       : (i < chargeLimits.size()) ? fmt::format("{} <= charge < {}", chargeLimits[i - 1], chargeLimits[i])
                                                   : fmt::format("charge >= {}", chargeLimits[i - 1]);
    histos.emplace_back(new TH2F(Form("hChargeAsymmVsWire%zu%s", i, extension),
                                 Form("cluster charge asymmetry vs distance to wire (%s);#Delta_x (cm);(NB-B)/(NB+B)", chargeRange.c_str()),
                                 25, -0.125, 0.125, 101, -1.01, 1.01));
    histos.back()->GetYaxis()->SetTitleOffset(1.5);
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillPreClusterInfoVsWire(double chargeNB, double chargeB, float dx, std::vector<TH1*>& histos)
{
  /// fill histograms of precluster characteristics versus distance to closest wire

  double charge = 0.5 * (chargeNB + chargeB);
  double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
  auto i = std::distance(chargeLimits.begin(), std::upper_bound(chargeLimits.begin(), chargeLimits.end(), charge));

  histos[0]->Fill(dx, chargeAsymm);
  histos[i + 1]->Fill(dx, chargeAsymm);
}

//_________________________________________________________________________________________________
TCanvas* DrawPreClusterInfoVsWire(std::vector<TH1*>& histos, const char* extension = "")
{
  /// draw histograms of precluster characteristics versus distance to closest wire

  TCanvas* c = new TCanvas(Form("preClusterInfoVsWire%s", extension),
                           Form("precluster characteristics versus distance to wire%s", extension),
                           10, 10, 1200, 600);
  c->Divide((chargeLimits.size() + 3) / 2, 2);

  int i = 0;
  for (auto h : histos) {
    c->cd(++i);
    gPad->SetLogz();
    h->Draw("colz");
  }

  return c;
}

//_________________________________________________________________________________________________
TH3* CreatePreClusterInfo3D(const char* extension = "")
{
  /// create 3D histogram of precluster characteristics

  auto h = new TH3F(Form("hChargeAsymm3D%s", extension),
                    "cluster charge asymmetry vs charge vs distance to wire;#Delta_x (cm);#sqrt{NB*B} (ADC);1/2 * ln(NB/B)",
                    25, -0.125, 0.125, 200, -0.25, 19999.75, 201, -1.005, 1.005);
  h->GetZaxis()->SetTitleOffset(1.5);
  h->SetDirectory(0);

  return h;
}

//_________________________________________________________________________________________________
void FillPreClusterInfo3D(double chargeNB, double chargeB, float dx, TH3* h)
{
  /// fill 3D histogram of precluster characteristics

  double charge = std::sqrt(chargeB * chargeNB);
  double chargeAsymm = 0.5 * std::log(chargeNB / chargeB);

  h->Fill(dx, charge, chargeAsymm);
}

#endif // PRECLUSTERUTILS_H_
