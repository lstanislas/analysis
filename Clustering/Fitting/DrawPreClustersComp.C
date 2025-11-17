#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TText.h>

#include <fmt/format.h>

#include "Framework/Logger.h"

template <class T>
T* GetObject(TFile& f, std::string objectName);
template <class T>
T* GetClone(TFile& f, std::string canvasName, std::string objectName, std::string extension);

//_________________________________________________________________________________________________
void DrawPreClustersComp(std::string file1 = "displays.root", std::string file2 = "displays.root",
                         bool norm = true)
{
  /// compare precluster and associated digits informations

  gStyle->SetOptStat(1);
  TH1::AddDirectory(false);

  TFile f[] = {{file1.c_str(), "read"}, {file2.c_str(), "read"}};
  int color[] = {2, 4};
  std::string hDrawOpt[] = {"", "sames"};

  auto l = new TPaveText(0.55, 0.75, 0.75, 0.85, "NDC");
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  for (int i = 0; i < 2; ++i) {
    auto t = l->AddText(fmt::format("file {}", i + 1).c_str());
    t->SetTextColor(color[i]);
  }

  auto hDraw = [norm, color, hDrawOpt](auto h, int i, double xMin, double xMax) {
    if (norm) {
      h->Scale(1. / h->GetEntries());
    }
    h->SetStats();
    h->SetLineColor(color[i]);
    if (xMax > xMin) {
      h->GetXaxis()->SetRangeUser(xMin, xMax);
    }
    h->Draw(hDrawOpt[i].c_str());
  };

  auto hDrawRat = [](auto h1, auto h2, double xMin, double xMax, double yMin, double yMax) {
    TH1* hRat = static_cast<TH1*>(h2->Clone());
    hRat->Divide(h1);
    hRat->SetStats(false);
    if (xMax > xMin) {
      hRat->GetXaxis()->SetRangeUser(xMin, xMax);
    }
    if (yMax > yMin) {
      hRat->GetYaxis()->SetRangeUser(yMin, yMax);
    }
    hRat->Draw();
  };

  std::string sStation[] = {"St1", "St2", "St345"};
  for (auto sSt : sStation) {

    auto cName = fmt::format("preClusterInfo{}", sSt);
    auto c = new TCanvas(fmt::format("{}Comp", cName).c_str(),
                         fmt::format("precluster characteristics {}", sSt).c_str(), 10, 10, 1200, 600);
    c->Divide(3, 2);
    auto cRat = new TCanvas(fmt::format("{}Rat", cName).c_str(),
                            fmt::format("precluster characteristics ratio {}", sSt).c_str(), 10, 10, 1200, 600);
    cRat->Divide(3, 2);
    std::string hNameBase[] = {"hCharge", "hChargeB", "hChargeNB", "hChargeAsymm2", "hDimX", "hDimY"};
    double xRange[6][2] = {{0., 10000.}, {0., 10000.}, {0., 10000.}, {1., -1.}, {0., 10.}, {0., 10.}};
    double xRangeRat[6][2] = {{0., 3000.}, {0., 3000.}, {0., 3000.}, {1., -1.}, {1., 8.}, {1., 8.}};
    double yRangeRat[6][2] = {{0.7, 1.3}, {0.7, 1.3}, {0.7, 1.3}, {0.5, 1.2}, {0.6, 1.4}, {0.6, 1.4}};
    for (int i = 0; i < 6; ++i) {
      auto hName = fmt::format("{}{}", hNameBase[i], sSt);
      TH1* h[2];
      for (int j = 0; j < 2; ++j) {
        h[j] = GetClone<TH1>(f[j], cName, hName, fmt::format("_{}", j + 1));
        c->cd(i + 1);
        gPad->SetLogy();
        hDraw(h[j], j, xRange[i][0], xRange[i][1]);
      }
      cRat->cd(i + 1);
      hDrawRat(h[0], h[1], xRangeRat[i][0], xRangeRat[i][1], yRangeRat[i][0], yRangeRat[i][1]);
    }
    c->cd(1);
    l->Clone()->Draw("same");

    cName = fmt::format("digitChargeInfo{}", sSt);
    c = new TCanvas(fmt::format("{}Comp", cName).c_str(),
                    fmt::format("digit characteristics {}", sSt).c_str(), 10, 10, 1200, 600);
    c->Divide(3, 2);
    cRat = new TCanvas(fmt::format("{}Rat", cName).c_str(),
                       fmt::format("digit characteristics ratio {}", sSt).c_str(), 10, 10, 1200, 600);
    cRat->Divide(3, 2);
    std::string hNameBase2[] = {"ADC", "ADCB", "ADCNB", "Samples", "SamplesB", "SamplesNB"};
    double xRange2[6][2] = {{0., 100.}, {0., 100.}, {0., 100.}, {0., 40.}, {0., 40.}, {0., 40.}};
    double xRangeRat2[6][2] = {{0., 1000.}, {0., 1000.}, {0., 1000.}, {0., 100.}, {0., 100.}, {0., 100.}};
    double yRangeRat2[6][2] = {{0.7, 1.3}, {0.7, 1.3}, {0.7, 1.3}, {0.6, 1.4}, {0.6, 1.4}, {0.6, 1.4}};
    for (int i = 0; i < 6; ++i) {
      auto hName = fmt::format("{}{}", hNameBase2[i], sSt);
      TH1* h[2];
      for (int j = 0; j < 2; ++j) {
        h[j] = GetClone<TH1>(f[j], cName, hName, fmt::format("_{}", j + 1));
        c->cd(i + 1);
        gPad->SetLogy();
        hDraw(h[j], j, xRange2[i][0], xRange2[i][1]);
      }
      cRat->cd(i + 1);
      hDrawRat(h[0], h[1], xRangeRat2[i][0], xRangeRat2[i][1], yRangeRat2[i][0], yRangeRat2[i][1]);
    }
    c->cd(1);
    l->Clone()->Draw("same");

    static const std::vector<std::pair<double, double>> chargeLimits{
      {-999999., 999999.},
      {0., 200.},
      {200., 400.},
      {400., 600.},
      {600., 800.},
      {800., 1000.},
      {1000., 1200.},
      {1200., 1500.},
      {1500., 2100.},
      {2100., 4000.},
      {4000., 8000.},
      {8000., 999999.}};
    c = new TCanvas(fmt::format("ChargeAsymm{}Comp", sSt).c_str(),
                    fmt::format("precluster charge asymmetry {}", sSt).c_str(), 10, 10, 1200, 900);
    c->Divide((chargeLimits.size() + 2) / 3, 3);
    cRat = new TCanvas(fmt::format("ChargeAsymm{}Rat", sSt).c_str(),
                       fmt::format("precluster charge asymmetry ratio {}", sSt).c_str(), 10, 10, 1200, 900);
    cRat->Divide((chargeLimits.size() + 2) / 3, 3);
    std::vector<TH1*> hAsymm[2];
    for (int j = 0; j < 2; ++j) {
      auto hAsymm3D = GetObject<TH3>(f[j], fmt::format("hChargeAsymm3D{}", sSt).c_str());
      if (hAsymm3D == nullptr) {
        continue;
      }
      for (size_t i = 0; i < chargeLimits.size(); ++i) {
        int firstBin = hAsymm3D->GetYaxis()->FindBin(chargeLimits[i].first + 1.);
        int lastBin = hAsymm3D->GetYaxis()->FindBin(chargeLimits[i].second - 1.);
        auto chargeMin = hAsymm3D->GetYaxis()->GetBinLowEdge(firstBin);
        auto chargeMax = hAsymm3D->GetYaxis()->GetBinLowEdge(lastBin + 1);
        auto hName = fmt::format("asymm{}_{}_{}", sSt, i, j);
        auto hTitle = fmt::format("cluster charge asymmetry {} ({} <= charge < {})", sSt,
                                  std::lround(chargeMin), std::lround(chargeMax));
        auto& h = hAsymm[j].emplace_back(hAsymm3D->ProjectionZ(hName.c_str(), 0, -1, firstBin, lastBin));
        h->SetTitle(hTitle.c_str());
        c->cd(i + 1);
        gPad->SetLogy();
        hDraw(h, j, 1., -1.);
      }
    }
    for (size_t i = 0; i < chargeLimits.size(); ++i) {
      cRat->cd(i + 1);
      hDrawRat(hAsymm[0][i], hAsymm[1][i], 1., -1., 0.5, 1.5);
    }
    c->cd(1);
    l->Clone()->Draw("same");
  }

  f[0].Close();
  f[1].Close();
}

//_________________________________________________________________________________________________
template <class T>
T* GetObject(TFile& f, std::string objectName)
{
  /// return the object

  auto o = f.FindObjectAny(objectName.c_str());
  if (o == nullptr) {
    LOGP(error, "object {} not found", objectName);
  }

  return static_cast<T*>(o);
}

//_________________________________________________________________________________________________
template <class T>
T* GetClone(TFile& f, std::string canvasName, std::string objectName, std::string extension)
{
  /// return a clone of the object from the canvas

  auto o = GetObject<TCanvas>(f, canvasName)->FindObject(objectName.c_str());
  if (o == nullptr) {
    LOGP(error, "object {} in canvas {} not found", objectName, canvasName);
    exit(-1);
  }

  return static_cast<T*>(o->Clone(fmt::format("{}{}", objectName, extension).c_str()));
}
