#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <gsl/span>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TParameter.h>

std::pair<int, int> LoadHistos(const char* fileName, std::vector<TH1*> histosAtVertex[2], TH1*& hmatchChi2,
                               TH1*& hNClustersPerCh, std::vector<TH1*>& chargeHistos, std::vector<TH1*>& multHistos);
void LoadHistosAtVertex(TFile* f, std::vector<TH1*>& histos, const char* extension);
void CompareHistosAtVertex(std::vector<TH1*> histos1, int nTF1, std::vector<TH1*> histos2, int nTF2, const char* extension);
void CompareMatchChi2(TH1* h1, int nTF1, TH1* h2, int nTF2);
void CompareNClustersPerCh(TH1* h1, TH1* h2);
void LoadChargeHistos(TFile* f, std::vector<TH1*>& histos, const char* extension);
void CompareChargeHistos(gsl::span<TH1*> histos1, int nTF1, gsl::span<TH1*> histos2, int nTF2, const char* extension);
void LoadMultHistos(TFile* f, std::vector<TH1*>& histos);
void CompareMultHistos(std::vector<TH1*>& histos1, int nROF1, int nTF1, std::vector<TH1*>& histos2, int nROF2, int nTF2, bool perROF);

//_________________________________________________________________________________________________
void DisplayTracksComparison(std::string inFileName1, std::string inFileName2, bool multPerROF = false)
{
  /// Compare histograms between the 2 input files

  std::vector<TH1*> histosAtVertex1[2] = {{}, {}};
  TH1* hmatchChi21 = nullptr;
  TH1* hNClustersPerCh1 = nullptr;
  std::vector<TH1*> chargeHistos1{};
  std::vector<TH1*> multHistos1{};
  auto [nTF1, nROF1] = LoadHistos(inFileName1.c_str(), histosAtVertex1, hmatchChi21, hNClustersPerCh1, chargeHistos1, multHistos1);

  std::vector<TH1*> histosAtVertex2[2] = {{}, {}};
  TH1* hmatchChi22 = nullptr;
  TH1* hNClustersPerCh2 = nullptr;
  std::vector<TH1*> chargeHistos2{};
  std::vector<TH1*> multHistos2{};
  auto [nTF2, nROF2] = LoadHistos(inFileName2.c_str(), histosAtVertex2, hmatchChi22, hNClustersPerCh2, chargeHistos2, multHistos2);

  // display histograms
  CompareHistosAtVertex(histosAtVertex1[0], nTF1, histosAtVertex2[0], nTF2, "mch");
  CompareHistosAtVertex(histosAtVertex1[1], nTF1, histosAtVertex2[1], nTF2, "muon");
  CompareMatchChi2(hmatchChi21, nTF1, hmatchChi22, nTF2);
  CompareNClustersPerCh(hNClustersPerCh1, hNClustersPerCh2);
  CompareChargeHistos({&chargeHistos1[0], 2}, nTF1, {&chargeHistos2[0], 2}, nTF2, "AllDigits_mch");
  CompareChargeHistos({&chargeHistos1[3], 2}, nTF1, {&chargeHistos2[3], 2}, nTF2, "AllDigits_muon");
  CompareChargeHistos({&chargeHistos1[6], 2}, nTF1, {&chargeHistos2[6], 2}, nTF2, "DigitsAtClusterPos_mch");
  CompareChargeHistos({&chargeHistos1[9], 2}, nTF1, {&chargeHistos2[9], 2}, nTF2, "DigitsAtClusterPos_muon");
  CompareMultHistos(multHistos1, nROF1, nTF1, multHistos2, nROF2, nTF2, multPerROF);
}

//_________________________________________________________________________________________________
std::pair<int, int> LoadHistos(const char* fileName, std::vector<TH1*> histosAtVertex[2], TH1*& hmatchChi2,
                               TH1*& hNClustersPerCh, std::vector<TH1*>& chargeHistos, std::vector<TH1*>& multHistos)
{
  /// load histograms from input file
  /// return the number of TF analyzed to produce them

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    std::cout << "opening file " << fileName << " failed" << std::endl;
    exit(-1);
  }

  auto p1 = static_cast<TParameter<int>*>(f->FindObjectAny("nTF"));
  if (!p1) {
    std::cout << "missing number of TF" << std::endl;
    exit(1);
  }

  auto p2 = static_cast<TParameter<int>*>(f->FindObjectAny("nROFSelected"));
  if (!p2) {
    std::cout << "missing number of selected ROF" << std::endl;
  }

  LoadHistosAtVertex(f, histosAtVertex[0], "mch");
  LoadHistosAtVertex(f, histosAtVertex[1], "muon");
  hmatchChi2 = static_cast<TH1*>(f->FindObjectAny("matchChi2"));
  hNClustersPerCh = static_cast<TH1*>(f->FindObjectAny("nClustersPerCh"));
  LoadChargeHistos(f, chargeHistos, "AllDig");
  LoadChargeHistos(f, chargeHistos, "AllDigMatch");
  LoadChargeHistos(f, chargeHistos, "");
  LoadChargeHistos(f, chargeHistos, "Match");
  LoadMultHistos(f, multHistos);

  return std::make_pair(p1->GetVal(), p2 ? p2->GetVal() : 1);
}

//_________________________________________________________________________________________________
void LoadHistosAtVertex(TFile* f, std::vector<TH1*>& histos, const char* extension)
{
  /// load single muon histograms at vertex

  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("pT%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("eta%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("phi%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("dca%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("pDCA23%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("pDCA310%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("rAbs%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("nClusters%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("chi2%s", extension))));

  for (auto h : histos) {
    if (!h) {
      std::cout << "missing histo at vertex" << std::endl;
      exit(1);
    }
  }
}

//_________________________________________________________________________________________________
void CompareHistosAtVertex(std::vector<TH1*> histos1, int nTF1, std::vector<TH1*> histos2, int nTF2, const char* extension)
{
  /// draw histograms at vertex

  double nTFref = (nTF1 > nTF2) ? nTF1 : nTF2;

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos1.size() > nPadsx * nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas* cHist = new TCanvas(Form("histos_%s", extension), Form("histos_%s", extension),
                               10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    gPad->SetLogy();
    histos1[i]->SetStats(false);
    histos1[i]->Scale(nTFref / nTF1);
    // histos1[i]->Scale(1. / histos1[i]->GetEntries());
    histos1[i]->SetLineColor(4);
    histos1[i]->Draw("hist");
    histos2[i]->Scale(nTFref / nTF2);
    // histos2[i]->Scale(1. / histos2[i]->GetEntries());
    histos2[i]->SetLineColor(2);
    histos2[i]->Draw("histsame");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], Form("file 1: %g tracks", histos1[0]->GetEntries()), "l");
  lHist->AddEntry(histos2[0], Form("file 2: %g tracks", histos2[0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");

  // draw ratios
  TCanvas* cRat = new TCanvas(Form("ratios_%s", extension), Form("histos2 / histos1 - %s", extension),
                              10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cRat->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cRat->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    TH1F* hRat = static_cast<TH1F*>(histos2[i]->Clone());
    hRat->Divide(histos1[i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  }

  // draw differences
  TCanvas* cDiff = new TCanvas(Form("differences_%s", extension), Form("histos2 - histos1 - %s", extension),
                               10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cDiff->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos1.size(); ++i) {
    cDiff->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    TH1F* hDiff = static_cast<TH1F*>(histos2[i]->Clone());
    hDiff->Add(histos1[i], -1.);
    hDiff->SetStats(false);
    hDiff->SetLineColor(2);
    hDiff->Draw("hist");
  }
}

//_________________________________________________________________________________________________
void CompareMatchChi2(TH1* h1, int nTF1, TH1* h2, int nTF2)
{
  /// compare matching chi2

  double nTFref = (nTF1 > nTF2) ? nTF1 : nTF2;

  TCanvas* cHist = new TCanvas("cMatch", "cMatch", 10, 10, 900, 300);
  cHist->Divide(3, 1);
  cHist->cd(1);
  gPad->SetLogy();
  h1->SetStats(false);
  h1->Scale(nTFref / nTF1);
  // h1->Scale(1. / h1->GetEntries());
  h1->SetLineColor(4);
  h1->Draw("hist");
  h2->Scale(nTFref / nTF2);
  // h2->Scale(1. / h2->GetEntries());
  h2->SetLineColor(2);
  h2->Draw("histsame");
  cHist->cd(2);
  TH1F* hRat = static_cast<TH1F*>(h2->Clone());
  hRat->SetTitle("h2 / h1");
  hRat->Divide(h1);
  hRat->SetStats(false);
  hRat->SetLineColor(2);
  hRat->Draw();
  cHist->cd(3);
  TH1F* hDiff = static_cast<TH1F*>(h2->Clone());
  hDiff->SetTitle("h2 - h1");
  hDiff->Add(h1, -1.);
  hDiff->SetStats(false);
  hDiff->SetLineColor(2);
  hDiff->Draw("hist");

  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(h1, Form("file 1: %d TF", nTF1), "l");
  lHist->AddEntry(h2, Form("file 2: %d TF", nTF2), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CompareNClustersPerCh(TH1* h1, TH1* h2)
{
  /// compare the number of clusters per muon per chamber

  if (!h1 || !h2) {
    return;
  }

  TCanvas* cHist = new TCanvas("cNClustersPerMuonPerCh", "cNClustersPerMuonPerCh", 10, 10, 600, 300);
  cHist->Divide(2, 1);
  cHist->cd(1);
  h1->SetStats(false);
  h1->SetLineColor(4);
  h1->Draw();
  h2->SetLineColor(2);
  h2->Draw("same");
  cHist->cd(2);
  TH1F* hRat = static_cast<TH1F*>(h2->Clone());
  hRat->SetTitle("h2 / h1");
  hRat->Divide(h1);
  hRat->SetStats(false);
  hRat->SetLineColor(2);
  hRat->Draw();

  TLegend* lHist = new TLegend(0.2, 0.7, 0.4, 0.85);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(h1, "file 1", "l");
  lHist->AddEntry(h2, "file 2", "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void LoadChargeHistos(TFile* f, std::vector<TH1*>& histos, const char* extension)
{
  /// load track charge histograms

  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("ADC%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("Samples%s", extension))));
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(Form("ADCvsSample%s", extension))));

  for (auto h : histos) {
    if (!h) {
      std::cout << "missing charge histo" << std::endl;
      exit(1);
    }
  }
}

//_________________________________________________________________________________________________
void CompareChargeHistos(gsl::span<TH1*> histos1, int nTF1, gsl::span<TH1*> histos2, int nTF2, const char* extension)
{
  /// draw track charge histograms

  double nTFref = (nTF1 > nTF2) ? nTF1 : nTF2;

  TCanvas* cHist = new TCanvas(Form("cCharge%s", extension), Form("cCharge%s", extension), 10, 10, 1200, 800);
  cHist->Divide(3, 2);
  for (int i = 0; i < 2; ++i) {
    cHist->cd(3 * i + 1);
    gPad->SetLogy();
    histos1[i]->SetStats(false);
    histos1[i]->Scale(nTFref / nTF1);
    // histos1[i]->Scale(1. / histos1[i]->GetEntries());
    histos1[i]->SetLineColor(4);
    histos1[i]->Draw("hist");
    histos2[i]->Scale(nTFref / nTF2);
    // histos2[i]->Scale(1. / histos2[i]->GetEntries());
    histos2[i]->SetLineColor(2);
    histos2[i]->Draw("histsame");
    cHist->cd(3 * i + 2);
    TH1F* hRat = static_cast<TH1F*>(histos2[i]->Clone());
    hRat->SetTitle("h2 / h1");
    hRat->Divide(histos1[i]);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw("hist");
    cHist->cd(3 * i + 3);
    TH1F* hDiff = static_cast<TH1F*>(histos2[i]->Clone());
    hDiff->SetTitle("h2 - h1");
    hDiff->Add(histos1[i], -1.);
    hDiff->SetStats(false);
    hDiff->SetLineColor(2);
    hDiff->Draw("hist");
  }

  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], Form("file 1: %d TF", nTF1), "l");
  lHist->AddEntry(histos2[0], Form("file 2: %d TF", nTF2), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void LoadMultHistos(TFile* f, std::vector<TH1*>& histos)
{
  /// load cluster multiplicity histograms

  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny("nClustersTotPerCh")));
  for (int i = 0; i < 10; ++i) {
    histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(fmt::format("nClustersTotPerDEinCh{}", i + 1).c_str())));
  }
  histos.emplace_back(static_cast<TH1*>(f->FindObjectAny("nClustersTot")));
  for (int i = 0; i < 10; ++i) {
    histos.emplace_back(static_cast<TH1*>(f->FindObjectAny(fmt::format("nClustersTotinCh{}", i + 1).c_str())));
  }
}

//_________________________________________________________________________________________________
void CompareMultHistos(std::vector<TH1*>& histos1, int nROF1, int nTF1, std::vector<TH1*>& histos2, int nROF2, int nTF2, bool perROF)
{
  /// draw cluster multiplicity histograms

  if (!histos1[0] || !histos2[0]) {
    std::cout << "missing multiplicity histo" << std::endl;
    return;
  }

  int norm1 = perROF ? nROF1 : nTF1;
  int norm2 = perROF ? nROF2 : nTF2;

  auto drawHistos = [&norm1, &norm2](TH1* h1, TH1* h2) {
    h1->SetStats(false);
    h1->Scale(1. / norm1);
    h1->SetLineColor(4);
    h1->Draw();
    h2->SetStats(false);
    h2->Scale(1. / norm2);
    h2->SetLineColor(2);
    h2->Draw("sames");
    double h1minmax[2] = {0., 0.};
    h1->GetMinimumAndMaximum(h1minmax[0], h1minmax[1]);
    double h2minmax[2] = {0., 0.};
    h2->GetMinimumAndMaximum(h2minmax[0], h2minmax[1]);
    h1minmax[0] = std::min(h1minmax[0], h2minmax[0]);
    h1minmax[1] = std::max(h1minmax[1], h2minmax[1]);
    double h1range = h1minmax[1] - h1minmax[0];
    h1->SetMinimum(h1minmax[0] - 0.1 * h1range);
    h1->SetMaximum(h1minmax[1] + 0.1 * h1range);
  };

  auto drawRatio = [](TH1* h1, TH1* h2, std::string title) {
    TH1D* hRat = static_cast<TH1D*>(h2->Clone());
    hRat->SetTitle(title.c_str());
    hRat->Divide(h1);
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
  };

  TLegend* lHist = new TLegend(0.15, 0.7, 0.75, 0.85);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], fmt::format("file 1: {:.2f} ROFs/TF", static_cast<float>(nROF1) / nTF1).c_str(), "l");
  lHist->AddEntry(histos2[0], fmt::format("file 2: {:.2f} ROFs/TF", static_cast<float>(nROF2) / nTF2).c_str(), "l");

  TCanvas* cCh = new TCanvas("cNClustersPerChSelected", "cNClustersPerChSelected", 10, 10, 600, 300);
  cCh->Divide(2, 1);
  cCh->cd(1);
  drawHistos(histos1[0], histos2[0]);
  lHist->Clone()->Draw("same");
  cCh->cd(2);
  drawRatio(histos1[0], histos2[0], "h2 / h1");

  TCanvas* cDE = new TCanvas("cNClustersPerDESelected", "cNClustersPerDESelected", 10, 10, 1500, 600);
  cDE->Divide(5, 2);
  TCanvas* cDERat = new TCanvas("cNClustersPerDERatioSelected", "cNClustersPerDERatioSelected", 10, 10, 1500, 600);
  cDERat->Divide(5, 2);
  for (int i = 1; i <= 10; ++i) {
    cDE->cd(i);
    drawHistos(histos1[i], histos2[i]);
    cDERat->cd(i);
    drawRatio(histos1[i], histos2[i], fmt::format("h2 / h1 in chamber {}", i));
  }
  cDE->cd(1);
  lHist->Clone()->Draw("same");

  TCanvas* cN = new TCanvas("cNClustersSelected", "cNClustersSelected", 10, 10, 600, 300);
  cN->Divide(2, 1);
  cN->cd(1);
  drawHistos(histos1[11], histos2[11]);
  histos1[11]->SetMinimum();
  gPad->SetLogy();
  lHist->Clone()->Draw("same");
  cN->cd(2);
  drawRatio(histos1[11], histos2[11], "h2 / h1");

  TCanvas* cNCh = new TCanvas("cNClustersPerChDistSelected", "cNClustersPerChDistSelected", 10, 10, 1500, 600);
  cNCh->Divide(5, 2);
  TCanvas* cNChRat = new TCanvas("cNClustersPerChRatioSelected", "cNClustersPerChRatioSelected", 10, 10, 1500, 600);
  cNChRat->Divide(5, 2);
  for (int i = 1; i <= 10; ++i) {
    cNCh->cd(i);
    drawHistos(histos1[11 + i], histos2[11 + i]);
    histos1[11 + i]->SetMinimum();
    gPad->SetLogy();
    cNChRat->cd(i);
    drawRatio(histos1[11 + i], histos2[11 + i], fmt::format("h2 / h1 in chamber {}", i));
  }
  cNCh->cd(1);
  lHist->Clone()->Draw("same");
}
