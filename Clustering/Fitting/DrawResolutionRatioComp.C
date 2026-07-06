#include <cmath>
#include <vector>
#include <sstream>
#include <string>

#include <fmt/format.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TKey.h>
#include <TH2D.h>
#include <TH1D.h>

#include "PlotsUtils.h"

//_________________________________________________________________________________________________
// This macro take the TList generated from ProjectionSparse.C for TMC and DATA
// and return 5 types of canvas for the different stations (4x3 + 3x2 in total) ->
// 1st type : plots of the std of the residuals distribution vs ADCfit
// 2nd type : plots of the ratio between TMC and DATA
// 3rd type : plots of the mean (the one from the fit of the std extraction) vs ADCfit for DATA and TMC
// 4th type : plots of the reduced chi2 (the one from the fit of the std extraction) vs ADCfit for DATA and TMC
// 5th type : 3 (station) x 2 (cathode) plots of the residuals distribution for different ADCfit binning range
// sparseTypes =  "Fit" , "Noise" or "Total"
// noise : "" (default) = sADC with alpha=1 ; "sADC_XpX" = alpha*sqrt(ADC) ; "MC_XpX" = 0.5*(sqrt(nSamples)+XpX) ; "MULT_XpX_XpX_XpX" = XpX*sqrt(ADC) + XpX*ADC + XpX*ADC*sqrt(ADC)

void DrawResolutionRatioComp(const std::string& file1 = "data_projection_sparse.root",
                             const std::string& file2 = "tmc_projection_sparse.root",
                             const std::string& outFile = "resolution_ratio.root",
                             const std::string& sparseType1 = "Fit",
                             const std::string& sparseType2 = "Fit",
                             const std::string& noise = "")
{
  TH1::AddDirectory(kFALSE);

  static const std::string sStation[3] = {"St1", "St2", "St345"};
  static const std::string sCathode[2] = {"Bend", "NBend"};

  // ADCfit binning range for the 5th type
  static const std::vector<std::pair<int, int>> chargeLimits{
    {20, 39}, {40, 59}, {60, 79}, {80, 119}, {120, 199}, {200, 399}, {400, 699}, {700, 1000}};

  auto extractGraphsAndHistos = [](TFile& f, const std::string& sparseType, std::vector<TGraphAsymmErrors*>& graphs,
                                   std::vector<TGraph*>& graphs1, std::vector<TGraph*>& graphs2,
                                   std::vector<TH1D*>& histograms) {
    // get canvas with 2D plots
    auto c2DName = fmt::format("c_2Dresiduals_{}", sparseType);
    TCanvas* c2D = dynamic_cast<TCanvas*>(f.Get(c2DName.c_str()));
    if (!c2D) {
      Warning("LoadCanvas", "Canvas '%s' not found in file!", c2DName.c_str());
      return;
    }

    // loop (St.1, St.2, St.345) x (B, NB)
    for (int i = 0; i < 6; ++i) {

      // --- get 2D histos and project them in ADCfit binning ranges ---
      auto h2DName = fmt::format("h2D_{}_{}_{}", sparseType, sStation[i / 2], sCathode[i % 2]);
      TH2D* h2D = dynamic_cast<TH2D*>(c2D->FindObject(h2DName.c_str()));
      if (!h2D) {
        Warning("LoadHistos", "h2D '%s' not found in file!", h2DName.c_str());
        continue;
      }

      for (const auto& [min, max] : chargeLimits) {
        histograms.push_back(h2D->ProjectionY(fmt::format("[{},{}]", min, max).c_str(),
                                              h2D->GetXaxis()->FindBin(min), h2D->GetXaxis()->FindBin(max)));
      }

      // --- get residuals and extract fit results ---
      auto lName = fmt::format("Residual_{}_{}_{}", sparseType, sStation[i / 2], sCathode[i % 2]);
      TList* list = dynamic_cast<TList*>(f.Get(lName.c_str()));
      if (!list) {
        Warning("LoadTLists", "TList '%s' not found in file!", lName.c_str());
        continue;
      }

      std::vector<double> x, exl, exr, y, ey, mean, Rchi2;
      for (TObject* obj : *list) {
        auto h = dynamic_cast<TH1*>(obj);
        if (!h) {
          std::cerr << "Error: Histogram not found or invalid." << std::endl;
          continue;
        }
        auto fit = h->GetFunction("fit");
        if (!fit) {
          std::cerr << "Error: Fit function not found or invalid" << std::endl;
          continue;
        }

        std::string name = h->GetName();
        std::stringstream ss(name);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, '_')) { // read the name of the elements from the TList which contain
          tokens.push_back(token);             // the charge interval in the name, e.g. "projY_h2D_Label_Station_Cathode_ADCmin_ADCmax"
        }
        if (tokens.size() < 7) {
          continue;
        }

        double X = std::stod(tokens[5]);     // lower ADC range
        double Y = std::stod(tokens[6]);     // higher ADC range
        double avg = 0.5 * (X + Y);          // mean of the ADC range
        double dx = std::abs(X - avg) + 0.5; // lower bin edge length
        double dy = std::abs(Y - avg) + 0.5; // higher bin edge length

        x.push_back(avg);
        exl.push_back(dx);
        exr.push_back(dy);
        y.push_back(fit->GetParameter(2));
        ey.push_back(fit->GetParError(2) / 2);
        mean.push_back(fit->GetParameter(1));
        Rchi2.push_back((fit->GetNDF() > 0) ? (fit->GetChisquare() / fit->GetNDF()) : 0);
      }

      // vector that contains graphs with asymm errors of the std of the residuals distribution vs ADCfit
      graphs.push_back(new TGraphAsymmErrors((int)x.size(), x.data(), y.data(), exl.data(), exr.data(), ey.data(), ey.data()));
      // vector for file1
      graphs1.push_back(new TGraph((int)x.size(), x.data(), mean.data()));
      // vector for file2
      graphs2.push_back(new TGraph((int)x.size(), x.data(), Rchi2.data()));
    }
  };

  TFile f1(file1.c_str(), "read");
  TFile f2(file2.c_str(), "read");
  TFile fout(outFile.c_str(), "recreate");

  std::vector<TGraphAsymmErrors*> g1, g2;
  std::vector<TGraph*> mean1, mean2, Rchi2_1, Rchi2_2;
  std::vector<TH1D*> h1, h2;
  extractGraphsAndHistos(f1, sparseType1, g1, mean1, Rchi2_1, h1);
  extractGraphsAndHistos(f2, sparseType2, g2, mean2, Rchi2_2, h2);

  // name, color, title, range, ...
  for (int i = 0; i < 6; i++) {

    std::string cathode = (i % 2 == 0) ? "(B)" : "(NB)";

    FillInfoGraphErrAsymm(g1[i], "ADC fit", "#sigma", sStation[i / 2], cathode, true);
    FillInfoGraph(mean1[i], "ADC fit", "#mu", sStation[i / 2], cathode, true);
    FillInfoGraph(Rchi2_1[i], "ADC fit", "#chi^{2}/ndf", sStation[i / 2], cathode, true);

    FillInfoGraphErrAsymm(g2[i], "ADC fit", "#sigma", sStation[i / 2], cathode, false);
    FillInfoGraph(mean2[i], "ADC fit", "#mu", sStation[i / 2], cathode, false);
    FillInfoGraph(Rchi2_2[i], "ADC fit", "#chi^{2}/ndf", sStation[i / 2], cathode, false);

    int nChargeLimits = chargeLimits.size();
    for (int j = 0; j < nChargeLimits; j++) {
      FillInfoHist(h1[j + nChargeLimits * i], "Residuals", "Counts", sStation[i / 2], cathode, true);
      FillInfoHist(h2[j + nChargeLimits * i], "Residuals", "Counts", sStation[i / 2], cathode, false);
    }
  }

  for (int i = 0; i < 3; i++) {
    std::cout << "Creating plots for " << sStation[i] << "..." << std::endl;
    // plot Resolution vs ADC from fitted pad
    auto gName = fmt::format("g_sigmaADC_{}", sStation[i]);
    tGraphErrAsymm(g1, g2, gName, i, noise);

    // plot Ratio
    auto rName = fmt::format("r_Ratio_{}", sStation[i]);
    tRatio(g1, g2, rName, i, noise);

    // scattered plot extracted mean vs ADC from fitted pad
    auto mName = fmt::format("m_MeanADC_{}", sStation[i]);
    tGraph(mean1, mean2, mName, i);

    // scattered plot extracted reduced chi2 vs ADC from fitted pad
    auto rcName = fmt::format("rc_Rchi2ADC_{}", sStation[i]);
    tGraph(Rchi2_1, Rchi2_2, rcName, i);

    // plot histograms of the residuals of file1 and file2 for different ADC binning.
    for (int j = 0; j < 2; j++) {
      std::string cathode = (j == 0) ? "(B)" : "(NB)";
      auto hName = fmt::format("h_BinResidualADC_{}_{}", sStation[i], cathode);
      tHist(h1, h2, hName, i, j, chargeLimits.size());
    }
  }

  std::cout << "Saving plots ..." << std::endl;
  for (int i = 0; i < 3; i++) {
    auto gName = fmt::format("g_sigmaADC_{}", sStation[i]);
    if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject(gName.c_str()))) {
      c->Write();
    } else {
      std::cerr << "Warning: Resolution vs ADC Canvas not found or invalid." << std::endl;
    }

    auto rName = fmt::format("r_Ratio_{}", sStation[i]);
    if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject(rName.c_str()))) {
      c->Write();
    } else {
      std::cerr << "Warning: Ratio Canvas not found or invalid." << std::endl;
    }

    auto mName = fmt::format("m_MeanADC_{}", sStation[i]);
    if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject(mName.c_str()))) {
      c->Write();
    } else {
      std::cerr << "Warning: Scattered Mean Canvas not found or invalid." << std::endl;
    }

    auto rcName = fmt::format("rc_Rchi2ADC_{}", sStation[i]);
    if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject(rcName.c_str()))) {
      c->Write();
    } else {
      std::cerr << "Warning: Scattered RChi2 Canvas not found or invalid." << std::endl;
    }

    for (int j = 0; j < 2; j++) {
      std::string cathode = (j == 0) ? "(B)" : "(NB)";
      auto hName = fmt::format("h_BinResidualADC_{}_{}", sStation[i], cathode);
      if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject(hName.c_str()))) {
        c->Write();
      } else {
        std::cerr << "Warning: Histogram residuals for ADC binning Canvas not found or invalid." << std::endl;
      }
    }
  }

  fout.Close();
}
