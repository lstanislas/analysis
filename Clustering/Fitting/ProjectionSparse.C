#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TROOT.h>
#include <fmt/format.h>

#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "PlotsUtils.h"
#include "ResolutionUtils.h"

//_________________________________________________________________________________________________
// This macro creates a root file containing TList for each (station x cathode x sparse type) combination
// sparse types : Fit (sigma_output), Noise (sigma_noise, ToyMC only), Total (sqrt(sigma_noise^2+sigma_Y^2), ToyMC only)
// + one TH2D which represents the residuals (ADC data - ADC fit) vs ADC fit
// cuts on asymmetry, total charge of the precluster, pvalue, NofSamples and the distance to the closest wire can be made
// cuts mean looking into a specific range (except for wire)
void ProjectionSparse(
  const std::string& inFile = "residuals_sparse.root",
  const std::string& outFile = "projection_sparse.root",
  const bool auto_bin = true,                 // adjust bin in the resolution extraction
  const std::pair<int, int>& projYX = {1, 2}, // default is Residual vs ADC from fitted pad
  const std::pair<std::optional<double>, std::optional<double>>& asym = {std::nullopt, std::nullopt},
  const std::pair<std::optional<double>, std::optional<double>>& chargetot = {std::nullopt, std::nullopt},
  const std::pair<std::optional<double>, std::optional<double>>& pvalue = {std::nullopt, std::nullopt},
  const std::pair<std::optional<double>, std::optional<double>>& nsamples = {std::nullopt, std::nullopt},
  const std::pair<std::optional<double>, std::optional<double>>& fraction = {std::nullopt, std::nullopt},
  const std::string& wire = "")
{
  auto warning = [](const std::string& label, const auto& range) {
    if (range.first && range.second) {
      std::cout << "-- WARNING -- : " << label << " selection is activated\n";
      std::cout << "SELECTION : " << *range.first << " < " << label << " < " << *range.second << "\n";
    }
  };

  warning("Asymmetry", asym);
  warning("Total Charge (ADC)", chargetot);
  warning("P-Value", pvalue);
  warning("NofSamples", nsamples);
  warning("Fraction", fraction);

  if (!wire.empty()) {
    std::cout << "-- WARNING -- : WIRE selection is activated\n";
    std::cout << "SELECTION : " << wire << std::endl;
  }

  std::cout << "loading data ..." << std::endl;
  TFile f(inFile.c_str(), "read");
  if (f.IsZombie()) {
    std::cerr << "Error: Cannot open file " << inFile << std::endl;
    return;
  }

  // TList created on-the-fly for each (station x cathode x sparse type) found in the file
  std::map<std::string, std::vector<TH2D*>> labelHistos;
  std::vector<TList*> allListResidual;
  std::string sStation[3] = {"St1", "St2", "St345"};

  // sparse types : {THnSparse name prefix, output label}
  // Noise and Total are only present in files produced from ToyMC input (see ResidualsSparse.C)
  const std::vector<std::pair<std::string, std::string>> sparseTypes = {
    {"MultiResolutionPreCluster", "Fit"},        // sigma_output (DATA or fitted TMC)
    {"MultiResolutionPreClusterNoise", "Noise"}, // sigma_noise (ToyMC only, asymm="copy" or "none")
    {"MultiResolutionPreClusterTotal", "Total"}  // sqrt(sigma_noise^2 + sigma_Y^2) (ToyMC only)
  };

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over the THnSparses ..." << std::endl;

  for (const auto& [prefix, label] : sparseTypes) {
    for (int i = 0; i < 6; i++) {
      auto sName = fmt::format("{}{}", prefix, sStation[i / 2]);

      // multi dimensional histogram whose axes are : {pvalue, residuals, ADC_fit, ADC_mes, ADC_cluster, nSamples, Asymm, Wire, Cathode, fraction}
      auto hSparse = dynamic_cast<THnSparse*>(f.Get(sName.c_str()));
      if (!hSparse) {
        continue; // skip if sparse not present in file (Noise/Total absent on DATA input)
      }

      // apply range cuts if needed :

      // P-Value
      if (pvalue.first && pvalue.second) {
        hSparse->GetAxis(0)->SetRangeUser(*pvalue.first, *pvalue.second);
      }

      // Asymmetry
      if (asym.first && asym.second) {
        hSparse->GetAxis(6)->SetRangeUser(*asym.first, *asym.second);
      }

      // Total Charge (ADC)
      if (chargetot.first && chargetot.second) {
        hSparse->GetAxis(4)->SetRangeUser(*chargetot.first, *chargetot.second);
      }

      // nSamples
      if (nsamples.first && nsamples.second) {
        hSparse->GetAxis(5)->SetRangeUser(*nsamples.first, *nsamples.second);
      }

      // fraction
      if (fraction.first && fraction.second) {
        hSparse->GetAxis(9)->SetRangeUser(*fraction.first, *fraction.second);
      }

      // Wire
      TAxis* axis7 = hSparse->GetAxis(7);
      if (wire == "top") {
        axis7->SetRange(1, 1);
      } else if (wire == "between") {
        axis7->SetRange(3, 3);
      } else if (wire == "crossover") {
        axis7->SetRange(2, 2);
      }

      // Cathode
      TAxis* axis8 = hSparse->GetAxis(8);
      std::string cathode = (i % 2 == 0) ? "Bend" : "NBend";
      if (i % 2 == 0) {
        axis8->SetRange(2, 2); // Bending bin
      } else {
        axis8->SetRange(1, 1); // NonBending bin
      }

      // 2D projection
      TH2D* h2D = dynamic_cast<TH2D*>(hSparse->Projection(projYX.first, projYX.second));
      auto hName = fmt::format("h2D_{}_{}_{}", label, sStation[i / 2], cathode);
      h2D->SetName(hName.c_str());
      h2D->SetTitle(hName.c_str());

      auto lName = fmt::format("Residual_{}_{}_{}", label, sStation[i / 2], cathode);
      TList* list = new TList();
      list->SetName(lName.c_str());

      int statistic = 2000;
      Resolution(list, h2D, statistic, auto_bin);
      allListResidual.push_back(list);

      // reduce the binning to avoid memory issue when saving the canvas
      h2D->RebinX(5);

      int firstBin = h2D->GetXaxis()->FindBin(20);
      int lastBin = h2D->FindLastBinAbove(0);
      if (lastBin > firstBin) {
        h2D->GetXaxis()->SetRange(firstBin, lastBin);
      }

      labelHistos[label].push_back(h2D);
    }
  }
  gStyle->SetOptStat(1);
  for (const auto& [prefix, label] : sparseTypes) {
    if (labelHistos.count(label) && !labelHistos[label].empty()) {
      auto cName = fmt::format("c_2Dresiduals_{}", label);
      plot2D(labelHistos[label], cName.c_str(), "residuals vs ADC", true);
    }
  }

  std::cout << "Saving plots ..." << std::endl;

  TFile fOut(outFile.c_str(), "recreate");
  if (fOut.IsZombie()) {
    std::cerr << "Error: Cannot open output file " << outFile << std::endl;
    return;
  }

  for (auto* list : allListResidual) {
    fOut.WriteTObject(list, list->GetName());
  }

  for (const auto& [prefix, label] : sparseTypes) {
    auto cName = fmt::format("c_2Dresiduals_{}", label);
    if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject(cName.c_str()))) {
      c->Write();
    }
  }

  fOut.Close();

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
}
