#include <cmath>
#include <optional>
#include <string>
#include <vector>
#include <fmt/format.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TROOT.h>

#include "PlotsUtils.h"
#include "ResolutionUtils.h"

//_________________________________________________________________________________________________
// require the MCH mapping to be loaded:
// gSystem->Load("libO2MCHGeometryTransformer"), gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

// This macro give the distribution of k3x and k3y for the station 1, 2 and 345
// cuts on asymmetry, total charge of the precluster, pvalue, total momentum, track angle and the distance to the closest wire can be made
// cuts mean looking into a specific range (except for wire)
void ProjectionK3Sparse(
  const std::string& inFile = "residuals_sparse.root",
  const std::string& outFile = "projection_sparse.root",
  const std::pair<std::optional<double>, std::optional<double>>& pvalue = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& asym = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& chargetot = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& p = { std::nullopt, std::nullopt },
  const std::pair<std::optional<double>, std::optional<double>>& phi = { std::nullopt, std::nullopt },
  const std::string& wire = "",
  const std::pair<std::optional<double>, std::optional<double>>& fraction = { std::nullopt, std::nullopt })
{
  auto warning = [](const std::string& label, const auto& range) {
    if (range.first && range.second) {
      std::cout << "-- WARNING -- : " << label << " selection is activated\n";
      std::cout << "SELECTION : " << *range.first << " < " << label << " < " << *range.second << "\n";
    }
  };

  warning("Asymmetry", asym);
  warning("Total Charge (ADC)", chargetot);
  warning("Track total momentum", p);
  warning("Track angle", phi);
  warning("P-Value", pvalue);
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

  // k3x and k3y distribution histograms
  std::vector<TH1D*> k3x, k3y;
  std::string sStation[3] = { "St1", "St2", "St345" };

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over the THnSparses ..." << std::endl;

  for (int i = 0; i < 3; i++) {
    auto sName = fmt::format("MultiK3PreCluster{}", sStation[i]);

    // multi dimensional histogram : contains as dim : {p-value, k3x, k3y, ADC_cluster, p , phi, nSamples, Asymm, Wire, Bending}
    auto hSparse = dynamic_cast<THnSparse*>(f.Get(sName.c_str()));
    if (!hSparse) {
      std::cerr << "Warning: Could not find THnSparse " << sName << std::endl;
      continue;
    }

    // Apply any range cuts if needed :

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
      hSparse->GetAxis(3)->SetRangeUser(*chargetot.first, *chargetot.second);
    }

    // Track total momentum (GeV/c)
    if (p.first && p.second) {
      hSparse->GetAxis(4)->SetRangeUser(*p.first, *p.second);
    }

    // Track angle (degrees)
    if (phi.first && phi.second) {
      hSparse->GetAxis(5)->SetRangeUser(*phi.first, *phi.second);
    }

    // fraction
    if (fraction.first && fraction.second) {
      hSparse->GetAxis(8)->SetRangeUser(*fraction.first, *fraction.second);
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

    // 1D projection
    TH1D* h1DX = dynamic_cast<TH1D*>(hSparse->Projection(1)); // for k3x
    TH1D* h1DY = dynamic_cast<TH1D*>(hSparse->Projection(2)); // for k3y
    auto hNameX = fmt::format("h1D_k3x_{}", sStation[i]);
    auto hTitleX = fmt::format("h1D_k3x_{}", sStation[i]);
    auto hNameY = fmt::format("h1D_k3y_{}", sStation[i]);
    auto hTitleY = fmt::format("h1D_k3y_{}", sStation[i]);
    h1DX->SetName(hNameX.c_str());
    h1DX->SetTitle(hTitleX.c_str());
    h1DY->SetName(hNameY.c_str());
    h1DY->SetTitle(hTitleY.c_str());

    k3x.push_back(h1DX);
    k3y.push_back(h1DY);
  }
  std::cout << "Saving plots ..." << std::endl;

  gStyle->SetOptStat(1);
  plot1D(k3x, "c_k3x", "K_{3X} per Station");
  plot1D(k3y, "c_k3y", "K_{3Y} per Station");

  TFile fOut(outFile.c_str(), "recreate");
  if (fOut.IsZombie()) {
    std::cerr << "Error: Cannot open output file " << outFile << std::endl;
    return;
  }

  if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject("c_k3x"))) {
    c->Write();
  } else {
    std::cerr << "Warning: Canvas c_k3x not found." << std::endl;
  }

  if (TCanvas* c = dynamic_cast<TCanvas*>(gROOT->FindObject("c_k3y"))) {
    c->Write();
  } else {
    std::cerr << "Warning: Canvas c_k3y not found." << std::endl;
  }

  fOut.Close();

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
}
