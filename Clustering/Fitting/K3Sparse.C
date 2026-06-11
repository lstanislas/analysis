#include <cmath>
#include <numeric>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include "TROOT.h"
#include <TTreeReader.h>
#include "TTreeReaderArray.h"
#include <TTreeReaderValue.h>

#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "DigitUtils.h"
#include "PlotsUtils.h"
#include "PreClusterUtils.h"
#include "ResolutionUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;

//_________________________________________________________________________________________________
// require the MCH mapping to be loaded:
// gSystem->Load("libO2MCHGeometryTransformer"), gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

// We fill a THnSparse (8 axis) with k3x, k3y and others information of the pre-cluster, to be treated later in ProjectionK3Sparse.C
void K3Sparse(int run, const char* inFile = "clusters.root", const char* outFile = "residuals_sparse.root",
              double minADCFit = 5)
{

  /// load CCDB objects
  InitFromCCDB(run, true, true, false);

  auto adcFitThreshold = std::abs(minADCFit);
  if (minADCFit != 0.) {
    std::cout << "-- WARNING -- : Fit ADC selection is activated" << std::endl;
    std::cout << "SELECTION : Fit ADC >= " << adcFitThreshold << std::endl;
  }

  // load histograms
  LoadHist();

  //________________________________________________________________________________________________
  // load input data and loop
  //________________________________________________________________________________________________
  std::cout << "loading data ..." << std::endl;

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");

  if (!dataReader->GetTree()->FindBranch("fitParameters")) {
    LOGP(error, "unable to load branch \"fitParameters\" from {}", inFile);
    exit(-1);
  }
  TTreeReaderArray<double> fitParameters(*dataReader, "fitParameters");

  if (!dataReader->GetTree()->FindBranch("pvalue")) {
    LOGP(error, "unable to load branch \"pvalue\" from {}", inFile);
    exit(-1);
  }
  TTreeReaderValue<double> pvalue(*dataReader, "pvalue");

  if (!dataReader->GetTree()->FindBranch("chi2")) {
    LOGP(error, "unable to load branch \"chi2\" from {}", inFile);
    exit(-1);
  }
  TTreeReaderValue<double> chi2(*dataReader, "chi2");

  int nClusters = dataReader->GetEntries(false);
  int iCluster(0);

  // multi dimensional histogram : contains as dim : {p-value, k3x, k3y, ADC_cluster, p , phi, Asymm, Wire, Bending}
  THnSparseD* hPreClusterInfoMULTIK3[3];
  hPreClusterInfoMULTIK3[0] = CreatePreClusterInfoMULTIK3("St1");
  hPreClusterInfoMULTIK3[1] = CreatePreClusterInfoMULTIK3("St2");
  hPreClusterInfoMULTIK3[2] = CreatePreClusterInfoMULTIK3("St345");

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over data ..." << std::endl;

  // loop precluster
  int discarded_cut_k3 = 0;
  int discarded_cut_ADC = 0;
  while (dataReader->Next()) {
    if (++iCluster % 10000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters
                << "..." << std::flush;
    }
    //___________________SELECTION__________________________
    // those 2 DE have lower HV for the run 529691
    if (run == 529691 && (cluster->getDEId() == 202 || cluster->getDEId() == 300)) {
      continue;
    }

    // cut on track angle at chamber
    double Track_angle = std::abs(std::atan2(trackParam->py, -trackParam->pz)) / pi * 180.;
    double Track_momentum = sqrt(trackParam->px * trackParam->px + trackParam->py * trackParam->py + trackParam->pz * trackParam->pz);
    if (Track_angle > 10.) {
      continue;
    }
    // cut on digit time
    std::vector<Digit> selectedDigits(*digits);
    selectedDigits.erase(
      std::remove_if(selectedDigits.begin(), selectedDigits.end(), [&trackTime](const auto& digit) {
        return std::abs(digit.getTime() + 1.5 - *trackTime) > 10.;
      }),
      selectedDigits.end());
    if (selectedDigits.empty()) {
      continue;
    }

    // reject mono-cathode preclusters after digit selection
    if (IsMonoCathode(selectedDigits)) {
      continue;
    }

    // reject composite preclusters
    if (IsComposite(selectedDigits, true)) {
      continue;
    }

    // check if precluster pass the fit selection if needed (N° pads, size, ...)
    if (!IsFittable(selectedDigits)) {
      continue;
    }

    // cut on precluster charge asymmetry
    auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
    double charge = sqrt(chargeNB * chargeB);
    if (std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    // cut on K3
    if ((fitParameters[2] < 1e-5) || (fitParameters[3] < 1e-5)) {
      discarded_cut_k3++;
      continue;
    }

    // cut on ADCfit
    if (minADCFit != 0.) {
      bool skip = false;
      for (const auto& digit : selectedDigits) {
        if (GetChargeIntegral(digit, {static_cast<double*>(fitParameters.GetAddress()), fitParameters.GetSize()}) < adcFitThreshold) {
          discarded_cut_ADC++;
          skip = true;
          break; // stop checking further digits if one fails
        }
      }
      if (skip && minADCFit > 0.) { // skip clusters with at least one ADCfit < minADCFit
        continue;
      }
      if (!skip && minADCFit < 0.) { // skip clusters where all ADCfit > |minADCFit| (keep only those rejected by the positive cut)
        continue;
      }
    }

    float dx_new = DistanceToClosestWire(cluster->getDEId(), fitParameters[0]); // use local X

    std::vector<double> parameters;
    for (int i = 0; i < 6; ++i) {
      parameters.push_back(fitParameters[i]);
    }
    parameters.push_back(charge);         // parameters[6]
    parameters.push_back(chargeAsymm);    // parameters[7]
    parameters.push_back(dx_new);         // parameters[8]
    parameters.push_back(*pvalue);        // parameters[9]
    parameters.push_back(Track_angle);    // parameters[10]
    parameters.push_back(Track_momentum); // parameters[11]

    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    FillK3Info(parameters, hPreClusterInfoMULTIK3[iSt]);

    // histograms that can't be in the THnSparse
    auto [nPadsNB, nPadsB] = GetNPads(selectedDigits);
    h2chi2_ndf[iSt]->Fill((nPadsNB + nPadsB - 4), *chi2);
    hprob[iSt]->Fill(*pvalue);
    hk3x[iSt]->Fill(parameters[2]);
    hk3y[iSt]->Fill(parameters[3]);
  }

  dataFileIn->Close();

  std::cout << "Creating fit status plots ..." << std::endl;
  gStyle->SetOptStat(1);
  plot2D(h2chi2_ndf, "c_chi2_ndf", "chi2 vs ndf");
  plot1D(hprob, "c_prob", "p-value");
  plot1D(hk3x, "c_k3x", "K_{3X}");
  plot1D(hk3y, "c_k3y", "K_{3Y}");

  // output
  std::cout << "Saving plots ..." << std::endl;
  TFile fOut(outFile, "recreate");
  for (THnSparseD* const& h : hPreClusterInfoMULTIK3) {
    if (h)
      h->Write();
  }
  std::vector<std::string> canvasNames = {
    "c_chi2_ndf",
    "c_prob",
    "c_k3x",
    "c_k3y",
  };
  for (const auto& name : canvasNames) {
    if (TCanvas* c = (TCanvas*)gROOT->FindObject(name.c_str())) {
      c->Write();
    } else {
      std::cerr << "Warning: Canvas " << name << " not found." << std::endl;
    }
  }

  fOut.Close();

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
  cout << "discarded clusters (cut on K3): " << discarded_cut_k3 << " / " << nClusters << endl;
  cout << "discarded clusters (cut on ADC): " << discarded_cut_ADC << " / " << nClusters << endl;
}
