#include <cmath>
#include <numeric>
#include <vector>

#include <fmt/format.h>

#include <TFile.h>
#include <TROOT.h>
#include <TTreeReaderArray.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <THnSparse.h>
#include <TH2D.h>
#include <TH1D.h>

#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "ResolutionUtils.h"
#include "PreClusterUtils.h"
#include "PlotsUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;
//_________________________________________________________________________________________________
// require the MCH mapping to be loaded:
// gSystem->Load("libO2MCHGeometryTransformer"),  gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

// Fill a THnSparse (10 axis) with the residuals and others information of the pre-cluster, to be treated later in ProjectionSparse.C
void ResidualsSparse(int run, const char* inFile = "clusters.root", const char* outFile = "residuals_sparse.root", int correctADCfit = 5, int correctADC = 15, bool correctCharge = false)
{

  // load CCDB objects
  InitFromCCDB(run, true, true, false);

  if (correctADCfit != 0) {
    std::cout << "-- WARNING -- : Fit ADC selection is activated " << std::endl;
    std::cout << "SELECTION : " << std::abs(correctADCfit) << " < FitADC" << std::endl;
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
  std::unique_ptr<TTreeReaderArray<double>> fitParameters{};
  std::unique_ptr<TTreeReaderArray<double>> trueParameters{};

  if (!dataReader->GetTree()->FindBranch("fitParameters")) {
    LOGP(error, "unable to load branch \"fitParameters\" from {}", inFile);
    exit(-1);
  }
  fitParameters = std::make_unique<TTreeReaderArray<double>>(*dataReader, "fitParameters");
  // detect ToyMC input : "parameters" carries true pre-noise values {X, Y, K3x, K3y, Qb_true, Qnb_true} (optional)
  if (dataReader->GetTree()->FindBranch("parameters")) {
    trueParameters = std::make_unique<TTreeReaderArray<double>>(*dataReader, "parameters");
  }
  bool isTMC = trueParameters != nullptr;

  std::unique_ptr<TTreeReaderValue<double>> pvalue{};
  if (!dataReader->GetTree()->FindBranch("pvalue")) {
    LOGP(error, "unable to load branch \"pvalue\" from {}", inFile);
    exit(-1);
  }
  pvalue = std::make_unique<TTreeReaderValue<double>>(*dataReader, "pvalue");

  std::unique_ptr<TTreeReaderValue<double>> chi2{};
  if (!dataReader->GetTree()->FindBranch("chi2")) {
    LOGP(error, "unable to load branch \"chi2\" from {}", inFile);
    exit(-1);
  }
  chi2 = std::make_unique<TTreeReaderValue<double>>(*dataReader, "chi2");

  int nClusters = dataReader->GetEntries(false);
  int iCluster(0);

  // multi dimensional histogram : contains as dim : {p-value, residuals, ADC_fit, ADC_mes, ADC_cluster, nSamples, Asymm, Wire, Bending, fraction}
  THnSparseD* hPreClusterInfoMULTI[3];
  hPreClusterInfoMULTI[0] = CreatePreClusterInfoMULTI("St1");
  hPreClusterInfoMULTI[1] = CreatePreClusterInfoMULTI("St2");
  hPreClusterInfoMULTI[2] = CreatePreClusterInfoMULTI("St345");
  // ToyMC only : sigma_noise (Qb_true/Qnb_true reference) and sigma_total = sqrt(sigma_noise^2 + sigma_Y^2) (Q_tot reference)
  THnSparseD* hPreClusterInfoNoise[3] = {nullptr, nullptr, nullptr};
  THnSparseD* hPreClusterInfoTotal[3] = {nullptr, nullptr, nullptr};
  if (isTMC) {
    hPreClusterInfoNoise[0] = CreatePreClusterInfoMULTI("NoiseSt1");
    hPreClusterInfoNoise[1] = CreatePreClusterInfoMULTI("NoiseSt2");
    hPreClusterInfoNoise[2] = CreatePreClusterInfoMULTI("NoiseSt345");
    hPreClusterInfoTotal[0] = CreatePreClusterInfoMULTI("TotalSt1");
    hPreClusterInfoTotal[1] = CreatePreClusterInfoMULTI("TotalSt2");
    hPreClusterInfoTotal[2] = CreatePreClusterInfoMULTI("TotalSt345");
  }

  static const char* sStationNames[3] = {"St1", "St2", "St345"};
  TH2D* h2ADCtrueVsADCfit[3] = {nullptr, nullptr, nullptr};
  if (isTMC) {
    for (int i = 0; i < 3; ++i) {
      auto hName = fmt::format("h2ADCtrueVsADCfit_{}", sStationNames[i]);
      h2ADCtrueVsADCfit[i] = new TH2D(hName.c_str(), sStationNames[i], 500, -0.5, 9999.5, 500, -0.5, 9999.5);
      h2ADCtrueVsADCfit[i]->SetDirectory(0);
      h2ADCtrueVsADCfit[i]->GetXaxis()->SetTitle("ADC_fit");
      h2ADCtrueVsADCfit[i]->GetYaxis()->SetTitle("ADC_true");
    }
  }

  auto tStart = std::chrono::high_resolution_clock::now();
  std::cout << "looping over data ..." << std::endl;

  // loop precluster data
  int discarded_cut_k3 = 0;
  int discarded_cut_ADC = 0;
  int discarded_cut_ADC_fit = 0;
  while (dataReader->Next()) {

    if (++iCluster % 10000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters << "..." << std::flush;
    }
    //___________________SELECTION__________________________
    // those 2 DE have lower HV for the run 529691
    if (run == 529691 && (cluster->getDEId() == 202 || cluster->getDEId() == 300)) {
      continue;
    }

    // cut on track angle at chamber
    if (std::abs(std::atan2(trackParam->py, -trackParam->pz)) / pi * 180. > 10.) {
      continue;
    }
    // cut on digit time
    std::vector<Digit> selectedDigits(*digits);
    selectedDigits.erase(
      std::remove_if(selectedDigits.begin(), selectedDigits.end(), [&trackTime](const auto& digit) { return std::abs(digit.getTime() + 1.5 - *trackTime) > 10.; }),
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

    // correct charge total bending and charge total nonbending
    auto local = GlobalToLocal(cluster->getDEId(), cluster->x, cluster->y, cluster->z, run < 300000);
    if (correctCharge) {
      auto [chargeFracNB, chargeFracB] = GetChargeFraction(selectedDigits, local.x(), local.y());
      chargeNB /= chargeFracNB;
      chargeB /= chargeFracB;
    }

    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
    double charge = sqrt(chargeNB * chargeB);

    if (std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    std::vector<double> parameters;
    for (int i = 0; i < 6; ++i) {
      parameters.push_back((*fitParameters)[i]);
    }

    // cut on K3
    if ((parameters[2] < 1e-5) || (parameters[3] < 1e-5)) {
      discarded_cut_k3++;
      continue;
    }

    // cut on ADC (protection) & ADCfit
    bool skip = false;
    bool skip2 = false;
    const uint32_t adcThreshold = static_cast<uint32_t>(std::abs(correctADC));
    for (auto digit : selectedDigits) {
      if (digit.getADC() < adcThreshold) {
        discarded_cut_ADC++;
        skip2 = true;
        break; // stop checking further digits if one fails
      }

      if (ADCFit(digit, parameters) < std::abs(correctADCfit)) {
        discarded_cut_ADC_fit++;
        skip = true;
        break; // stop checking further digits if one fails
      }
    }

    if (skip2 && correctADC > 0) { // skip clusters with at least one ADC < correctADC
      continue;
    }
    if (!skip2 && correctADC < 0) { // skip clusters where all ADC > |correctADC| (keep only those rejected by the positive cut)
      continue;
    }

    if (skip && correctADCfit > 0) { // skip clusters with at least one ADCfit < correctADCfit
      continue;
    }
    if (!skip && correctADCfit < 0) { // skip clusters where all ADCfit > |correctADCfit| (keep only those rejected by the positive cut)
      continue;
    }

    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    float dx_new = DistanceToClosestWire(cluster->getDEId(), parameters[0]); // use local X

    parameters.push_back(charge);      // parameters[6]
    parameters.push_back(chargeAsymm); // parameters[7]
    parameters.push_back(dx_new);      // parameters[8]
    parameters.push_back(**pvalue);    // parameters[9]
    parameters.push_back(chargeB);     // parameters[10]
    parameters.push_back(chargeNB);    // parameters[11]

    for (auto digit : selectedDigits) {
      FillResolutionInfo(digit, parameters, hPreClusterInfoMULTI[iSt]);
    }

    // fill true-charge THnSparses for ToyMC (Measured - True)
    // use fitted TMC position/K3/wire/pvalue for axes, override only charges with true pre-noise values
    if (isTMC) {
      double Qb_true = (*trueParameters)[4];
      double Qnb_true = (*trueParameters)[5];
      double Q_tot = std::sqrt(Qb_true * Qnb_true);

      // parNoise: residual = ADC_measured - rho(X_true, K3_true) * Q_true_cathode = pure noise
      std::vector<double> parNoise = parameters; // start from fitted TMC parameters (all 12 elements)
      parNoise[0] = (*trueParameters)[0];        // X_true
      parNoise[1] = (*trueParameters)[1];        // Y_true
      parNoise[2] = (*trueParameters)[2];        // K3x_true
      parNoise[3] = (*trueParameters)[3];        // K3y_true
      parNoise[4] = Qb_true;
      parNoise[5] = Qnb_true;
      parNoise[10] = Qb_true;
      parNoise[11] = Qnb_true;

      // parTotal: residual = ADC_measured - rho(X_true, K3_true) * Q_tot = noise + asymmetry (non-zero when Qb_true != Qnb_true)
      std::vector<double> parTotal = parameters;
      parTotal[0] = (*trueParameters)[0]; // X_true
      parTotal[1] = (*trueParameters)[1]; // Y_true
      parTotal[2] = (*trueParameters)[2]; // K3x_true
      parTotal[3] = (*trueParameters)[3]; // K3y_true
      parTotal[4] = Q_tot;
      parTotal[5] = Q_tot;
      parTotal[10] = Q_tot;
      parTotal[11] = Q_tot;

      for (auto digit : selectedDigits) {
        FillResolutionInfo(digit, parNoise, hPreClusterInfoNoise[iSt]);
        FillResolutionInfo(digit, parTotal, hPreClusterInfoTotal[iSt]);
      }
      {
        std::vector<double> pre_fit(parameters.begin(), parameters.begin() + 6);
        std::vector<double> pre_true(parNoise.begin(), parNoise.begin() + 6);
        for (auto digit : selectedDigits) {
          h2ADCtrueVsADCfit[iSt]->Fill(ADCFit(digit, pre_fit), ADCFit(digit, pre_true));
        }
      }
    }

    // histograms that can't be in the THnSparse
    hAsymm[2 * iSt]->Fill(chargeAsymm);
    hAsymm[2 * iSt + 1]->Fill((parameters[5] - parameters[4]) / (parameters[4] + parameters[5]));
    auto [nPadsNB, nPadsB] = GetNPads(selectedDigits);
    h2chi2_ndf[iSt]->Fill((nPadsNB + nPadsB - 4), **chi2);
    hprob[iSt]->Fill(**pvalue);
  }

  dataFileIn->Close();

  std::cout << "Creating fit status plots ..." << std::endl;
  gStyle->SetOptStat(1);
  plot2D(h2chi2_ndf, "c_chi2_ndf", "chi2 vs ndf");
  plot1D(hprob, "c_prob", "p-value");
  plotSAME(hAsymm, "c_asymm", "Asymmetry");
  if (isTMC) {
    TCanvas* cADC = new TCanvas("c_ADCtrueVsADCfit", "ADC_true vs ADC_fit", 1800, 600);
    cADC->Divide(3, 1);
    for (int i = 0; i < 3; ++i) {
      cADC->cd(i + 1);
      gPad->SetLogz();
      h2ADCtrueVsADCfit[i]->Draw("colz");
    }
  }

  // output
  std::cout << "Saving plots ..." << std::endl;
  TFile fOut(outFile, "recreate");
  for (THnSparseD* const& h : hPreClusterInfoMULTI) {
    if (h)
      h->Write();
  }
  if (isTMC) {
    for (THnSparseD* const& h : hPreClusterInfoNoise) {
      if (h)
        h->Write();
    }
    for (THnSparseD* const& h : hPreClusterInfoTotal) {
      if (h)
        h->Write();
    }
  }
  if (isTMC) {
    for (int i = 0; i < 3; ++i) {
      if (h2ADCtrueVsADCfit[i])
        h2ADCtrueVsADCfit[i]->Write();
    }
    if (TCanvas* c = (TCanvas*)gROOT->FindObject("c_ADCtrueVsADCfit"))
      c->Write();
  }

  std::vector<std::string> canvasNames = {
    "c_chi2_ndf",
    "c_prob",
    "c_asymm"};
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
  cout << "discarded clusters (cut on ADCfit): " << discarded_cut_ADC_fit << " / " << nClusters << endl;
  cout << "TOTAL discarded clusters : " << (discarded_cut_ADC_fit + discarded_cut_ADC + discarded_cut_k3) << " / " << nClusters << endl;
}
