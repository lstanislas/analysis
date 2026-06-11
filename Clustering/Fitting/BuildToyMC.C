#include <array>
#include <cmath>
#include <chrono>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TSystem.h>
#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "Framework/Logger.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "PreClusterUtils.h"
#include "ToyMCUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;

//_________________________________________________________________________________________________
// run : run number
// inFile : root data file
// mode : "full" = use all clusters and do toyMC ; "cut" = use clusters that passed the selection and do toyMC
// fit : "none" = use cluster parameters from data ; "fit" = use clusters parameters from fit

// in XpX, p mean point (e.g. 2p34 == 2.34 , 0p4 == 0.4), useful for writting file name

// asymm : "none" = no asymmetry ; "copy" = copy the asymmetry from the data or from the fit; "gaus_XpX" = default asymm function in MC * XpX; "tripleGaus" = triple gaussians
// noise : "none" = no noise ; "MC_XpX" = gaussian noise with sigma = 0.5 * (sqrt(nSamples) + XpX) ; "MULT_XpX_XpX_XpX" = gaussian noise with sigma = XpX * sqrt(ADC) + XpX * ADC + XpX * sqrt(ADC) * ADC ; "predefined" = tuning per-station and cathode
// threshold : "none" = no threshold ; "gaus" = gaussian threshold ; "uniform" = static threshold
// k3x and k3y : change K3 values if positive (superseed other settings)
// try_tmc : redo ToyMC if the cluster isnt in the correct subspace (default = 50)
//_________________________________________________________________________________________________

/// store the new clusters together with the corresponding input data in outFile
/// require the MCH mapping to be loaded: gSystem->Load("libO2MCHGeometryTransformer"),  gSystem->Load("libO2MCHMappingImpl4"), gSystem->Load("libO2MCHTracking")

void BuildToyMC(int run, std::string inFile, std::string mode, std::string fit,
                std::string asymm, std::string noise, std::string threshold,
                double k3x = -1., double k3y = -1., int try_tmc = 50)
{

  if (mode != "full" && mode != "cut") {
    LOGP(error, "unknown simulation mode. Must be \"full\" or \"cut\"");
    exit(-1);
  }

  if (fit != "none" && fit != "fit") {
    LOGP(error, "unknown fit mode. Must be \"none\" or \"fit\"");
    exit(-1);
  }

  if (mode == "full" && fit == "fit") {
    LOGP(error, "fit mode incompatible with simulation mode \"full\"");
    exit(-1);
  }

  if (asymm != "none" && asymm != "copy" && !asymm.starts_with("gaus_") && asymm != "tripleGaus") {
    LOGP(error, "unknown asymmetry mode. Must be \"none\", \"copy\", \"gaus_XpX\" or \"tripleGaus\"");
    exit(-1);
  }

  if (noise != "none" && !noise.starts_with("MC_") && !noise.starts_with("MULT_") && noise != "predefined") {
    LOGP(error, "unknown noise mode. Must be \"none\", \"MC_XpX\", \"MULT_XpX_XpX_XpX\" or \"predefined\"");
    exit(-1);
  }

  if (threshold != "none" && threshold != "gaus" && threshold != "uniform") {
    LOGP(error, "unknown threshold mode. Must be \"none\", \"gaus\" or \"uniform\"");
    exit(-1);
  }

  // load CCDB objects
  InitFromCCDB(run, true, true, false);

  // load input data
  auto [dataFileIn, dataReader] = LoadData(inFile.c_str(), "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");
  std::unique_ptr<TTreeReaderArray<double>> fitParameters{};
  if (fit == "fit") {
    if (!dataReader->GetTree()->FindBranch("fitParameters")) {
      LOGP(error, "unable to load branch \"fitParameters\" from {}", inFile);
      exit(-1);
    }
    fitParameters = std::make_unique<TTreeReaderArray<double>>(*dataReader, "fitParameters");
  }

  // Create the output directory if it doesn't exist
  std::string outDir = "production/tmc";
  gSystem->mkdir(outDir.c_str(), true);

  // setup the output
  auto outFile = fmt::format("production/tmc/tmc_{}_{}_{}_{}_{}_{}_{}_{}.root", run, k3x, k3y, mode, fit, asymm, noise, threshold);
  TFile dataFileOut(outFile.c_str(), "recreate");
  TTree* dataTreeOut = new TTree("data", "tree tmc data");
  TrackParamStruct etrackParam;
  dataTreeOut->Branch("trackParameters", &etrackParam);
  int etrackTime;
  dataTreeOut->Branch("trackTime", &etrackTime);
  Cluster ecluster;
  dataTreeOut->Branch("clusters", &ecluster);
  std::vector<Digit> edigits;
  dataTreeOut->Branch("digits", &edigits);
  std::array<double, 6> parameters;
  dataTreeOut->Branch("parameters", &parameters);

  int nClusters = dataReader->GetEntries(false);
  int iCluster = 0;
  int selected = 0;
  int discarded = 0;
  int totalTries = 0;
  int retriedClusters = 0;
  auto tStart = std::chrono::high_resolution_clock::now();

  while (dataReader->Next()) {

    if (++iCluster % 10000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters << "..." << std::flush;
    }

    //___________________PRE-SELECTION__________________________
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

    // cut on precluster charge asymmetry
    auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
    if (std::abs(chargeAsymm) > 0.5) {
      continue;
    }

    // check if precluster pass the fit selection if needed
    if (mode == "cut" && !IsFittable(selectedDigits)) {
      continue;
    }

    ++selected;

    //___________________INIT PARAMETERS___________________________
    if (fit == "fit") {
      for (int i = 0; i < 6; ++i) {
        parameters[i] = (*fitParameters)[i];
      }
    } else {
      auto local = GlobalToLocal(cluster->getDEId(), cluster->x, cluster->y, cluster->z, run < 300000);
      // add charge fraction correction if no use of fit
      auto [chargeFracNB, chargeFracB] = GetChargeFraction(selectedDigits, local.x(), local.y());
      chargeNB /= chargeFracNB;
      chargeB /= chargeFracB;

      parameters[0] = local.x(); // X
      parameters[1] = local.y(); // Y
      parameters[2] = 0.3;       // K3X
      parameters[3] = 0.3;       // K3Y
      parameters[4] = chargeB;   // Qb_tot
      parameters[5] = chargeNB;  // Qnb_tot
    }

    // set K3X and K3Y for the predefined noise since the study was done with these values
    if (noise == "predefined") {
      int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
      if (iSt == 0) {
        parameters[2] = 0.3;
        parameters[3] = 0.32;
      } else if (iSt == 1) {
        parameters[2] = 0.46;
        parameters[3] = 0.43;
      } else if (iSt == 2) {
        parameters[2] = 0.35;
        parameters[3] = 0.33;
      }
    }

    // set K3X and K3Y to requested values (superseed other settings)
    if (k3x > 0.) {
      parameters[2] = k3x;
    }
    if (k3y > 0.) {
      parameters[3] = k3y;
    }

    //___________________RUN MC___________________________
    if (mode == "full") {

      edigits.clear();
      TMC(edigits, *trackTime, cluster->getDEId(), parameters, asymm, noise, threshold);
      if (edigits.empty() || IsMonoCathode(edigits)) {
        ++discarded;
        continue;
      }
    } else {

      int tries = 0;
      do {
        edigits.clear();
        TMC(edigits, *trackTime, cluster->getDEId(), parameters, asymm, noise, threshold);
        ++tries;
      } while (!IsFittable(edigits) && tries < try_tmc);
      if (!IsFittable(edigits)) {
        ++discarded;
        continue;
      }
      totalTries += tries;
      if (tries > 1) {
        ++retriedClusters;
      }
    }

    //___________________SAVE OUTPUT___________________________
    // TMC process will associate the same tracks and time as data on the corresponding precluster
    etrackParam = *trackParam;
    etrackTime = *trackTime;
    ecluster = MakeCluster(cluster->uid, parameters[0], parameters[1]);
    dataTreeOut->Fill();
  }

  auto tEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> timer = tEnd - tStart;
  cout << "\r\033[Kprocessing completed. Duration = " << timer.count() << " s" << endl;
  cout << "selected clusters = " << selected << " / " << nClusters << endl;
  cout << "discarded clusters : " << discarded << " / " << selected << endl;
  if (mode == "cut") {
    int accepted = selected - discarded;
    cout << "clusters needing at least 1 retry : " << retriedClusters << " / " << accepted
         << fmt::format(" ({:.1f}%)", 100. * retriedClusters / accepted) << endl;
    cout << fmt::format("average TMC attempts per accepted cluster : {:.2f}", static_cast<double>(totalTries) / accepted) << endl;
  }
  dataFileOut.Write("", TObject::kOverwrite);
  dataFileOut.Close();
  dataFileIn->Close();
  cout << "input file : " << inFile << endl;
  cout << "output file : " << outFile << endl;
}
