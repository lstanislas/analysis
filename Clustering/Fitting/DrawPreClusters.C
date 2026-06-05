#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "CommonUtils/ConfigurableParam.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/Logger.h"
#include "MCHBase/TrackBlock.h"

#include "CCDBUtils.h"
#include "ClusterUtils.h"
#include "DataUtils.h"
#include "DigitUtils.h"
#include "PreClusterUtils.h"

using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::TrackParamStruct;

static constexpr double pi = 3.14159265358979323846;

void SetupRun2Mathieson();

//_________________________________________________________________________________________________
void DrawPreClusters(int run, bool applyTrackSelection = false, bool applyClusterSelection = false,
                     bool applyTimeSelection = false, bool correctCharge = false,
                     bool useFitPos = false, bool useFitCharge = false,
                     const char* inFile = "clusters.root", const char* outFile = "displays.root")
{
  /// draw precluster and associated digits informations
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  if (run < 300000) {
    o2::base::GeometryManager::loadGeometry("O2geometry.root");
    SetupRun2Mathieson();
  } else {
    InitFromCCDB(run, true, true, false);
  }

  auto [dataFileIn, dataReader] = LoadData(inFile, "data");
  TTreeReaderValue<TrackParamStruct> trackParam(*dataReader, "trackParameters");
  TTreeReaderValue<int> trackTime(*dataReader, "trackTime");
  TTreeReaderValue<Cluster> cluster(*dataReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> digits(*dataReader, "digits");
  std::unique_ptr<TTreeReaderArray<double>> fitParameters{};
  if (useFitPos || useFitCharge) {
    if (!dataReader->GetTree()->FindBranch("fitParameters")) {
      LOGP(error, "unable to load branch \"fitParameters\" from {}", inFile);
      exit(-1);
    }
    fitParameters = std::make_unique<TTreeReaderArray<double>>(*dataReader, "fitParameters");
  }

  std::vector<TH1*> preClusterInfo{};
  CreatePreClusterInfo(preClusterInfo);
  std::vector<TH1*> preClusterInfoSt[3] = {{}, {}, {}};
  CreatePreClusterInfo(preClusterInfoSt[0], "St1");
  CreatePreClusterInfo(preClusterInfoSt[1], "St2");
  CreatePreClusterInfo(preClusterInfoSt[2], "St345");

  std::vector<TH1*> preClusterInfoVsWireSt[3] = {{}, {}, {}};
  CreatePreClusterInfoVsWire(preClusterInfoVsWireSt[0], "St1");
  CreatePreClusterInfoVsWire(preClusterInfoVsWireSt[1], "St2");
  CreatePreClusterInfoVsWire(preClusterInfoVsWireSt[2], "St345");

  TH3* hPreClusterInfo3D[4] = {nullptr, nullptr, nullptr, nullptr};
  hPreClusterInfo3D[0] = CreatePreClusterInfo3D("St1");
  hPreClusterInfo3D[1] = CreatePreClusterInfo3D("St2");
  hPreClusterInfo3D[2] = CreatePreClusterInfo3D("St345");
  hPreClusterInfo3D[3] = CreatePreClusterInfo3D();

  std::vector<TH1*> digitTimeInfo{};
  CreateDigitTimeInfo(digitTimeInfo);
  std::vector<TH1*> digitChargeInfo{};
  CreateDigitChargeInfo(digitChargeInfo);
  std::vector<TH1*> digitChargeInfoSt[3] = {{}, {}, {}};
  CreateDigitChargeInfo(digitChargeInfoSt[0], "St1");
  CreateDigitChargeInfo(digitChargeInfoSt[1], "St2");
  CreateDigitChargeInfo(digitChargeInfoSt[2], "St345");

  int nClusters = dataReader->GetEntries(false);
  int iCluster = 0;
  while (dataReader->Next()) {
    if (++iCluster % 100000 == 0) {
      std::cout << "\rprocessing cluster " << iCluster << " / " << nClusters << "..." << std::flush;
    }

    // those 2 DE have lower HV for the run 529691
    if (run == 529691 && (cluster->getDEId() == 202 || cluster->getDEId() == 300)) {
      continue;
    }

    // cut on track angle at chamber
    if (applyTrackSelection && std::abs(std::atan2(trackParam->py, -trackParam->pz)) / pi * 180. > 10.) {
      continue;
    }

    // cut on track sign
    // if (applyTrackSelection && trackParam->sign > 0) {
    //   continue;
    // }

    // cut on digit time
    std::vector<Digit> selectedDigits(*digits);
    if (applyTimeSelection) {
      selectedDigits.erase(
        std::remove_if(selectedDigits.begin(), selectedDigits.end(), [&trackTime](const auto& digit) {
          return std::abs(digit.getTime() + 1.5 - *trackTime) > 10.;
        }),
        selectedDigits.end());
      if (selectedDigits.empty()) {
        continue;
      }
    }

    // reject mono-cathode clusters after digit selection
    const auto [sizeX, sizeY] = GetSize(selectedDigits);
    if (sizeX == 0 || sizeY == 0) {
      continue;
    }

    // reject composite preclusters
    if (IsComposite(selectedDigits, true)) {
      continue;
    }

    // cut on distance to closest wire
    double localX, localY;
    if (useFitPos) {
      localX = (*fitParameters)[0];
      localY = (*fitParameters)[1];
    } else {
      auto local = GlobalToLocal(cluster->getDEId(), cluster->x, cluster->y, cluster->z, run < 300000);
      localX = local.x();
      localY = local.y();
    }
    float dx = DistanceToClosestWire(cluster->getDEId(), localX);
    // float dx = DistanceToClosestWire(cluster->getDEId(), trackParam->x, trackParam->y, trackParam->z, run < 300000);
    // if (applyClusterSelection && std::abs(dx) > 0.015) {
    // if (applyClusterSelection && std::abs(dx) < 0.085) {
    //   continue;
    // }

    // cut on cluster charge
    auto [chargeNB, chargeB] = GetCharge(selectedDigits, run < 300000);
    // double charge = 0.5 * (chargeNB + chargeB);
    // if (applyClusterSelection && charge < 4000.) {
    //   continue;
    // }

    // cut on cluster charge asymmetry
    double chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
    if (applyClusterSelection && std::abs(chargeAsymm) > 0.5) {
      continue;
    }
    // if (applyClusterSelection && (chargeAsymm >= -0.2 || chargeAsymm < -0.3)) {
    //   continue;
    // }

    // use fitted charge or correct pad charge and re-cut on cluster charge asymmetry
    if (correctCharge || useFitCharge) {
      if (useFitCharge) {
        chargeB = (*fitParameters)[4];
        chargeNB = (*fitParameters)[5];
      } else {
        auto [chargeFracNB, chargeFracB] = GetChargeFraction(selectedDigits, localX, localY);
        chargeNB /= chargeFracNB;
        chargeB /= chargeFracB;
      }
      chargeAsymm = (chargeNB - chargeB) / (chargeNB + chargeB);
      if (applyClusterSelection && std::abs(chargeAsymm) > 0.5) {
        continue;
      }
    }

    // cut on cluster size
    // if (applyClusterSelection && sizeX < 2) {
    // if (applyClusterSelection && (sizeX < 2 || sizeY < 2)) {
    //   continue;
    // }

    // cut on cluster size asymmetry
    // if (applyClusterSelection && (sizeY > sizeX + 3 || sizeX > sizeY + 2)) {
    //   continue;
    // }

    // check if precluster pass the fit selection
    // if (applyClusterSelection && !IsFittable(selectedDigits)) {
    //   continue;
    // }

    // double maxChargeB = -1;
    // Digit* digitMaxB = nullptr;
    // double maxChargeNB = -1;
    // Digit* digitMaxNB = nullptr;
    // for (auto& digit : selectedDigits) {
    //   const auto& segmentation = o2::mch::mapping::segmentation(digit.getDetID());
    //   if (segmentation.isBendingPad(digit.getPadID())) {
    //     if (digit.getADC() > maxChargeB) {
    //       digitMaxB = &digit;
    //       maxChargeB = digit.getADC();
    //     }
    //   } else {
    //     if (digit.getADC() > maxChargeNB) {
    //       digitMaxNB = &digit;
    //       maxChargeNB = digit.getADC();
    //     }
    //   }
    // }

    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfo);
    FillPreClusterInfo3D(chargeNB, chargeB, dx, hPreClusterInfo3D[3]);
    int iSt = (cluster->getChamberId() < 4) ? cluster->getChamberId() / 2 : 2;
    FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfoSt[iSt]);
    FillPreClusterInfo3D(chargeNB, chargeB, dx, hPreClusterInfo3D[iSt]);
    FillPreClusterInfoVsWire(chargeNB, chargeB, dx, preClusterInfoVsWireSt[iSt]);

    for (const auto& digit : selectedDigits) {
      // if (!(&digit == digitMaxB || &digit == digitMaxNB)) {
      //   continue;
      // }
      FillDigitTimeInfo(digit, *trackTime, digitTimeInfo);
      FillDigitChargeInfo(digit, digitChargeInfo, chargeAsymm, run < 300000);
      FillDigitChargeInfo(digit, digitChargeInfoSt[iSt], chargeAsymm, run < 300000);
    }
  }
  cout << "\r\033[Kprocessing completed" << endl;

  dataFileIn->Close();

  // display
  gStyle->SetOptStat(1);

  auto c = DrawPreClusterInfo(preClusterInfo);
  auto cSt1 = DrawPreClusterInfo(preClusterInfoSt[0], "St1");
  auto cSt2 = DrawPreClusterInfo(preClusterInfoSt[1], "St2");
  auto cSt345 = DrawPreClusterInfo(preClusterInfoSt[2], "St345");

  auto cwSt1 = DrawPreClusterInfoVsWire(preClusterInfoVsWireSt[0], "St1");
  auto cwSt2 = DrawPreClusterInfoVsWire(preClusterInfoVsWireSt[1], "St2");
  auto cwSt345 = DrawPreClusterInfoVsWire(preClusterInfoVsWireSt[2], "St345");

  auto ct = DrawDigitTimeInfo(digitTimeInfo);
  auto cc = DrawDigitChargeInfo(digitChargeInfo);
  auto ccSt1 = DrawDigitChargeInfo(digitChargeInfoSt[0], "St1");
  auto ccSt2 = DrawDigitChargeInfo(digitChargeInfoSt[1], "St2");
  auto ccSt345 = DrawDigitChargeInfo(digitChargeInfoSt[2], "St345");

  // output
  TFile fOut(outFile, "recreate");
  for (auto h : hPreClusterInfo3D) {
    h->Write();
  }
  c->Write();
  cSt1->Write();
  cSt2->Write();
  cSt345->Write();
  cwSt1->Write();
  cwSt2->Write();
  cwSt345->Write();
  ct->Write();
  cc->Write();
  ccSt1->Write();
  ccSt2->Write();
  ccSt345->Write();
  fOut.Close();
}

//_________________________________________________________________________________________________
void SetupRun2Mathieson()
{
  /// use run2 Mathieson parameterizations

  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKx3St1", "0.7000");
  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKy3St1", "0.7550");

  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKx3St2345", "0.7131");
  o2::conf::ConfigurableParam::setValue("MCHResponse.mathiesonSqrtKy3St2345", "0.7642");
}
