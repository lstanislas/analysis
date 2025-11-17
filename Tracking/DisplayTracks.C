#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include <unordered_map>

#include <fmt/format.h>

#include <gsl/span>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TDatabasePDG.h>
#include <TGeoManager.h>
#include <Math/Vector4D.h>

#include "CCDB/BasicCCDBManager.h"
#include "Framework/Logger.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CommonUtils/NameConf.h"
#include "CommonConstants/LHCConstants.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"

using namespace o2;

uint32_t nOrbitsPerTF = 128;

// first orbit of the first TF of each run
const std::unordered_map<uint32_t, uint32_t> firstTForbit0perRun{
  {505207, 133875},
  {505217, 14225007},
  {505278, 1349340},
  {505285, 1488862},
  {505303, 2615411},
  {505397, 5093945},
  {505404, 19196217},
  {505405, 28537913},
  {505406, 41107641},
  {505413, 452530},
  {505440, 13320708},
  {505443, 26546564},
  {505446, 177711},
  {505548, 88037114},
  {505582, 295044346},
  {505600, 417241082},
  {505623, 10445984},
  {505629, 126979},
  {505637, 338969},
  {505645, 188222},
  {505658, 81044},
  {505669, 328291},
  {505673, 30988},
  {505713, 620506},
  {505720, 5359903},
  {514053, 10919},
  {521520, 84669824},
  {521684, 49607462},
  {529270, 68051456},
  {529414, 860032},
  {529691, 74406400},
  {544490, 26693600},
  {549586, 67400512},
  {559410, 527962048},
  {562849, 179331744},
  {562967, 47243328},
  {562968, 62076224},
  {562970, 86353952},
  {562971, 97188864},
  {562975, 146579680},
  {562976, 160383072}};

struct TrackInfo {
  TrackInfo(const mch::TrackMCH& mch) : mchTrack(mch) {}

  const mch::TrackMCH& mchTrack;
  mch::TrackParam paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  std::vector<const mch::Digit*> digits{};
  double mchTime = -1.;
  double mchTimeRMS = 0.;
  double mchTimeSt12 = -1.;
  double mchTimeRMSSt12 = 0.;
  double mchTimeSt345 = -1.;
  double mchTimeRMSSt345 = 0.;
  std::vector<const mch::Digit*> digitsAtClusterPos{};
  double mchTimeAtClusterPos = -1.;
  double mchTimeRMSAtClusterPos = 0.;
  double mchTimeAtClusterPosSt12 = -1.;
  double mchTimeRMSAtClusterPosSt12 = 0.;
  double mchTimeAtClusterPosSt345 = -1.;
  double mchTimeRMSAtClusterPosSt345 = 0.;
  int midTime = -1;
};

o2::mch::geo::TransformationCreator transformation;
static const double muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
uint16_t minNSamplesSignal = 17;
double signalParam[4] = {80., 16., 12., 1.2};
uint16_t minNSamplesBackground = 14;
double backgroundParam[4] = {18., 24., -20., 7.};
int bcIntegrationRange = 10; // time window ([-range, range]) to integrate digits
int minNDigitsSignal = 10; // minimum number of digits passing the signal cuts to select signal events

constexpr double pi() { return 3.14159265358979323846; }
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
void LoadDigits(TrackInfo& trackInfo, const std::vector<mch::Cluster>& clusters, const std::vector<mch::Digit>& digits,
                bool selectSignal, bool rejectBackground);
uint32_t IRtoBCinTF(InteractionRecord ir, uint32_t firstTForbit0, bool isMC);
void computeMCHTime(const std::vector<const mch::Digit*>& digits, double& mean, double& rms, int deMin = 100, int deMax = 1025);
const dataformats::TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<dataformats::TrackMCHMID>& muonTracks);
bool ExtrapToVertex(TrackInfo& trackInfo);
bool IsSelected(TrackInfo& trackInfo);
bool IsSignal(TrackInfo& trackInfo);
bool IsReconstructible(TrackInfo& trackInfo);
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension);
void FillHistosAtVertex(const TrackInfo& trackInfo, std::vector<TH1*>& histos);
void DrawHistosAtVertex(std::vector<TH1*> histos[2]);
void CreateTimeHistos(std::vector<TH1*>& histos, const char* extension);
void FillTimeHistos(const std::vector<const mch::Digit*>& digits, double mchTime, double mchTimeRMS, int midTime,
                    gsl::span<TH1*> histos, int deMin = 100, int deMax = 1025);
void DrawTimeHistos(gsl::span<TH1*> histos, const char* extension);
void CreateROFTimeHistos(std::vector<TH1*>& histos);
void FillROFTimeHistos(int rofTime[2], double mchTime, int midTime, gsl::span<TH1*> histos);
void DrawROFTimeHistos(gsl::span<TH1*> histos);
void CreateChargeHistos(std::vector<TH1*>& histos, const char* extension);
void FillChargeHistos(const std::vector<const mch::Digit*>& digits, gsl::span<TH1*> histos, int deMin = 100, int deMax = 1025);
void DrawChargeHistos(gsl::span<TH1*> histos, const char* extension);
void CreateCorrelationHistos(std::vector<TH1*>& histos);
void FillCorrelationHistos(const std::vector<const mch::Digit*>& digits, TH1* hist, double timeRef, int deMin = 100, int deMax = 1025);
void DrawCorrelationHistos(std::vector<TH1*>& histos);
void CreateMultHistos(std::vector<TH1*>& histos);
void FillMultHistos(const mch::ROFRecord& rof, const std::vector<mch::Cluster>& clusters, std::vector<TH1*>& histos);
void DrawMultHistos(std::vector<TH1*>& histos);
double signalCut(double* x, double* p);
double backgroundCut(double* x, double* p);
void WriteHistos(TFile* f, const char* dirName, const std::vector<TH1*>& histos);

//_________________________________________________________________________________________________
void DisplayTracks(int runNumber, std::string clusterFileName = "",
                   std::string mchFileName = "mchtracks.root", std::string muonFileName = "muontracks.root",
                   bool applyTrackSelection = false, bool selectSignal = false, bool rejectBackground = false,
                   std::string outFileName = "", bool isMC = false)
{
  /// show the characteristics of the reconstructed tracks
  /// store the ouput histograms in outFileName if any

  // make sure the correct mapping is loaded
  auto& segmentation = mch::mapping::segmentation(300);
  if (segmentation.nofPads() != 27873) {
    LOG(error) << "wrong mapping implementation";
    LOG(error) << "do gSystem->Load(\"libO2MCHMappingImpl4\") before compiling this macro";
    exit(-1);
  }

  // load magnetic field (from o2sim_grp.root) and geometry (from o2sim_geometry.root)
  // const auto grp = parameters::GRPObject::loadFrom(base::NameConf::getGRPFileName());
  // base::GeometryManager::loadGeometry();

  // load magnetic field and geometry from CCDB
  auto& ccdb = ccdb::BasicCCDBManager::instance();
  auto [tStart, tEnd] = ccdb.getRunDuration(runNumber);
  ccdb.setTimestamp(tEnd);
  // ccdb.setURL("file://ccdb");
  auto grp = ccdb.get<parameters::GRPMagField>("GLO/Config/GRPMagField");
  // ccdb.setURL("http://localhost:6060");
  auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");

  // and prepare track extrapolation to vertex (0,0,0)
  base::Propagator::initFieldFromGRP(grp);
  mch::TrackExtrap::setField();
  mch::TrackExtrap::useExtrapV2();
  transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

  // find the first orbit of the first TF of this run
  uint32_t firstTForbit0(0);
  auto itOrbit0 = firstTForbit0perRun.find(runNumber);
  if (itOrbit0 != firstTForbit0perRun.end()) {
    firstTForbit0 = itOrbit0->second;
  } else {
    LOG(warning) << "first orbit not found for this run";
  }

  // find the number of orbits per TF for this run
  if (runNumber >= 544490) {
    nOrbitsPerTF = 32;
  }

  // load tracks
  auto [fMCH, mchReader] = LoadData(mchFileName.c_str(), "o2sim");
  TTreeReaderValue<std::vector<mch::ROFRecord>> mchROFs = {*mchReader, "trackrofs"};
  TTreeReaderValue<std::vector<mch::TrackMCH>> mchTracks = {*mchReader, "tracks"};
  TTreeReaderValue<std::vector<mch::Cluster>> mchTrackClusters = {*mchReader, "trackclusters"};
  TFile* fClusters = nullptr;
  TTreeReader* clusterReader = nullptr;
  std::unique_ptr<TTreeReaderValue<std::vector<mch::Digit>>> mchDigits{};
  std::unique_ptr<TTreeReaderValue<std::vector<mch::Cluster>>> mchClusters{};
  std::unique_ptr<TTreeReaderValue<std::vector<mch::ROFRecord>>> mchClusterROFs{};
  if (clusterFileName.empty()) {
    //read digits from the mchtracks.root
    if (mchReader->GetTree()->FindBranch("trackdigits")) {
      mchDigits = std::make_unique<TTreeReaderValue<std::vector<mch::Digit>>>(*mchReader, "trackdigits");
    } else {
      LOG(error) << "track clusters point to digits in the cluster file, which is not provided";
    }
  } else {
    // read digits from the mchclusters.root
    if (mchReader->GetTree()->FindBranch("trackdigits")) {
      LOG(error) << "track clusters point to digits in the track file";
      exit(-1);
    }
    std::tie(fClusters, clusterReader) = LoadData(clusterFileName.c_str(), "o2sim");
    mchDigits = std::make_unique<TTreeReaderValue<std::vector<mch::Digit>>>(*clusterReader, "clusterdigits");
    mchClusters = std::make_unique<TTreeReaderValue<std::vector<mch::Cluster>>>(*clusterReader, "clusters");
    mchClusterROFs = std::make_unique<TTreeReaderValue<std::vector<mch::ROFRecord>>>(*clusterReader, "clusterrofs");
  }
  auto [fMUON, muonReader] = LoadData(muonFileName.c_str(), "o2sim");
  TTreeReaderValue<std::vector<dataformats::TrackMCHMID>> muonTracks = {*muonReader, "tracks"};
  int nTF = muonReader->GetEntries(false);
  if (mchReader->GetEntries(false) != nTF) {
    LOG(error) << mchFileName << " and " << muonFileName << " do not contain the same number of TF";
    exit(-1);
  }

  // open output file if any
  TFile* fOut = nullptr;
  if (!outFileName.empty()) {
    fOut = TFile::Open(outFileName.c_str(), "RECREATE");
  }

  // create histograms
  std::vector<TH1*> histosAtVertex[2] = {{}, {}};
  CreateHistosAtVertex(histosAtVertex[0], "mch");
  CreateHistosAtVertex(histosAtVertex[1], "muon");
  TH1* hmatchChi2 = new TH1F("matchChi2", "normalized matching #chi^{2};#chi^{2} / ndf", 500, 0., 50.);
  hmatchChi2->SetDirectory(0);
  TH1* hNClustersPerCh = new TH1F("nClustersPerCh", "number of clusters per muon per chamber;chamber", 10, 0.5, 10.5);
  hNClustersPerCh->SetDirectory(0);
  std::vector<TH1*> timeHistos{};
  CreateTimeHistos(timeHistos, "AllDig");
  CreateTimeHistos(timeHistos, "");
  std::vector<TH1*> timeHistosSt12{};
  CreateTimeHistos(timeHistosSt12, "St12AllDig");
  CreateTimeHistos(timeHistosSt12, "St12");
  std::vector<TH1*> timeHistosSt345{};
  CreateTimeHistos(timeHistosSt345, "St345AllDig");
  CreateTimeHistos(timeHistosSt345, "St345");
  std::vector<TH1*> rofTimeHistos{};
  CreateROFTimeHistos(rofTimeHistos);
  std::vector<TH1*> chargeHistos{};
  CreateChargeHistos(chargeHistos, "AllDig");
  CreateChargeHistos(chargeHistos, "");
  std::vector<TH1*> chargeHistosSt1{};
  CreateChargeHistos(chargeHistosSt1, "St1AllDig");
  CreateChargeHistos(chargeHistosSt1, "St1");
  std::vector<TH1*> chargeHistosSt2{};
  CreateChargeHistos(chargeHistosSt2, "St2AllDig");
  CreateChargeHistos(chargeHistosSt2, "St2");
  std::vector<TH1*> chargeHistosSt345{};
  CreateChargeHistos(chargeHistosSt345, "St345AllDig");
  CreateChargeHistos(chargeHistosSt345, "St345");
  std::vector<TH1*> corrHistos{};
  CreateCorrelationHistos(corrHistos);
  std::vector<TH1*> multHistos{};
  CreateMultHistos(multHistos);
  TH1F* hMass = new TH1F("mass", "#mu^{+}#mu^{-} invariant mass;mass (GeV/c^{2})", 1600, 0., 20.);
  hMass->SetDirectory(0);

  int iTF = 0;
  int nROFs[4] = {0, 0, 0, 0};
  while (mchReader->Next() && muonReader->Next() && (clusterReader == nullptr || clusterReader->Next())) {
    if (++iTF % 1000 == 0) {
      std::cout << "\rprocessing TF " << iTF << " / " << nTF << "..." << std::flush;
    }

    if (mchClusterROFs && (*mchClusterROFs)->size() != mchROFs->size()) {
      LOGP(error, "{} and {} do not contain the same number of ROFs in TF {}", mchFileName, clusterFileName, iTF);
      exit(-1);
    }

    // count total number of ROFs
    nROFs[0] += mchROFs->size();

    for (size_t iROF = 0; iROF < mchROFs->size(); ++iROF) {
      const auto& mchROF = (*mchROFs)[iROF];

      // count number of ROFs with MCH track(s)
      if (mchROF.getNEntries() > 0) {
        nROFs[1]++;
      }

      int rofTime[2] = {0, 0};
      rofTime[0] = IRtoBCinTF(mchROF.getBCData(), firstTForbit0, isMC);
      rofTime[1] = rofTime[0] + mchROF.getBCWidth() - 1;

      // if (mchROF.getBCWidth() != 1) {
      //   continue;
      // }

      // if (mchROF.getBCData() != InteractionRecord(3293, 51384998)) {
      //   continue;
      // }

      std::vector<ROOT::Math::PxPyPzMVector> muVector{};
      std::vector<int> muSign{};
      bool containsMuon = false;
      bool containsSelectedMuon = false;
      for (int iMCHTrack = mchROF.getFirstIdx(); iMCHTrack <= mchROF.getLastIdx(); ++iMCHTrack) {

        // find the corresponding MCH-MID matched track
        auto muon = FindMuon(iMCHTrack, *muonTracks);
        if (muon) {
          containsMuon = true;
        }

        TrackInfo trackInfo((*mchTracks)[iMCHTrack]);
        // if (trackInfo.mchTrack.getNClusters() != 10) {
        //   continue;
        // }

        // compute the track parameters at vertex
        if (!ExtrapToVertex(trackInfo)) {
          LOG(error) << "track extrapolation to vertex failed";
          continue;
        }

        // apply track selection if requested
        if (applyTrackSelection && !IsSelected(trackInfo)) {
          continue;
        }

        // fill MID timing info
        if (muon) {
          containsSelectedMuon = true;
          trackInfo.midTime = IRtoBCinTF(muon->getIR(), firstTForbit0, isMC);
          // check that MID time is within MCH ROF boundaries
          if (trackInfo.midTime < rofTime[0] || trackInfo.midTime > rofTime[1]) {
            LOGP(error, "MID time {} outside ROF boundaries [{},{}]", trackInfo.midTime, rofTime[0], rofTime[1]);
          }
        }

        // fill digit info
        if (mchDigits) {
          LoadDigits(trackInfo, *mchTrackClusters, **mchDigits, selectSignal, rejectBackground);
          computeMCHTime(trackInfo.digits, trackInfo.mchTime, trackInfo.mchTimeRMS);
          // float dt = (*mchTracks)[iMCHTrack].getTimeMUS().getTimeStamp() - static_cast<float>(o2::constants::lhc::LHCBunchSpacingMUS * trackInfo.mchTime);
          // if (dt != 0.f) {
          //   LOGP(info, "dt = {}", dt);
          // }
          computeMCHTime(trackInfo.digitsAtClusterPos, trackInfo.mchTimeAtClusterPos, trackInfo.mchTimeRMSAtClusterPos);
          computeMCHTime(trackInfo.digits, trackInfo.mchTimeSt12, trackInfo.mchTimeRMSSt12, 100, 403);
          computeMCHTime(trackInfo.digitsAtClusterPos, trackInfo.mchTimeAtClusterPosSt12, trackInfo.mchTimeRMSAtClusterPosSt12, 100, 403);
          computeMCHTime(trackInfo.digits, trackInfo.mchTimeSt345, trackInfo.mchTimeRMSSt345, 500, 1025);
          computeMCHTime(trackInfo.digitsAtClusterPos, trackInfo.mchTimeAtClusterPosSt345, trackInfo.mchTimeRMSAtClusterPosSt345, 500, 1025);
        }

        // apply digit selection if requested
        if (selectSignal && !IsSignal(trackInfo)) {
          continue;
        }
        if (rejectBackground && !IsReconstructible(trackInfo)) {
          continue;
        }

        // check that digits' time is within ROF boundaries (not valid in trigger mode)
        for (const auto digit : trackInfo.digits) {
          if (digit->getTime() < rofTime[0] || digit->getTime() > rofTime[1]) {
            LOGP(error, "digit time {} outside ROF boundaries [{},{}]", digit->getTime(), rofTime[0], rofTime[1]);
          }
        }

        // fill histograms
        FillHistosAtVertex(trackInfo, histosAtVertex[0]);
        FillChargeHistos(trackInfo.digits, {&chargeHistos[0], 3});
        FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistos[6], 3});
        FillChargeHistos(trackInfo.digits, {&chargeHistosSt1[0], 3}, 100, 203);
        FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistosSt1[6], 3}, 100, 203);
        FillChargeHistos(trackInfo.digits, {&chargeHistosSt2[0], 3}, 300, 403);
        FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistosSt2[6], 3}, 300, 403);
        FillChargeHistos(trackInfo.digits, {&chargeHistosSt345[0], 3}, 500, 1025);
        FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistosSt345[6], 3}, 500, 1025);
        if (muon) {
          FillHistosAtVertex(trackInfo, histosAtVertex[1]);
          hmatchChi2->Fill(muon->getMatchChi2OverNDF());
          for (int iCl = trackInfo.mchTrack.getFirstClusterIdx(); iCl <= trackInfo.mchTrack.getLastClusterIdx(); ++iCl) {
            hNClustersPerCh->Fill((*mchTrackClusters)[iCl].getChamberId() + 1);
          }
          FillChargeHistos(trackInfo.digits, {&chargeHistos[3], 3});
          FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistos[9], 3});
          FillChargeHistos(trackInfo.digits, {&chargeHistosSt1[3], 3}, 100, 203);
          FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistosSt1[9], 3}, 100, 203);
          FillChargeHistos(trackInfo.digits, {&chargeHistosSt2[3], 3}, 300, 403);
          FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistosSt2[9], 3}, 300, 403);
          FillChargeHistos(trackInfo.digits, {&chargeHistosSt345[3], 3}, 500, 1025);
          FillChargeHistos(trackInfo.digitsAtClusterPos, {&chargeHistosSt345[9], 3}, 500, 1025);
          FillCorrelationHistos(trackInfo.digits, corrHistos[0], trackInfo.mchTime);
          FillCorrelationHistos(trackInfo.digits, corrHistos[1], trackInfo.mchTime, 100, 203);
          FillCorrelationHistos(trackInfo.digits, corrHistos[2], trackInfo.mchTime, 300, 403);
          FillCorrelationHistos(trackInfo.digits, corrHistos[3], trackInfo.mchTime, 500, 1025);
          muVector.emplace_back(trackInfo.paramAtVertex.px(), trackInfo.paramAtVertex.py(), trackInfo.paramAtVertex.pz(), muMass);
          muSign.emplace_back((trackInfo.paramAtVertex.getCharge() > 0) ? 1 : -1);
        }
        FillTimeHistos(trackInfo.digits, trackInfo.mchTime, trackInfo.mchTimeRMS,
                       trackInfo.midTime, {&timeHistos[0], 6});
        FillTimeHistos(trackInfo.digitsAtClusterPos, trackInfo.mchTimeAtClusterPos, trackInfo.mchTimeRMSAtClusterPos,
                       trackInfo.midTime, {&timeHistos[6], 6});
        FillTimeHistos(trackInfo.digits, trackInfo.mchTimeSt12, trackInfo.mchTimeRMSSt12,
                       trackInfo.midTime, {&timeHistosSt12[0], 6}, 100, 403);
        FillTimeHistos(trackInfo.digitsAtClusterPos, trackInfo.mchTimeAtClusterPosSt12, trackInfo.mchTimeRMSAtClusterPosSt12,
                       trackInfo.midTime, {&timeHistosSt12[6], 6}, 100, 403);
        FillTimeHistos(trackInfo.digits, trackInfo.mchTimeSt345, trackInfo.mchTimeRMSSt345,
                       trackInfo.midTime, {&timeHistosSt345[0], 6}, 500, 1025);
        FillTimeHistos(trackInfo.digitsAtClusterPos, trackInfo.mchTimeAtClusterPosSt345, trackInfo.mchTimeRMSAtClusterPosSt345,
                       trackInfo.midTime, {&timeHistosSt345[6], 6}, 500, 1025);
        FillROFTimeHistos(rofTime, trackInfo.mchTime, trackInfo.midTime, rofTimeHistos);
      }

      // count number of ROFs with muon track(s)
      if (containsMuon) {
        nROFs[2]++;
      }

      // count number of ROFs with selected muon track(s) and fill cluster multiplicity histograms for them
      if (containsSelectedMuon) {
        nROFs[3]++;
        if (mchClusterROFs) {
          FillMultHistos((**mchClusterROFs)[iROF], **mchClusters, multHistos);
        }
      }

      if (muVector.size() > 1) {
        for (size_t i = 0; i < muVector.size(); ++i) {
          for (size_t j = i + 1; j < muVector.size(); ++j) {
            if (muSign[i] * muSign[j] < 0) {
              ROOT::Math::PxPyPzMVector dimu = muVector[i] + muVector[j];
              hMass->Fill(dimu.M());
            }
          }
        }
      }
    }
  }
  std::cout << "\r\033[Kprocessing " << nTF << " TF completed" << std::endl;

  std::cout << "total number of ROFs = " << nROFs[0] << std::endl;
  std::cout << "with tracks          = " << nROFs[1] << " (" << 100. * nROFs[1] / nROFs[0] << " %)" << std::endl;
  std::cout << "with muons           = " << nROFs[2] << " (" << 100. * nROFs[2] / nROFs[0] << " %)" << std::endl;
  std::cout << "with selected muons  = " << nROFs[3] << " (" << 100. * nROFs[3] / nROFs[0] << " %)" << std::endl;

  fMCH->Close();
  if (fClusters != nullptr) {
    fClusters->Close();
  }
  fMUON->Close();

  hNClustersPerCh->Scale(1. / histosAtVertex[1][0]->GetEntries());

  // store histograms if requested
  if (fOut) {
    fOut->cd();
    TParameter<int> p1("nTF", nTF);
    p1.Write();
    TParameter<int> p2("nROFSelected", nROFs[3]);
    p2.Write();
    WriteHistos(fOut, "general", histosAtVertex[0]);
    WriteHistos(fOut, "general", histosAtVertex[1]);
    hmatchChi2->Write();
    hNClustersPerCh->Write();
    WriteHistos(fOut, "time", timeHistos);
    WriteHistos(fOut, "time", timeHistosSt12);
    WriteHistos(fOut, "time", timeHistosSt345);
    WriteHistos(fOut, "time", rofTimeHistos);
    WriteHistos(fOut, "charge", chargeHistos);
    WriteHistos(fOut, "charge", chargeHistosSt1);
    WriteHistos(fOut, "charge", chargeHistosSt2);
    WriteHistos(fOut, "charge", chargeHistosSt345);
    WriteHistos(fOut, "correlations", corrHistos);
    WriteHistos(fOut, "multiplicity", multHistos);
    fOut->Close();
  }

  // display histograms
  DrawHistosAtVertex(histosAtVertex);
  TCanvas* cMatch = new TCanvas();
  gPad->SetLogy();
  hmatchChi2->Draw();
  TCanvas* cNClustersPerCh = new TCanvas();
  hNClustersPerCh->Draw();
  DrawTimeHistos({&timeHistos[0], 6}, "AllDigits");
  DrawTimeHistos({&timeHistos[6], 6}, "DigitsAtClusterPos");
  DrawTimeHistos({&timeHistosSt12[0], 6}, "St12_AllDigits");
  DrawTimeHistos({&timeHistosSt12[6], 6}, "St12_DigitsAtClusterPos");
  DrawTimeHistos({&timeHistosSt345[0], 6}, "St345_AllDigits");
  DrawTimeHistos({&timeHistosSt345[6], 6}, "St345_DigitsAtClusterPos");
  DrawROFTimeHistos({&rofTimeHistos[0], 6});
  DrawChargeHistos({&chargeHistos[0], 6}, "AllDigits");
  DrawChargeHistos({&chargeHistos[6], 6}, "DigitsAtClusterPos");
  DrawChargeHistos({&chargeHistosSt1[0], 6}, "St1_AllDigits");
  DrawChargeHistos({&chargeHistosSt1[6], 6}, "St1_DigitsAtClusterPos");
  DrawChargeHistos({&chargeHistosSt2[0], 6}, "St2_AllDigits");
  DrawChargeHistos({&chargeHistosSt2[6], 6}, "St2_DigitsAtClusterPos");
  DrawChargeHistos({&chargeHistosSt345[0], 6}, "St345_AllDigits");
  DrawChargeHistos({&chargeHistosSt345[6], 6}, "St345_DigitsAtClusterPos");
  DrawCorrelationHistos(corrHistos);
  DrawMultHistos(multHistos);
  TCanvas* cMass = new TCanvas();
  gPad->SetLogy();
  hMass->Draw();
}

//_________________________________________________________________________________________________
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName)
{
  /// open the input file and get the intput tree

  TFile* f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    LOG(error) << "opening file " << fileName << " failed";
    exit(-1);
  }

  TTreeReader* r = new TTreeReader(treeName, f);
  if (r->IsZombie()) {
    LOG(error) << "tree " << treeName << " not found";
    exit(-1);
  }

  return std::make_tuple(f, r);
}

//_________________________________________________________________________________________________
void LoadDigits(TrackInfo& trackInfo, const std::vector<mch::Cluster>& clusters, const std::vector<mch::Digit>& digits,
                bool selectSignal, bool rejectBackground)
{
  /// fill the lists of digits associated to the track

  int nClusterOnTopOfNoDigit(0);

  for (int iCl = trackInfo.mchTrack.getFirstClusterIdx(); iCl <= trackInfo.mchTrack.getLastClusterIdx(); ++iCl) {

    const auto& cluster = clusters[iCl];

    // get the pads at the cluster position
    math_utils::Point3D<float> global{cluster.x, cluster.y, cluster.z};
    auto t = transformation(cluster.getDEId());
    auto local = t^(global);
    int padIDNB(-1), padIDB(-1);
    auto& segmentation = mch::mapping::segmentation(cluster.getDEId());
    bool padsFound = segmentation.findPadPairByPosition(local.x(), local.y(), padIDB, padIDNB);
    bool padFoundNB = padsFound || segmentation.isValid(padIDNB);
    bool padFoundB = padsFound || segmentation.isValid(padIDB);
    if (!padFoundNB && !padFoundB) {
      LOG(warning) << "cluster on top of no pad";
    }

    bool digitFound(false);
    for (uint32_t iDig = 0; iDig < cluster.nDigits; ++iDig) {
      const auto& digit = digits[cluster.firstDigit + iDig];
      if (selectSignal) {
        double nSample = digit.getNofSamples();
        if (digit.getNofSamples() < minNSamplesSignal || digit.getADC() < signalCut(&nSample, signalParam)) {
          continue;
        }
      }
      if (rejectBackground) {
        double nSample = digit.getNofSamples();
        if (digit.getNofSamples() < minNSamplesBackground || digit.getADC() < backgroundCut(&nSample, backgroundParam)) {
          continue;
        }
      }
      trackInfo.digits.push_back(&digit);
      if ((padFoundNB && digit.getPadID() == padIDNB) || (padFoundB && digit.getPadID() == padIDB)) {
        trackInfo.digitsAtClusterPos.push_back(&digit);
        digitFound = true;
      }
    }
    if (!digitFound) {
      ++nClusterOnTopOfNoDigit;
    }
  }

  if (nClusterOnTopOfNoDigit > 0 && trackInfo.midTime >= 0) {
    LOG(warning) << "matched track with " << nClusterOnTopOfNoDigit << "/"
                 << trackInfo.mchTrack.getNClusters() << " clusters on top of no digit";
  }
}

//_________________________________________________________________________________________________
uint32_t IRtoBCinTF(InteractionRecord ir, uint32_t firstTForbit0, bool isMC)
{
  /// convert IR to BC from TF start
  if (ir.orbit < firstTForbit0) {
    LOGP(error, "IR orbit {} < first orbit of first TF {} !?", ir.orbit, firstTForbit0);
  }
  uint32_t orbitInTF = ir.orbit - firstTForbit0;
  if (!isMC) {
    orbitInTF = orbitInTF % nOrbitsPerTF;
  }
  return orbitInTF * constants::lhc::LHCMaxBunches + ir.bc;
}

//_________________________________________________________________________________________________
void computeMCHTime(const std::vector<const mch::Digit*>& digits, double& mean, double& rms, int deMin, int deMax)
{
  /// compute the average time and time dispersion of MCH digits

  if (digits.empty()) {
    LOG(error) << "cannot compute mch time";
    return;
  }

  mean = 0.;
  double t2 = 0.;
  double n = 0.;
  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      mean += digit->getTime();
      t2 += static_cast<double>(digit->getTime()) * digit->getTime();
      n += 1.;
    }
  }

  if (n == 0.) {
    LOG(error) << "cannot compute mch time";
    mean = -1.;
    return;
  }

  mean /= n;
  t2 /= n;
  rms = sqrt(t2 - mean * mean);
  // add 1.5 BC to account for the fact that the actual digit time in BC units
  // can be between t and t+3, hence t+1.5 in average
  mean += 1.5;
}

//_________________________________________________________________________________________________
const dataformats::TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<dataformats::TrackMCHMID>& muonTracks)
{
  /// find the MCH-MID matched track corresponding to this MCH track
  for (const auto& muon : muonTracks) {
    if (muon.getMCHRef().getIndex() == iMCHTrack) {
      return &muon;
    }
  }
  return nullptr;
}

//_________________________________________________________________________________________________
bool ExtrapToVertex(TrackInfo& trackInfo)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber
  /// return false if the propagation fails

  // extrapolate to vertex
  trackInfo.paramAtVertex.setZ(trackInfo.mchTrack.getZ());
  trackInfo.paramAtVertex.setParameters(trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToVertex(trackInfo.paramAtVertex, 0., 0., 0., 0., 0.)) {
    return false;
  }

  // extrapolate to DCA
  mch::TrackParam trackParamAtDCA(trackInfo.mchTrack.getZ(), trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, 0.)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor();
  double dcaY = trackParamAtDCA.getBendingCoor();
  trackInfo.dca = sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  mch::TrackParam trackParamAtRAbs(trackInfo.mchTrack.getZ(), trackInfo.mchTrack.getParameters());
  if (!mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
    return false;
  }
  double xAbs = trackParamAtRAbs.getNonBendingCoor();
  double yAbs = trackParamAtRAbs.getBendingCoor();
  trackInfo.rAbs = sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(TrackInfo& trackInfo)
{
  /// apply standard track selections + pDCA

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(trackInfo.rAbs / 505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double p = trackInfo.paramAtVertex.p();
  double eta = 0.5 * log((p + trackInfo.paramAtVertex.pz()) / (p - trackInfo.paramAtVertex.pz()));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = trackInfo.mchTrack.getP() * trackInfo.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double nrp = nSigmaPDCA * relPRes * p;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * p;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  // double pT = sqrt(trackInfo.paramAtVertex.px() * trackInfo.paramAtVertex.px() +
  //                  trackInfo.paramAtVertex.py() * trackInfo.paramAtVertex.py());
  // if (pT < 1.) {
  //   return false;
  // }

  return true;
}

//_________________________________________________________________________________________________
bool IsSignal(TrackInfo& trackInfo)
{
  /// check if the track has still enough digits in time to pass the signal selection

  int nDigits(0);

  for (const auto digit : trackInfo.digits) {
    if (digit->getTime() + 1.5 >= trackInfo.mchTime - bcIntegrationRange && digit->getTime() + 1.5 <= trackInfo.mchTime + bcIntegrationRange) {
      ++nDigits;
    }
  }

  return nDigits > minNDigitsSignal;
}

//_________________________________________________________________________________________________
bool IsReconstructible(TrackInfo& trackInfo)
{
  /// check if the track has still enough digits to be reconstructible

  bool hasDigits[10] = {false, false, false, false, false, false, false, false, false, false};
  for (const auto digit : trackInfo.digits) {
    hasDigits[digit->getDetID() / 100 - 1] = true;
  }

  int nFiredChambersSt45 = 0;
  for (int i = 6; i < 10; ++i) {
    if (hasDigits[i]) {
      ++nFiredChambersSt45;
    }
  }

  return (hasDigits[0] || hasDigits[1]) && (hasDigits[2] || hasDigits[3]) && (hasDigits[4] || hasDigits[5]) && nFiredChambersSt45 >= 3;
}

//_________________________________________________________________________________________________
void CreateHistosAtVertex(std::vector<TH1*>& histos, const char* extension)
{
  /// create single muon histograms at vertex

  histos.emplace_back(new TH1F(Form("pT%s", extension), "pT;p_{T} (GeV/c)", 300, 0., 30.));
  histos.emplace_back(new TH1F(Form("eta%s", extension), "eta;eta", 200, -4.5, -2.));
  histos.emplace_back(new TH1F(Form("phi%s", extension), "phi;phi", 360, 0., 360.));
  histos.emplace_back(new TH1F(Form("dca%s", extension), "DCA;DCA (cm)", 500, 0., 500.));
  histos.emplace_back(new TH1F(Form("pDCA23%s", extension), "pDCA for #theta_{abs} < 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("pDCA310%s", extension), "pDCA for #theta_{abs} > 3#circ;pDCA (GeV.cm/c)", 2500, 0., 5000.));
  histos.emplace_back(new TH1F(Form("rAbs%s", extension), "rAbs;R_{abs} (cm)", 1000, 0., 100.));
  histos.emplace_back(new TH1F(Form("nClusters%s", extension), "number of clusters per track;n_{clusters}", 20, 0., 20.));
  histos.emplace_back(new TH1F(Form("chi2%s", extension), "normalized #chi^{2};#chi^{2} / ndf", 500, 0., 50.));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillHistosAtVertex(const TrackInfo& trackInfo, std::vector<TH1*>& histos)
{
  /// fill single muon histograms at vertex

  double thetaAbs = TMath::ATan(trackInfo.rAbs / 505.) * TMath::RadToDeg();
  double pDCA = trackInfo.mchTrack.getP() * trackInfo.dca;
  double pT = sqrt(trackInfo.paramAtVertex.px() * trackInfo.paramAtVertex.px() +
                   trackInfo.paramAtVertex.py() * trackInfo.paramAtVertex.py());
  double p = trackInfo.paramAtVertex.p();
  double eta = 0.5 * log((p + trackInfo.paramAtVertex.pz()) / (p - trackInfo.paramAtVertex.pz()));
  double phi = 180. + atan2(-trackInfo.paramAtVertex.px(), -trackInfo.paramAtVertex.py()) / pi() * 180.;

  histos[0]->Fill(pT);
  histos[1]->Fill(eta);
  histos[2]->Fill(phi);
  histos[3]->Fill(trackInfo.dca);
  if (thetaAbs < 3) {
    histos[4]->Fill(pDCA);
  } else {
    histos[5]->Fill(pDCA);
  }
  histos[6]->Fill(trackInfo.rAbs);
  histos[7]->Fill(trackInfo.mchTrack.getNClusters());
  histos[8]->Fill(trackInfo.mchTrack.getChi2OverNDF());
}

//_________________________________________________________________________________________________
void DrawHistosAtVertex(std::vector<TH1*> histos[2])
{
  /// draw histograms at vertex

  // find the optimal number of pads
  int nPadsx(1), nPadsy(1);
  while ((int)histos[0].size() > nPadsx * nPadsy) {
    if (nPadsx == nPadsy) {
      ++nPadsx;
    } else {
      ++nPadsy;
    }
  }

  // draw histograms
  TCanvas* cHist = new TCanvas("histos", "histos", 10, 10, TMath::Max(nPadsx * 300, 1200), TMath::Max(nPadsy * 300, 900));
  cHist->Divide(nPadsx, nPadsy);
  for (int i = 0; i < (int)histos[0].size(); ++i) {
    cHist->cd((i / nPadsx) * nPadsx + i % nPadsx + 1);
    gPad->SetLogy();
    histos[0][i]->SetStats(false);
    histos[0][i]->SetLineColor(4);
    histos[0][i]->Draw();
    histos[1][i]->SetLineColor(2);
    histos[1][i]->Draw("same");
  }

  // add a legend
  TLegend* lHist = new TLegend(0.5, 0.65, 0.9, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0][0], Form("%g mch tracks", histos[0][0]->GetEntries()), "l");
  lHist->AddEntry(histos[1][0], Form("%g muon tracks", histos[1][0]->GetEntries()), "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateTimeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create track time histograms

  histos.emplace_back(new TH1F(Form("timeResVsMCH%s", extension), "#Deltat vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeResVsMCH%sMatch", extension), "#Deltat vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeResVsMID%s", extension), "#Deltat vs MID time (matched tracks);#Deltat (BC)", 4001, -2000.5, 2000.5));
  histos.emplace_back(new TH1F(Form("timeRMS%s", extension), "MCH time dispersion;#sigmat (BC)", 4001, -0.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeRMS%sMatch", extension), "MCH time dispersion;#sigmat (BC)", 4001, -0.25, 2000.25));
  histos.emplace_back(new TH1F(Form("timeDiffMCHMID%s", extension), "<MCH time> - MID time (matched tracks);#Deltat (BC)", 8001, -2000.25, 2000.25));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillTimeHistos(const std::vector<const mch::Digit*>& digits, double mchTime, double mchTimeRMS, int midTime,
                    gsl::span<TH1*> histos, int deMin, int deMax)
{
  /// fill track time histograms

  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      histos[0]->Fill(digit->getTime() + 1.5 - mchTime);
    }
  }
  histos[3]->Fill(mchTimeRMS);

  if (midTime >= 0) {
    for (const auto digit : digits) {
      if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
        histos[1]->Fill(digit->getTime() + 1.5 - mchTime);
        histos[2]->Fill(digit->getTime() + 1.5 - midTime);
      }
    }
    histos[4]->Fill(mchTimeRMS);
    histos[5]->Fill(mchTime - midTime);
  }
}

//_________________________________________________________________________________________________
void DrawTimeHistos(gsl::span<TH1*> histos, const char* extension)
{
  /// draw track time histograms

  TCanvas* cHist = new TCanvas(Form("cTime%s", extension), Form("cTime%s", extension), 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[1]->SetLineColor(2);
  histos[1]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(2);
  histos[2]->Draw();
  cHist->cd(3);
  gPad->SetLogy();
  histos[3]->SetStats(false);
  histos[3]->SetLineColor(4);
  histos[3]->Draw();
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  cHist->cd(4);
  gPad->SetLogy();
  histos[5]->SetStats(false);
  histos[5]->SetLineColor(2);
  histos[5]->Draw();

  TLegend* lHist = new TLegend(0.1, 0.8, 0.5, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0], "all tracks", "l");
  lHist->AddEntry(histos[1], "matched tracks", "l");
  cHist->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void CreateROFTimeHistos(std::vector<TH1*>& histos)
{
  /// create histograms to compare the track time with ROF time range

  histos.emplace_back(new TH1F("rofMinVsMCHTime", "ROF boundaries vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F("rofMaxVsMCHTime", "ROF boundaries vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F("rofMinVsMIDTime", "ROF boundaries vs MID time;#Deltat (BC)", 4001, -2000.5, 2000.5));
  histos.emplace_back(new TH1F("rofMaxVsMIDTime", "ROF boundaries vs MID time;#Deltat (BC)", 4001, -2000.5, 2000.5));
  histos.emplace_back(new TH1F("rofCenterVsMCHTime", "ROF center vs <MCH time>;#Deltat (BC)", 8001, -2000.25, 2000.25));
  histos.emplace_back(new TH1F("rofCenterVsMIDTime", "ROF center vs MID time;#Deltat (BC)", 8001, -2000.25, 2000.25));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}
//_________________________________________________________________________________________________
void FillROFTimeHistos(int rofTime[2], double mchTime, int midTime, gsl::span<TH1*> histos)
{
  /// fill histograms to compare the track time with ROF time range

  double center = 0.5 * (rofTime[1] + rofTime[0]);

  histos[0]->Fill(rofTime[0] - mchTime);
  histos[1]->Fill(rofTime[1] - mchTime);
  histos[4]->Fill(center - mchTime);

  if (midTime >= 0) {
    histos[2]->Fill(rofTime[0] - midTime);
    histos[3]->Fill(rofTime[1] - midTime);
    histos[5]->Fill(center - midTime);
  }
}

//_________________________________________________________________________________________________
void DrawROFTimeHistos(gsl::span<TH1*> histos)
{
  /// draw histograms to compare the track time with ROF time range
  TCanvas* cHist = new TCanvas("cROFTime", "cROFTime", 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(3);
  histos[0]->Draw();
  histos[1]->SetLineColor(6);
  histos[1]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(3);
  histos[2]->Draw();
  histos[3]->SetLineColor(6);
  histos[3]->Draw("same");
  cHist->cd(3);
  gPad->SetLogy();
  histos[4]->SetStats(false);
  histos[4]->Draw();
  cHist->cd(4);
  gPad->SetLogy();
  histos[5]->SetStats(false);
  histos[5]->Draw();
}

//_________________________________________________________________________________________________
void CreateChargeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create track charge histograms

  histos.emplace_back(new TH1F(Form("ADC%s", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%s", extension), "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSample%s", extension), "ADC vs N samples (all tracks);N samples;ADC", 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  histos.emplace_back(new TH1F(Form("ADC%sMatch", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%sMatch", extension), "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSample%sMatch", extension), "ADC vs N samples (matched tracks);N samples;ADC", 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillChargeHistos(const std::vector<const mch::Digit*>& digits, gsl::span<TH1*> histos, int deMin, int deMax)
{
  /// fill track charge histograms
  for (const auto digit : digits) {
    if (digit->getDetID() >= deMin && digit->getDetID() <= deMax) {
      histos[0]->Fill(digit->getADC());
      histos[1]->Fill(digit->getNofSamples());
      histos[2]->Fill(digit->getNofSamples(), digit->getADC());
    }
  }
}

//_________________________________________________________________________________________________
void DrawChargeHistos(gsl::span<TH1*> histos, const char* extension)
{
  /// draw track charge histograms

  TCanvas* cHist = new TCanvas(Form("cCharge%s", extension), Form("cCharge%s", extension), 10, 10, 800, 800);
  cHist->Divide(2, 2);
  cHist->cd(1);
  gPad->SetLogy();
  histos[0]->SetStats(false);
  histos[0]->SetLineColor(4);
  histos[0]->Draw();
  histos[3]->SetLineColor(2);
  histos[3]->Draw("same");
  cHist->cd(2);
  gPad->SetLogy();
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(4);
  histos[1]->Draw();
  histos[4]->SetLineColor(2);
  histos[4]->Draw("same");
  cHist->cd(3);
  gPad->SetLogz();
  histos[2]->SetStats(false);
  histos[2]->Draw("colz");
  cHist->cd(4);
  gPad->SetLogz();
  histos[5]->SetStats(false);
  histos[5]->Draw("colz");

  TLegend* lHist = new TLegend(0.5, 0.8, 0.9, 0.9);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos[0], "all tracks", "l");
  lHist->AddEntry(histos[3], "matched tracks", "l");
  cHist->cd(1);
  lHist->Draw("same");

  static TF1* fSignal = new TF1("fSignal", signalCut, 0, 1023, 4);
  fSignal->SetParameters(signalParam);
  fSignal->SetLineColor(2);
  static TF1* fBackground = new TF1("fBackground", backgroundCut, 0, 1023, 4);
  fBackground->SetParameters(backgroundParam);
  fBackground->SetLineColor(4);
  cHist->cd(3);
  fSignal->Draw("same");
  fBackground->Draw("same");
  cHist->cd(4);
  fSignal->Draw("same");
  fBackground->Draw("same");
}

//_________________________________________________________________________________________________
void CreateCorrelationHistos(std::vector<TH1*>& histos)
{
  /// create correlation histograms between number of digits and total charge

  histos.emplace_back(new TH2F("ChargevsNDigits", "Charge vs N digits;N digits;ADC", 100, 0, 100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt1", "Charge vs N digits (St1);N digits;ADC", 100, 0, 100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt2", "Charge vs N digits (St2);N digits;ADC", 100, 0, 100, 10000, 0, 100000));
  histos.emplace_back(new TH2F("ChargevsNDigitsSt345", "Charge vs N digits (St345);N digits;ADC", 100, 0, 100, 10000, 0, 100000));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillCorrelationHistos(const std::vector<const mch::Digit*>& digits, TH1* hist, double timeRef, int deMin, int deMax)
{
  /// fill correlation histograms between number of digits and total charge

  uint32_t charge(0);
  int nDigits(0);

  for (const auto digit : digits) {
    if (digit->getDetID() < deMin || digit->getDetID() > deMax) {
      continue;
    }
    if (digit->getTime() + 1.5 < timeRef - bcIntegrationRange || digit->getTime() + 1.5 > timeRef + bcIntegrationRange) {
      continue;
    }
    charge += digit->getADC();
    ++nDigits;
  }

  hist->Fill(nDigits, charge);
}

//_________________________________________________________________________________________________
void DrawCorrelationHistos(std::vector<TH1*>& histos)
{
  /// draw correlation histograms between number of digits and total charge

  TCanvas* cCorr = new TCanvas("cCorr", "cCorr", 10, 10, 800, 800);
  cCorr->Divide(2, 2);
  for (int i = 0; i < 4; ++i) {
    cCorr->cd(i + 1);
    gPad->SetLogz();
    histos[i]->Draw("boxcolz");
  }
}

//_________________________________________________________________________________________________
void CreateMultHistos(std::vector<TH1*>& histos)
{
  /// create multiplicity histograms

  static const int nDE[10] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};

  histos.emplace_back(new TH1D("nClustersTotPerCh", "number of clusters per chamber;chamber", 10, 0.5, 10.5));

  for (int i = 0; i < 10; ++i) {
    histos.emplace_back(new TH1D(fmt::format("nClustersTotPerDEinCh{}", i + 1).c_str(),
                                 fmt::format("number of clusters per DE in chamber {};DE", i + 1).c_str(),
                                 nDE[i], 100 * (i + 1) - 0.5, 100 * (i + 1) + nDE[i] - 0.5));
  }

  histos.emplace_back(new TH1D("nClustersTot", "number of clusters;#clusters", 10001, -0.5, 10000.5));

  for (int i = 0; i < 10; ++i) {
    histos.emplace_back(new TH1D(fmt::format("nClustersTotinCh{}", i + 1).c_str(),
                                 fmt::format("number of clusters in chamber {};#clusters", i + 1).c_str(),
                                 2001, -0.5, 2000.5));
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillMultHistos(const mch::ROFRecord& rof, const std::vector<mch::Cluster>& clusters, std::vector<TH1*>& histos)
{
  /// fill multiplicity histograms

  std::array<int, 10> nClusters{};

  for (int i = rof.getFirstIdx(); i <= rof.getLastIdx(); ++i) {

    const auto& cluster = clusters[i];
    int chId = cluster.getChamberId();

    histos[0]->Fill(chId + 1);
    histos[chId + 1]->Fill(cluster.getDEId());

    nClusters[chId]++;
  }

  histos[11]->Fill(rof.getNEntries());

  for (int i = 0; i < 10; ++i) {
    histos[12 + i]->Fill(nClusters[i]);
  }
}

//_________________________________________________________________________________________________
void DrawMultHistos(std::vector<TH1*>& histos)
{
  /// draw multiplicity histograms

  TCanvas* cN = new TCanvas("cNClustersSelected", "cNClustersSelected", 10, 10, 300, 300);
  cN->cd();
  histos[11]->Draw();
  gPad->SetLogy();

  TCanvas* cNCh = new TCanvas("cNClustersPerChDistSelected", "cNClustersPerChDistSelected", 10, 10, 1500, 600);
  cNCh->Divide(5, 2);
  for (int i = 1; i <= 10; ++i) {
    cNCh->cd(i);
    histos[11 + i]->Draw();
    gPad->SetLogy();
  }
}

//_________________________________________________________________________________________________
double signalCut(double* x, double* p)
{
  /// function used to select the signal
  double x0 = pow(p[0] / p[2], 1. / p[3]) + p[1];
  if (x[0] < x0) {
    return p[0];
  } else {
    return p[2] * pow(x[0] - p[1], p[3]);
  }
}

//_________________________________________________________________________________________________
double backgroundCut(double* x, double* p)
{
  /// function used to select the signal
  double x0 = (p[3] * p[2] - p[1] * p[0]) / (p[3] - p[1]);
  if (x[0] < x0) {
    return p[1] * (x[0] - p[0]);
  } else {
    return p[3] * (x[0] - p[2]);
  }
}

//_________________________________________________________________________________________________
void WriteHistos(TFile* f, const char* dirName, const std::vector<TH1*>& histos)
{
  /// write histograms in the subdirectory dirName

  f->mkdir(dirName, dirName, true);
  f->cd(dirName);

  for (auto h : histos) {
    h->Write();
  }
}
