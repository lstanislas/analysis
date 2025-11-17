#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TDatabasePDG.h>
#include <Math/Vector4D.h>

#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/Cluster.h"
#include "MCHBase/TrackBlock.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHTracking/TrackExtrap.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

#include "/Users/PILLOT/Work/Alice/Macros/Tracking/TrackMCHv1.h"

using namespace std;
using namespace ROOT::Math;
using o2::dataformats::TrackMCHMID;
using o2::mch::Cluster;
using o2::mch::ROFRecord;
using o2::mch::Track;
using o2::mch::TrackExtrap;
using o2::mch::TrackFitter;
using o2::mch::TrackMCH;
using o2::mch::TrackMCHv1;
using o2::mch::TrackParam;
using o2::mch::TrackParamStruct;

static const double muMass = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
double chi2Max = 2. * 4. * 4.;
TrackFitter trackFitter{};

//_________________________________________________________________________________________________
struct TrackAtVtxStruct {
  TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

//_________________________________________________________________________________________________
struct TrackStruct {
  PxPyPzMVector pxpypzm{};
  short sign = 0;
  double dca = 0.;
  double rAbs = 0.;
  double chi2 = 0.;
  TrackParamStruct param{};
  std::vector<Cluster> clusters{};
  Track track{};
  bool matchFound = false;
  bool matchIdentical = false;

  //_______________________________________________________________________________________________
  bool operator==(const TrackStruct& track) const
  {
    /// tracks are considered identical when all their clusters match within chi2Max
    if (this->clusters.size() != track.clusters.size()) {
      return false;
    }
    for (size_t iCl = 0; iCl != this->clusters.size(); ++iCl) {
      auto& cl1 = this->clusters[iCl];
      auto& cl2 = track.clusters[iCl];
      if (cl1.getDEId() != cl2.getDEId()) {
        return false;
      }
      double dx = cl1.getX() - cl2.getX();
      double dy = cl1.getY() - cl2.getY();
      double chi2 = dx * dx / (cl1.getEx2() + cl2.getEx2()) + dy * dy / (cl1.getEy2() + cl2.getEy2());
      if (chi2 > chi2Max) {
        return false;
      }
    }
    return true;
  }

  //_______________________________________________________________________________________________
  bool match(const TrackStruct& track) const
  {
    /// Try to match this track with the given track. Matching conditions:
    /// - more than 50% of clusters from one of the two tracks matched with clusters from the other
    /// - at least 1 cluster matched before and 1 cluster matched after the dipole

    size_t nMatchClusters(0);
    bool matchCluster[10] = {false, false, false, false, false, false, false, false, false, false};

    for (const auto& cl1 : this->clusters) {
      for(const auto& cl2 : track.clusters) {
        if (cl1.getDEId() == cl2.getDEId()) {
          double dx = cl1.getX() - cl2.getX();
          double dy = cl1.getY() - cl2.getY();
          double chi2 = dx * dx / (cl1.getEx2() + cl2.getEx2()) + dy * dy / (cl1.getEy2() + cl2.getEy2());
          if (chi2 <= chi2Max) {
            matchCluster[cl1.getChamberId()] = true;
            ++nMatchClusters;
            break;
          }
        }
      }
    }

    return ((matchCluster[0] || matchCluster[1] || matchCluster[2] || matchCluster[3]) &&
            (matchCluster[6] || matchCluster[7] || matchCluster[8] || matchCluster[9]) &&
            (2 * nMatchClusters > this->clusters.size() || 2 * nMatchClusters > track.clusters.size()));
  }
};

//_________________________________________________________________________________________________
void LoadCCDB(int run);
int ReadNextEvent(ifstream& inFile, bool selectTracks, std::vector<TrackStruct>& tracks);
bool ReadTrack(ifstream& inFile, TrackStruct& track);
template <class T>
void ReadNextEventV5(ifstream& inFile, int& event, bool selectTracks, std::vector<TrackStruct>& tracks);
std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
template <typename T>
bool FillTrack(TrackStruct& track, const TrackAtVtxStruct* trackAtVtx, const T& mchTrack, std::vector<Cluster>& clusters);
const TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<TrackMCHMID>& muonTracks);
bool ExtrapToVertex(TrackStruct& track);
bool IsSelected(TrackStruct& track);
void CompareTracks(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, std::vector<TH1*>& histos);
void CreateResiduals(std::vector<TH1*>& histos, const char* extension, double range);
void FillResiduals(std::vector<TrackStruct>& tracks, std::vector<TH1*>& histos, bool rejectMonoCathodes, bool matched = false);
void FillResiduals(TrackStruct& track1, TrackStruct& track2, std::vector<TH1*>& histos);
void NormResiduals(std::vector<TH1*>& histos);
void DrawResiduals(std::vector<TH1*>& histos, const char* extension);
void DrawResiduals(std::vector<TH1*>& histos1, std::vector<TH1*>& histos2, const char* extension);
void DrawRatios(std::vector<TH1*>& histos1, std::vector<TH1*>& histos2, const char* extension);
pair<double, double> GetSigma(TH1* h, int color);
double CrystalBallSymmetric(double* xx, double* par);

//_________________________________________________________________________________________________
void CompareTrackResolution(int run, float l3Current, float dipoleCurrent,
                            string inFileName1, int versionFile1,
                            string inFileName2, int versionFile2,
                            bool selectTracks = false,
                            bool rejectMonoCathodes = true)
{
  /// Compare the cluster-track residuals between the tracks stored in the 2 binary files
  /// file version 4: param at vertex + dca + rAbs + chi2 + param at 1st cluster + clusters (v2)
  /// file version 5: nTracksAtVtx + nMCHTracks + nClusters + list of Tracks at vertex (param at vertex + dca + rAbs + mchTrackIdx) + list of MCH tracks (v1) + list of clusters (v2)
  /// file version 6: nTracksAtVtx + nMCHTracks + nClusters + list of Tracks at vertex (param at vertex + dca + rAbs + mchTrackIdx) + list of MCH tracks (v2) + list of clusters (v2)

  // open files
  ifstream inFile1(inFileName1, ios::binary);
  ifstream inFile2(inFileName2, ios::binary);
  if (!inFile1.is_open() || !inFile2.is_open()) {
    cout << "fail opening files" << endl;
    return;
  }

  // prepare the track fitter
  if (run > 0) {
    LoadCCDB(run);
  } else {
    trackFitter.initField(l3Current, dipoleCurrent);
  }
  trackFitter.smoothTracks(true);

  int event1(-1), event2(-1);
  std::vector<TrackStruct> tracks1{};
  std::vector<TrackStruct> tracks2{};
  bool readNextEvent1(true), readNextEvent2(true);
  std::vector<TH1*> residuals[5] = {{}, {}, {}, {}, {}};
  CreateResiduals(residuals[0], "All1", 2.);
  CreateResiduals(residuals[1], "All2", 2.);
  CreateResiduals(residuals[2], "Matched1", 2.);
  CreateResiduals(residuals[3], "Matched2", 2.);
  CreateResiduals(residuals[4], "ClCl", 0.2);

  while (true) {

    cout << "\rprocessing event " << event1 + 1 << "..." << flush;

    if (readNextEvent1) {
      if (versionFile1 < 5) {
        event1 = ReadNextEvent(inFile1, selectTracks, tracks1);
      } else if (versionFile1 == 5) {
        ReadNextEventV5<TrackMCHv1>(inFile1, event1, selectTracks, tracks1);
      } else {
        ReadNextEventV5<TrackMCH>(inFile1, event1, selectTracks, tracks1);
      }
      FillResiduals(tracks1, residuals[0], rejectMonoCathodes);
    }

    if (readNextEvent2) {
      if (versionFile2 < 5) {
        event2 = ReadNextEvent(inFile2, selectTracks, tracks2);
      } else if (versionFile2 == 5) {
        ReadNextEventV5<TrackMCHv1>(inFile2, event2, selectTracks, tracks2);
      } else {
        ReadNextEventV5<TrackMCH>(inFile2, event2, selectTracks, tracks2);
      }
      FillResiduals(tracks2, residuals[1], rejectMonoCathodes);
    }

    // reaching end of both files
    if (event1 < 0 && event2 < 0) {
      break;
    }

    if (event1 == event2) {
      // reading the same event --> we can compare tracks
      CompareTracks(tracks1, tracks2, residuals[4]);
      FillResiduals(tracks1, residuals[2], rejectMonoCathodes, true);
      FillResiduals(tracks2, residuals[3], rejectMonoCathodes, true);
      readNextEvent1 = true;
      readNextEvent2 = true;
    } else if (event2 < 0 || (event1 >= 0 && event1 < event2)) {
      // event 2 > event 1 or reaching end of file 2
      readNextEvent1 = true;
      readNextEvent2 = false;
    } else {
      // event 1 > event 2 or reaching end of file 1
      readNextEvent1 = false;
      readNextEvent2 = true;
    }
  }

  cout << "\r\033[Kprocessing completed" << endl;

  gStyle->SetOptStat(1);
  DrawResiduals(residuals[4], "ClCl");
  DrawResiduals(residuals[0], residuals[1], "All");
  DrawResiduals(residuals[2], residuals[3], "Matched");
  DrawRatios(residuals[0], residuals[1], "All");
  DrawRatios(residuals[2], residuals[3], "Matched");

  inFile1.close();
  inFile2.close();
}

//_________________________________________________________________________________________________
void CompareTrackResolution(int run,
                            string mchFileName1, string muonFileName1,
                            string mchFileName2, string muonFileName2,
                            bool selectTracks = false, bool selectMatched = false,
                            bool rejectMonoCathodes = true)
{
  /// Compare the cluster-track residuals between the tracks stored in the 2 sets of root files

  /// access CCDB and prepare track extrapolation
  // o2::base::GeometryManager::loadGeometry("O2geometry.root");
  // trackFitter.initField(29999.998047, 5999.966797);
  LoadCCDB(run);
  trackFitter.smoothTracks(true);

  // load tracks
  auto [fMCH1, mchReader1] = LoadData(mchFileName1.c_str(), "o2sim");
  TTreeReaderValue<std::vector<ROFRecord>> mchROFs1 = {*mchReader1, "trackrofs"};
  TTreeReaderValue<std::vector<TrackMCH>> mchTracks1 = {*mchReader1, "tracks"};
  TTreeReaderValue<std::vector<Cluster>> mchClusters1 = {*mchReader1, "trackclusters"};
  auto [fMUON1, muonReader1] = LoadData(muonFileName1.c_str(), "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> muonTracks1 = {*muonReader1, "tracks"};
  auto [fMCH2, mchReader2] = LoadData(mchFileName2.c_str(), "o2sim");
  TTreeReaderValue<std::vector<ROFRecord>> mchROFs2 = {*mchReader2, "trackrofs"};
  TTreeReaderValue<std::vector<TrackMCH>> mchTracks2 = {*mchReader2, "tracks"};
  TTreeReaderValue<std::vector<Cluster>> mchClusters2 = {*mchReader2, "trackclusters"};
  auto [fMUON2, muonReader2] = LoadData(muonFileName2.c_str(), "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> muonTracks2 = {*muonReader2, "tracks"};
  int nTF1 = mchReader1->GetEntries(false);
  if (muonReader1->GetEntries(false) != nTF1) {
    LOGP(error, " {} and {} do not contain the same number of TF", mchFileName1, muonFileName1);
    exit(-1);
  }
  int nTF2 = mchReader2->GetEntries(false);
  if (muonReader2->GetEntries(false) != nTF2) {
    LOGP(error, " {} and {} do not contain the same number of TF", mchFileName2, muonFileName2);
    exit(-1);
  }
  if (nTF1 != nTF2) {
    LOGP(warning, "the two set of files do not contain the same number of TF ({} vs {})", nTF1, nTF2);
  }

  std::vector<TrackStruct> tracks1{};
  std::vector<TrackStruct> tracks2{};
  std::vector<TH1*> residuals[5] = {{}, {}, {}, {}, {}};
  CreateResiduals(residuals[0], "All1", 2.);
  CreateResiduals(residuals[1], "All2", 2.);
  CreateResiduals(residuals[2], "Matched1", 2.);
  CreateResiduals(residuals[3], "Matched2", 2.);
  CreateResiduals(residuals[4], "ClCl", 0.2);

  int nTF = TMath::Max(nTF1, nTF2);
  for (int iTF = 0; iTF < nTF; ++iTF) {
    cout << "\rprocessing TF " << iTF << "..." << flush;

    size_t nROFs1 = 0;
    if (iTF < nTF1) {
      if (!mchReader1->Next() || !muonReader1->Next()) {
        LOGP(error, "cannot read the next TF in the first set of files");
        exit(-1);
      }
      nROFs1 = mchROFs1->size();
    }

    size_t nROFs2 = 0;
    if (iTF < nTF2) {
      if (!mchReader2->Next() || !muonReader2->Next()) {
        LOGP(error, "cannot read the next TF in the second set of files");
        exit(-1);
      }
      nROFs2 = mchROFs2->size();
    }

    auto nROFs = TMath::Max(nROFs1, nROFs2);
    int iROF1(-1), iROF2(-1);
    for (size_t iROF = 0; iROF < nROFs; ++iROF) {

      if (iROF < nROFs1) {
        ++iROF1;
        const auto& mchROF1 = (*mchROFs1)[iROF];
        tracks1.clear();
        tracks1.reserve(mchROF1.getNEntries());
        for (int iMCHTrack = mchROF1.getFirstIdx(); iMCHTrack <= mchROF1.getLastIdx(); ++iMCHTrack) {
          if (selectMatched && !FindMuon(iMCHTrack, *muonTracks1)) {
            continue;
          }
          tracks1.emplace_back();
          if (!FillTrack(tracks1.back(), nullptr, (*mchTracks1)[iMCHTrack], *mchClusters1)) {
            tracks1.pop_back();
          } else if (selectTracks && !(ExtrapToVertex(tracks1.back()) && IsSelected(tracks1.back()))) {
            tracks1.pop_back();
          }
        }
        FillResiduals(tracks1, residuals[0], rejectMonoCathodes);
      }

      if (iROF < nROFs2) {
        ++iROF2;
        const auto& mchROF2 = (*mchROFs2)[iROF];
        tracks2.clear();
        tracks2.reserve(mchROF2.getNEntries());
        for (int iMCHTrack = mchROF2.getFirstIdx(); iMCHTrack <= mchROF2.getLastIdx(); ++iMCHTrack) {
          if (selectMatched && !FindMuon(iMCHTrack, *muonTracks2)) {
            continue;
          }
          tracks2.emplace_back();
          if (!FillTrack(tracks2.back(), nullptr, (*mchTracks2)[iMCHTrack], *mchClusters2)) {
            tracks2.pop_back();
          } else if (selectTracks && !(ExtrapToVertex(tracks2.back()) && IsSelected(tracks2.back()))) {
            tracks2.pop_back();
          }
        }
        FillResiduals(tracks2, residuals[1], rejectMonoCathodes);
      }

      if (iROF1 == iROF2) {
        // reading the same event --> we can compare tracks
        CompareTracks(tracks1, tracks2, residuals[4]);
        FillResiduals(tracks1, residuals[2], rejectMonoCathodes, true);
        FillResiduals(tracks2, residuals[3], rejectMonoCathodes, true);
      }
    }
  }

  cout << "\r\033[Kprocessing completed" << endl;

  if (nTF1 != nTF2) {
    NormResiduals(residuals[0]);
    NormResiduals(residuals[1]);
  }

  gStyle->SetOptStat(1);
  DrawResiduals(residuals[4], "ClCl");
  DrawResiduals(residuals[0], residuals[1], "All");
  DrawResiduals(residuals[2], residuals[3], "Matched");
  DrawRatios(residuals[0], residuals[1], "All");
  DrawRatios(residuals[2], residuals[3], "Matched");

  fMCH1->Close();
  fMUON1->Close();
  fMCH2->Close();
  fMUON2->Close();
}

//_________________________________________________________________________________________________
void LoadCCDB(int run)
{
  /// access CCDB and prepare track extrapolation

  // load magnetic field and geometry from CCDB
  auto& ccdb = o2::ccdb::BasicCCDBManager::instance();
  auto [tStart, tEnd] = ccdb.getRunDuration(run);
  ccdb.setTimestamp(tEnd);
  auto grp = ccdb.get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  auto geom = ccdb.get<TGeoManager>("GLO/Config/GeometryAligned");

  // prepare track extrapolation
  o2::base::Propagator::initFieldFromGRP(grp);
  TrackExtrap::setField();
}

//_________________________________________________________________________________________________
int ReadNextEvent(ifstream& inFile, bool selectTracks, std::vector<TrackStruct>& tracks)
{
  /// read the next event in the input file

  tracks.clear();
  
  int event(-1);
  if (!inFile.read(reinterpret_cast<char*>(&event), sizeof(int))) {
    // reaching end of file
    return -1;
  }
  
  int size(0);
  inFile.read(reinterpret_cast<char*>(&size), sizeof(int));
  if (size == 0) {
    // empty event
    return event;
  }

  // read the tracks
  int nTracks(0);
  inFile.read(reinterpret_cast<char*>(&nTracks), sizeof(int));
  tracks.reserve(nTracks);
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    tracks.emplace_back();
    if (!ReadTrack(inFile, tracks.back())) {
      tracks.pop_back();
    }
    if (selectTracks && !IsSelected(tracks.back())) {
      tracks.pop_back();
    }
  }

  return event;
}

//_________________________________________________________________________________________________
bool ReadTrack(ifstream& inFile, TrackStruct& track)
{
  /// read one track from the input file and refit it to get the parameters at each cluster
  /// return false if the refitting fails

  TrackParamStruct paramAtVtx{};
  inFile.read(reinterpret_cast<char*>(&paramAtVtx), sizeof(TrackParamStruct));
  track.pxpypzm.SetPx(paramAtVtx.px);
  track.pxpypzm.SetPy(paramAtVtx.py);
  track.pxpypzm.SetPz(paramAtVtx.pz);
  track.pxpypzm.SetM(muMass);
  track.sign = paramAtVtx.sign;
  inFile.read(reinterpret_cast<char*>(&(track.dca)), sizeof(double));
  inFile.read(reinterpret_cast<char*>(&(track.rAbs)), sizeof(double));

  inFile.read(reinterpret_cast<char*>(&(track.param)), sizeof(TrackParamStruct));

  inFile.read(reinterpret_cast<char*>(&(track.chi2)), sizeof(double));

  int nClusters(0);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));
  track.clusters.resize(nClusters);
  inFile.read(reinterpret_cast<char*>(track.clusters.data()), nClusters * sizeof(Cluster));
  for (const auto& cluster : track.clusters) {
    track.track.createParamAtCluster(cluster);
  }

  try {
    TrackExtrap::useExtrapV2();
    trackFitter.fit(track.track);
    return true;
  } catch (exception const& e) {
    return false;
  }
}

//_________________________________________________________________________________________________
template <class T>
void ReadNextEventV5(ifstream& inFile, int& event, bool selectTracks, std::vector<TrackStruct>& tracks)
{
  /// read the next event in the input file

  tracks.clear();

  int nTracksAtVtx(-1);
  if (!inFile.read(reinterpret_cast<char*>(&nTracksAtVtx), sizeof(int))) {
    event = -1;
    return; // reaching end of file ...
  } else {
    ++event; // ... or reading the next event
  }
  int nMCHTracks(-1);
  inFile.read(reinterpret_cast<char*>(&nMCHTracks), sizeof(int));
  int nClusters(-1);
  inFile.read(reinterpret_cast<char*>(&nClusters), sizeof(int));

  // read the tracks at vertex, the MCH tracks and the attached clusters
  std::vector<TrackAtVtxStruct> tracksAtVtx(nTracksAtVtx);
  inFile.read(reinterpret_cast<char*>(tracksAtVtx.data()), nTracksAtVtx * sizeof(TrackAtVtxStruct));
  std::vector<T> mchTracks(nMCHTracks);
  inFile.read(reinterpret_cast<char*>(mchTracks.data()), nMCHTracks * sizeof(T));
  std::vector<Cluster> clusters(nClusters);
  inFile.read(reinterpret_cast<char*>(clusters.data()), nClusters * sizeof(Cluster));

  if (nMCHTracks == 0) {
    return; // need MCH tracks to study the resolution
  }

  if (nTracksAtVtx > 0) {

    // fill the internal track structure based on the provided tracks at vertex
    tracks.reserve(nTracksAtVtx);
    for (const auto& trackAtVtx : tracksAtVtx) {
      tracks.emplace_back();
      if (!FillTrack(tracks.back(), &trackAtVtx, mchTracks[trackAtVtx.mchTrackIdx], clusters)) {
        tracks.pop_back();
      }
      if (selectTracks && !IsSelected(tracks.back())) {
        tracks.pop_back();
      }
    }

  } else {

    // fill the internal track structure based on the MCH tracks
    tracks.reserve(nMCHTracks);
    for (const auto& mchTrack : mchTracks) {
      tracks.emplace_back();
      if (!FillTrack(tracks.back(), nullptr, mchTrack, clusters)) {
        tracks.pop_back();
      }
      if (selectTracks && !(ExtrapToVertex(tracks.back()) && IsSelected(tracks.back()))) {
        tracks.pop_back();
      }
    }
  }
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
template <typename T>
bool FillTrack(TrackStruct& track, const TrackAtVtxStruct* trackAtVtx, const T& mchTrack, std::vector<Cluster>& clusters)
{
  /// fill the internal track structure from the provided informations
  /// return false if the refitting fails

  if (trackAtVtx) {
    track.pxpypzm.SetPx(trackAtVtx->paramAtVertex.px);
    track.pxpypzm.SetPy(trackAtVtx->paramAtVertex.py);
    track.pxpypzm.SetPz(trackAtVtx->paramAtVertex.pz);
    track.pxpypzm.SetM(muMass);
    track.sign = trackAtVtx->paramAtVertex.sign;
    track.dca = trackAtVtx->dca;
    track.rAbs = trackAtVtx->rAbs;
  }

  track.chi2 = mchTrack.getChi2();
  track.param.x = mchTrack.getX();
  track.param.y = mchTrack.getY();
  track.param.z = mchTrack.getZ();
  track.param.px = mchTrack.getPx();
  track.param.py = mchTrack.getPy();
  track.param.pz = mchTrack.getPz();
  track.param.sign = mchTrack.getSign();

  track.clusters.insert(track.clusters.begin(), clusters.begin() + mchTrack.getFirstClusterIdx(),
                        clusters.begin() + mchTrack.getLastClusterIdx() + 1);
  for (const auto& cluster : track.clusters) {
    track.track.createParamAtCluster(cluster);
  }

  try {
    TrackExtrap::useExtrapV2();
    trackFitter.fit(track.track);
    return true;
  } catch (exception const& e) {
    return false;
  }
}

//_________________________________________________________________________________________________
const TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<TrackMCHMID>& muonTracks)
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
bool ExtrapToVertex(TrackStruct& track)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber

  if (!gGeoManager) {
    cout << "Cannot extrapolate to vertex without the geometry" << endl;
    exit(-1);
  }

  // same extrapolation method as in the AOD producer
  TrackExtrap::useExtrapV2(false);

  // convert parameters at first cluster in TrackParam format
  TrackParam trackParam;
  trackParam.setNonBendingCoor(track.param.x);
  trackParam.setBendingCoor(track.param.y);
  trackParam.setZ(track.param.z);
  trackParam.setNonBendingSlope(track.param.px / track.param.pz);
  trackParam.setBendingSlope(track.param.py / track.param.pz);
  trackParam.setInverseBendingMomentum(track.param.sign / TMath::Sqrt(track.param.py * track.param.py + track.param.pz * track.param.pz));

  // extrapolate to vertex
  TrackParam trackParamAtVertex(trackParam);
  if (!TrackExtrap::extrapToVertex(trackParamAtVertex, 0., 0., 0., 0., 0.)) {
    return false;
  }
  track.pxpypzm.SetPx(trackParamAtVertex.px());
  track.pxpypzm.SetPy(trackParamAtVertex.py());
  track.pxpypzm.SetPz(trackParamAtVertex.pz());
  track.pxpypzm.SetM(muMass);
  track.sign = trackParamAtVertex.getCharge();

  // extrapolate to DCA
  TrackParam trackParamAtDCA(trackParam);
  if (!TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, 0.)) {
    return false;
  }
  double dcaX = trackParamAtDCA.getNonBendingCoor();
  double dcaY = trackParamAtDCA.getBendingCoor();
  track.dca = TMath::Sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  if (!TrackExtrap::extrapToZ(trackParam, -505.)) {
    return false;
  }
  double xAbs = trackParam.getNonBendingCoor();
  double yAbs = trackParam.getBendingCoor();
  track.rAbs = TMath::Sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(TrackStruct& track)
{
  /// apply standard track selections + pDCA

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;
/*
  if (track.pxpypzm.Pt() < 1.) {
    return false;
  }
*/
  double thetaAbs = TMath::ATan(track.rAbs/505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double eta = track.pxpypzm.Eta();
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pUncorr = TMath::Sqrt(track.param.px*track.param.px + track.param.py*track.param.py + track.param.pz*track.param.pz);
  double pDCA = pUncorr * track.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double pTot = track.pxpypzm.P();
  double nrp = nSigmaPDCA * relPRes * pTot;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * pTot;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect*pResEffect + slopeResEffect*slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
void CompareTracks(std::vector<TrackStruct>& tracks1, std::vector<TrackStruct>& tracks2, std::vector<TH1*>& histos)
{
  /// compare the tracks between the 2 events

  // first look for identical tracks in the 2 events
  for (auto& track1 : tracks1) {
    // find a track in the second event identical to track1 and not already matched
    auto itTrack2 = tracks2.begin();
    do {
      itTrack2 = find(itTrack2, tracks2.end(), track1);
    } while (itTrack2 != tracks2.end() && itTrack2->matchFound && ++itTrack2 != tracks2.end());
    if (itTrack2 != tracks2.end()) {
      track1.matchFound = true;
      itTrack2->matchFound = true;
      track1.matchIdentical = true;
      itTrack2->matchIdentical = true;
      FillResiduals(track1, *itTrack2, histos);
    }
  }

  // then look for similar tracks in the 2 events
  for (auto& track1 : tracks1) {
    // skip already matched tracks
    if (track1.matchFound) {
      continue;
    }
    // find a track in the second event similar to track1 and not already matched
    for (auto& track2 : tracks2) {
      if (!track2.matchFound && track2.match(track1)) {
        track1.matchFound = true;
        track2.matchFound = true;
        FillResiduals(track1, track2, histos);
        break;
      }
    }
  }
}

//_________________________________________________________________________________________________
void CreateResiduals(std::vector<TH1*>& histos, const char* extension, double range)
{
  /// create histograms of cluster-track residuals
  if (histos.size() == 0) {
    for (int iSt = 1; iSt <= 5; ++iSt) {
      histos.emplace_back(new TH1F(Form("resX%sSt%d", extension, iSt),
                                   Form("#DeltaX Station %d;#DeltaX (cm)", iSt), 2000, -range, range));
      histos.emplace_back(new TH1F(Form("resY%sSt%d", extension, iSt),
                                   Form("#DeltaY Station %d;#DeltaY (cm)", iSt), 2000, -range, range));
    }
    histos.emplace_back(new TH1F(Form("resX%s", extension), "#DeltaX;#DeltaX (cm)", 2000, -range, range));
    histos.emplace_back(new TH1F(Form("resY%s", extension), "#DeltaY;#DeltaY (cm)", 2000, -range, range));
  }
  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillResiduals(TrackStruct& track1, TrackStruct& track2, std::vector<TH1*>& histos)
{
  /// fill histograms of cluster-cluster residuals
  for (const auto& cl1 : track1.clusters) {
    for (const auto& cl2 : track2.clusters) {
      if (cl1.getDEId() == cl2.getDEId()) {
        histos[cl1.getChamberId() / 2 * 2]->Fill(cl2.getX() - cl1.getX());
        histos[cl1.getChamberId() / 2 * 2 + 1]->Fill(cl2.getY() - cl1.getY());
        histos[10]->Fill(cl2.getX() - cl1.getX());
        histos[11]->Fill(cl2.getY() - cl1.getY());
      }
    }
  }
}

//_________________________________________________________________________________________________
void FillResiduals(std::vector<TrackStruct>& tracks, std::vector<TH1*>& histos, bool rejectMonoCathodes, bool matched)
{
  /// fill histograms of cluster-track residuals
  /// if rejectMonoCathodes, skip the residual in the direction that is not measured
  for (const auto& track : tracks) {
    if (!matched || track.matchFound) {
      for (const auto& param : track.track) {
        if (!rejectMonoCathodes || param.getClusterPtr()->getEx() < 5.) {
          histos[param.getClusterPtr()->getChamberId() / 2 * 2]->Fill(param.getClusterPtr()->getX() - param.getNonBendingCoor());
          histos[10]->Fill(param.getClusterPtr()->getX() - param.getNonBendingCoor());
        }
        if (!rejectMonoCathodes || param.getClusterPtr()->getEy() < 5.) {
          histos[param.getClusterPtr()->getChamberId() / 2 * 2 + 1]->Fill(param.getClusterPtr()->getY() - param.getBendingCoor());
          histos[11]->Fill(param.getClusterPtr()->getY() - param.getBendingCoor());
        }
      }
    }
  }
}

//_________________________________________________________________________________________________
void NormResiduals(std::vector<TH1*>& histos)
{
  /// normalize histograms

  for (auto& h : histos) {
    h->Scale(1. / h->GetEntries());
  }
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos, const char* extension)
{
  /// draw cluster-cluster residuals

  int nPadsx = (histos.size() + 1) / 2;
  TCanvas* c = new TCanvas(Form("residual%s", extension), Form("residual%s", extension), 10, 10, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : histos) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    gPad->SetLogy();
    h->Draw();
    h->GetXaxis()->SetRangeUser(-0.02, 0.02);
    ++i;
  }
}

//_________________________________________________________________________________________________
void DrawResiduals(std::vector<TH1*>& histos1, std::vector<TH1*>& histos2, const char* extension)
{
  /// draw cluster-track residuals

  TGraphErrors* g[2][2] = {nullptr};
  const char* dir[2] = {"X", "Y"};
  int color[2] = {4, 2};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      g[i][j] = new TGraphErrors(6);
      g[i][j]->SetName(Form("sigma%s%s%d", dir[i], extension, j));
      g[i][j]->SetTitle(Form("#sigma%s per station;station ID;#sigma%s (cm)", dir[i], dir[i]));
      g[i][j]->SetMarkerStyle(kFullDotLarge);
      g[i][j]->SetMarkerColor(color[j]);
      g[i][j]->SetLineColor(color[j]);
    }
  }

  int nPadsx = (histos1.size() + 1) / 2;
  TCanvas* c = new TCanvas(Form("residual%s", extension), Form("residual%s", extension), 10, 10, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : histos1) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    gPad->SetLogy();
    h->SetLineColor(color[0]);
    h->Draw();
    auto res1 = GetSigma(h, color[0]);
    g[i % 2][0]->SetPoint(i / 2, i / 2 + 1., res1.first);
    g[i % 2][0]->SetPointError(i / 2, 0., res1.second);
    histos2[i]->SetLineColor(color[1]);
    histos2[i]->Draw("sames");
    auto res2 = GetSigma(histos2[i], color[1]);
    g[i % 2][1]->SetPoint(i / 2, i / 2 + 1., res2.first);
    g[i % 2][1]->SetPointError(i / 2, 0., res2.second);
    printf("%s: %f ± %f cm --> %f ± %f cm\n", h->GetName(), res1.first, res1.second, res2.first, res2.second);
    ++i;
  }

  TCanvas* c2 = new TCanvas(Form("sigma%s", extension), Form("sigma%s", extension), 10, 10, 600, 300);
  c2->Divide(2, 1);
  c2->cd(1);
  g[0][0]->Draw("ap");
  g[0][1]->Draw("p");
  c2->cd(2);
  g[1][0]->Draw("ap");
  g[1][1]->Draw("p");

  // add a legend
  TLegend* lHist = new TLegend(0.2, 0.65, 0.4, 0.8);
  lHist->SetFillStyle(0);
  lHist->SetBorderSize(0);
  lHist->AddEntry(histos1[0], "file 1", "l");
  lHist->AddEntry(histos2[0], "file 2", "l");
  c->cd(1);
  lHist->Draw("same");
}

//_________________________________________________________________________________________________
void DrawRatios(std::vector<TH1*>& histos1, std::vector<TH1*>& histos2, const char* extension)
{
  /// draw ratios of cluster-track residuals

  int nPadsx = (histos1.size() + 1) / 2;
  TCanvas* c = new TCanvas(Form("ratio%s", extension), Form("ratio%s", extension), 10, 10, nPadsx * 300, 600);
  c->Divide(nPadsx, 2);
  int i(0);
  for (const auto& h : histos2) {
    c->cd((i % 2) * nPadsx + i / 2 + 1);
    TH1* hRat = new TH1F(*static_cast<TH1F*>(h));
    hRat->SetDirectory(0);
    hRat->Rebin(2);
    auto* h1 = static_cast<TH1F*>(histos1[i]->Clone());
    h1->Rebin(2);
    hRat->Divide(h1);
    delete h1;
    hRat->SetStats(false);
    hRat->SetLineColor(2);
    hRat->Draw();
    hRat->GetXaxis()->SetRangeUser(-0.5, 0.5);
    hRat->GetYaxis()->SetRangeUser(0.8, 1.2);
    ++i;
  }
}

//_________________________________________________________________________________________________
pair<double, double> GetSigma(TH1* h, int color)
{
  /// get the dispersion of the histo

  static TF1* fCrystalBall = new TF1("CrystalBall", CrystalBallSymmetric, -2., 2., 5);
  fCrystalBall->SetLineColor(color);

  if (h->GetEntries() < 10.) {
    return make_pair(0., 0.);
  }

  double sigmaTrk = 0.2; // 2 mm
  double sigmaTrkCut = 4.; // 4 sigma

  // first fit
  double xMin = -0.5 * sigmaTrkCut * sigmaTrk;
  double xMax = 0.5 * sigmaTrkCut * sigmaTrk;
  fCrystalBall->SetRange(xMin, xMax);
  fCrystalBall->SetParameters(h->GetEntries(), 0., 0.1, 2., 1.5);
  fCrystalBall->SetParLimits(1, xMin, xMax);
  fCrystalBall->SetParLimits(2, 0., 1.);
  fCrystalBall->FixParameter(3, 1.e6);
  h->Fit(fCrystalBall, "RNQ");

  // rebin histo
  // int rebin = static_cast<Int_t>(TMath::Min(0.1 * h->GetNbinsX(), TMath::Max(0.3 * fCrystalBall->GetParameter(2) / h->GetBinWidth(1), 1.)));
  // while (h->GetNbinsX() % rebin != 0) {
  //   rebin--;
  // }
  int rebin = 2;
  h->Rebin(rebin);

  // second fit
  fCrystalBall->SetParameter(0, fCrystalBall->GetParameter(0) * rebin);
  fCrystalBall->ReleaseParameter(3);
  fCrystalBall->SetParameter(3, 2.);
  fCrystalBall->SetParameter(4, 1.5);
  h->Fit(fCrystalBall, "RQ");

  return make_pair(fCrystalBall->GetParameter(2), fCrystalBall->GetParError(2));
}

//_________________________________________________________________________________________________
double CrystalBallSymmetric(double* xx, double* par)
{
  /// Crystal Ball definition

  ///par[0] = Normalization
  ///par[1] = mean
  ///par[2] = sigma
  ///par[3] = alpha = alpha'
  ///par[4] = n = n'

  double tp = fabs((xx[0] - par[1]) / par[2]);

  double absAlpha = fabs(par[3]);
  double ap = pow(par[4] / absAlpha, par[4]) * exp(-0.5 * absAlpha * absAlpha);
  double bp = par[4] / absAlpha - absAlpha;

  if (tp < absAlpha) {
    return par[0] * (exp(-0.5 * tp * tp)); // gaussian core
  } else {
    return par[0] * (ap / pow(bp + tp, par[4])); //left and right tails
  }
}
