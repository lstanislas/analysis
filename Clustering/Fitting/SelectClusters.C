#include <cmath>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <gsl/span>

#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
// #include <TRandom.h>

#include "CommonConstants/LHCConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/Digit.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DetectorsBase/GeometryManager.h"
#include "Framework/Logger.h"
#include "MCHBase/TrackBlock.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHTracking/TrackParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"

#include "DataUtils.h"
#include "CCDBUtils.h"
#include "PreClusterUtils.h"
#include "TrackUtils.h"

using o2::dataformats::TrackMCHMID;
using o2::mch::Cluster;
using o2::mch::Digit;
using o2::mch::Track;
using o2::mch::TrackExtrap;
using o2::mch::TrackFitter;
using o2::mch::TrackMCH;
using o2::mch::TrackParam;
using o2::mch::TrackParamStruct;

struct TrackStruct {
  Track track{};
  TrackParam paramAtVtx{};
  double dca = 0.;
  double rAbs = 0.;
};

const TrackMCHMID* FindMuon(uint32_t iMCHTrack, const std::vector<TrackMCHMID>& muonTracks);
bool ClustersToTrack(gsl::span<Cluster> trackClusters, TrackStruct& track);
bool ExtrapToVertex(TrackStruct& track);
bool IsSelected(const TrackStruct& track);
int IRtoBCinTF(o2::InteractionRecord ir, int run, bool isMC);

//_________________________________________________________________________________________________
void SelectClusters(int run, const char* clusterFile, const char* trackFile, const char* muonFile,
                    bool applyTrackSelection = false, bool isMC = false,
                    const char* outFile = "clusters.root")
{
  /// select the isolated clusters attached to a (selected) muon track
  /// store them with the associated track params and digits in outFile
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  /// load CCDB objects and prepare track fitting
  if (run < 300000) {
    o2::base::GeometryManager::loadGeometry("O2geometry.root");
    TrackFitter trackFitter{};
    trackFitter.initField(29999.998047, 5999.966797);
  } else {
    InitFromCCDB(run, false, true, true);
  }

  // load clusters and tracks
  auto [fClusters, clusterReader] = LoadData(clusterFile, "o2sim");
  TTreeReaderValue<std::vector<Cluster>> allClusters(*clusterReader, "clusters");
  TTreeReaderValue<std::vector<Digit>> allClusterDigits(*clusterReader, "clusterdigits");
  auto [fTracks, trackReader] = LoadData(trackFile, "o2sim");
  TTreeReaderValue<std::vector<TrackMCH>> allTracks(*trackReader, "tracks");
  TTreeReaderValue<std::vector<Cluster>> allTrackClusters(*trackReader, "trackclusters");
  auto [fMuons, muonReader] = LoadData(muonFile, "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> allMuons(*muonReader, "tracks");

  // make sure the track clusters point to digits in the cluster file
  if (trackReader->GetTree()->FindBranch("trackdigits")) {
    LOG(error) << " track clusters point to digits in the track file. Re-run the MCH reco without the option --digits";
    exit(-1);
  }

  int nTF = clusterReader->GetEntries(false);
  if (trackReader->GetEntries(false) != nTF ||
      muonReader->GetEntries(false) != nTF) {
    LOG(error) << " not all files contain the same number of TF";
    exit(-1);
  }

  // setup the output
  TFile dataFile(outFile, "recreate");
  TTree dataTree("data", "tree with input data");
  TrackParamStruct trackParamOut;
  dataTree.Branch("trackParameters", &trackParamOut);
  int trackTime;
  dataTree.Branch("trackTime", &trackTime);
  Cluster clusterOut;
  dataTree.Branch("clusters", &clusterOut);
  std::vector<Digit> digitsOut{};
  dataTree.Branch("digits", &digitsOut);

  std::vector<TH1*> preClusterInfo{};
  CreatePreClusterInfo(preClusterInfo);
  std::vector<TH1*> preClusterInfoSt[3] = {{}, {}, {}};
  CreatePreClusterInfo(preClusterInfoSt[0], "St1");
  CreatePreClusterInfo(preClusterInfoSt[1], "St2");
  CreatePreClusterInfo(preClusterInfoSt[2], "St345");

  int iTF(0);
  while (clusterReader->Next() && trackReader->Next() && muonReader->Next()) {
    std::cout << "\rprocessing TF " << ++iTF << " / " << nTF << "..." << std::flush;

    // count the number of clusters associated to each precluster (identified by the first digit index)
    std::unordered_map<uint32_t, int> nClusters{};
    for (const auto& cluster : *allClusters) {
      nClusters[cluster.firstDigit]++;
    }

    for (const auto& muon : *allMuons) {
      const auto& mchTrack = (*allTracks)[muon.getMCHRef().getIndex()];
      gsl::span<Cluster> trackClusters(&(*allTrackClusters)[mchTrack.getFirstClusterIdx()], mchTrack.getNClusters());

      // create an internal track and compute its parameters at each attached cluster
      TrackStruct track{};
      if (!ClustersToTrack(trackClusters, track)) {
        continue;
      }

      // propagate the track to the vertex (0,0,0) and apply selections if requested
      if (applyTrackSelection && !(ExtrapToVertex(track) && IsSelected(track))) {
        continue;
      }

      // get the track time wrt the beginning of the TF
      trackTime = IRtoBCinTF(muon.getIR(), run, isMC);

      for (const auto& param : track.track) {
        const auto& cluster = *(param.getClusterPtr());
        const gsl::span<const Digit> digits(&(*allClusterDigits)[cluster.firstDigit], cluster.nDigits);

        // select clusters of interest (isolated, bi-cathode, ...)
        if (nClusters[cluster.firstDigit] != 1 || IsMonoCathode(digits)) {
          continue;
        }

        // fill precluster characteristics (charge, size, ...)
        const auto [sizeX, sizeY] = GetSize(digits);
        const auto [chargeNB, chargeB] = GetCharge(digits, run < 300000);
        FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfo);
        int iSt = (cluster.getChamberId() < 4) ? cluster.getChamberId() / 2 : 2;
        FillPreClusterInfo(chargeNB, chargeB, sizeX, sizeY, preClusterInfoSt[iSt]);

        // fill output tree
        trackParamOut = param.getTrackParamStruct();
        clusterOut = cluster;
        digitsOut.clear();
        digitsOut.insert(digitsOut.end(), digits.begin(), digits.end());
        dataTree.Fill();
      }
    }
  }
  cout << "\r\033[Kprocessing completed" << endl;

  gStyle->SetOptStat(1);

  auto c = DrawPreClusterInfo(preClusterInfo);
  auto cSt1 = DrawPreClusterInfo(preClusterInfoSt[0], "St1");
  auto cSt2 = DrawPreClusterInfo(preClusterInfoSt[1], "St2");
  auto cSt345 = DrawPreClusterInfo(preClusterInfoSt[2], "St345");

  dataFile.Write();
  c->Write();
  cSt1->Write();
  cSt2->Write();
  cSt345->Write();
  dataFile.Close();
  fClusters->Close();
  fTracks->Close();
  fMuons->Close();
}

//_________________________________________________________________________________________________
void SelectClusters(int run, const char* trackFile, const char* muonFile,
                    bool applyTrackSelection = false, bool selectMatched = false,
                    const char* outClFile = "clusterslite.root", const char* outTrFile = "trackslite.root")
{
  /// select all clusters attached to a (selected) mch/muon track
  /// store the clusters with the associated track params in outClFile
  /// store the track parameters at vertex and dca in outTrFile
  /// require the MCH mapping to be loaded: gSystem->Load("libO2MCHMappingImpl4")

  /// load CCDB objects and prepare track fitting
  if (run < 300000) {
    o2::base::GeometryManager::loadGeometry("O2geometry.root");
    TrackFitter trackFitter{};
    trackFitter.initField(29999.998047, 5999.966797);
  } else {
    InitFromCCDB(run, false, true, true);
  }

  // load clusters and tracks
  auto [fTracks, trackReader] = LoadData(trackFile, "o2sim");
  TTreeReaderValue<std::vector<TrackMCH>> allTracks(*trackReader, "tracks");
  TTreeReaderValue<std::vector<Cluster>> allTrackClusters(*trackReader, "trackclusters");
  auto [fMuons, muonReader] = LoadData(muonFile, "o2sim");
  TTreeReaderValue<std::vector<TrackMCHMID>> allMuons(*muonReader, "tracks");

  int nTF = trackReader->GetEntries(false);
  if (muonReader->GetEntries(false) != nTF) {
    LOG(error) << " not all files contain the same number of TF";
    exit(-1);
  }

  // setup the cluster output
  TFile dataClFile(outClFile, "recreate");
  TTree dataClTree("data", "tree with input data");
  TrackParamStruct trackParamOut;
  dataClTree.Branch("trackParameters", &trackParamOut);
  Cluster clusterOut;
  dataClTree.Branch("clusters", &clusterOut);

  // setup the track output
  TFile dataTrFile(outTrFile, "recreate");
  TTree dataTrTree("data", "tree with input data");
  TrackLite trackOut;
  dataTrTree.Branch("track", &trackOut);

  int iTF(0);
  while (trackReader->Next() && muonReader->Next()) {
    std::cout << "\rprocessing TF " << ++iTF << " / " << nTF << "..." << std::flush;

    for (uint32_t iMCHTrack = 0; iMCHTrack < allTracks->size(); ++iMCHTrack) {

      // select muon tracks if requested
      auto muon = FindMuon(iMCHTrack, *allMuons);
      if (selectMatched && !muon) {
        continue;
      }

      const auto& mchTrack = (*allTracks)[iMCHTrack];
      gsl::span<Cluster> trackClusters(&(*allTrackClusters)[mchTrack.getFirstClusterIdx()], mchTrack.getNClusters());

      // create an internal track and compute its parameters at each attached cluster
      TrackStruct track{};
      if (!ClustersToTrack(trackClusters, track)) {
        continue;
      }

      // propagate the track to the vertex (0,0,0)
      if (!ExtrapToVertex(track)) {
        continue;
      }

      // apply selections if requested
      if (applyTrackSelection && !IsSelected(track)) {
        continue;
      }

      // fill output track tree
      trackOut.param = track.paramAtVtx.getTrackParamStruct();
      trackOut.rAbs = track.rAbs;
      trackOut.dca = track.dca;
      trackOut.pUncorr = mchTrack.getP();
      trackOut.nClusters = mchTrack.getNClusters();
      trackOut.chi2 = mchTrack.getChi2OverNDF();
      trackOut.matchChi2 = muon ? muon->getMatchChi2OverNDF() : 0.;
      dataTrTree.Fill();

      // fill output cluster tree
      for (const auto& param : track.track) {
        trackParamOut = param.getTrackParamStruct();
        clusterOut = *(param.getClusterPtr());
        dataClTree.Fill();
      }
    }
  }
  cout << "\r\033[Kprocessing completed" << endl;

  dataTrFile.Write();
  dataTrFile.Close();
  dataClFile.Write();
  dataClFile.Close();
  fTracks->Close();
  fMuons->Close();
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
bool ClustersToTrack(gsl::span<Cluster> trackClusters, TrackStruct& track)
{
  /// fill the internal track with the associated clusters and fit it
  /// return false if the fit fails

  static TrackFitter trackFitter{};
  trackFitter.smoothTracks(true);
  TrackExtrap::useExtrapV2();

  for (auto& cluster : trackClusters) {
    // if (cluster.ey < 1.) {
    //   cluster.ey = 0.05;
    // }
    // if (cluster.getChamberId() == 2 || cluster.getChamberId() == 3) {
    //   cluster.y += gRandom->Gaus(0., 0.02);
    //   cluster.ey *= sqrt(2.);
    // }
    track.track.createParamAtCluster(cluster);
  }

  try {
    trackFitter.fit(track.track);
    return true;
  } catch (exception const& e) {
    return false;
  }
}

//_________________________________________________________________________________________________
bool ExtrapToVertex(TrackStruct& track)
{
  /// compute the track parameters at vertex, at DCA and at the end of the absorber

  TrackExtrap::useExtrapV2(false);

  // extrapolate to vertex
  track.paramAtVtx = track.track.first();
  if (!TrackExtrap::extrapToVertex(track.paramAtVtx, 0., 0., 0., 0., 0.)) {
    return false;
  }

  // extrapolate to DCA
  TrackParam paramAtDCA(track.track.first());
  if (!TrackExtrap::extrapToVertexWithoutBranson(paramAtDCA, 0.)) {
    return false;
  }
  double dcaX = paramAtDCA.getNonBendingCoor();
  double dcaY = paramAtDCA.getBendingCoor();
  track.dca = TMath::Sqrt(dcaX * dcaX + dcaY * dcaY);

  // extrapolate to the end of the absorber
  TrackParam paramAtAbs(track.track.first());
  if (!TrackExtrap::extrapToZ(paramAtAbs, -505.)) {
    return false;
  }
  double xAbs = paramAtAbs.getNonBendingCoor();
  double yAbs = paramAtAbs.getBendingCoor();
  track.rAbs = TMath::Sqrt(xAbs * xAbs + yAbs * yAbs);

  return true;
}

//_________________________________________________________________________________________________
bool IsSelected(const TrackStruct& track)
{
  /// apply standard track selections + pDCA + p > 10 GeV/c

  static const double sigmaPDCA23 = 80.;
  static const double sigmaPDCA310 = 54.;
  static const double nSigmaPDCA = 6.;
  static const double relPRes = 0.0004;
  static const double slopeRes = 0.0005;

  double thetaAbs = TMath::ATan(track.rAbs / 505.) * TMath::RadToDeg();
  if (thetaAbs < 2. || thetaAbs > 10.) {
    return false;
  }

  double p = track.paramAtVtx.p();
  if (p < 10.) {
    return false;
  }

  double px = track.paramAtVtx.px();
  double py = track.paramAtVtx.py();
  double pT = TMath::Sqrt(px * px + py * py);
  if (pT < 1.) {
    return false;
  }

  double pz = track.paramAtVtx.pz();
  double eta = 0.5 * log((p + pz) / (p - pz));
  if (eta < -4. || eta > -2.5) {
    return false;
  }

  double pDCA = track.track.first().p() * track.dca;
  double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
  double nrp = nSigmaPDCA * relPRes * p;
  double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
  double slopeResEffect = 535. * slopeRes * p;
  double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
  if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
int IRtoBCinTF(o2::InteractionRecord ir, int run, bool isMC)
{
  /// convert IR to BC from TF start

  static const std::unordered_map<int, uint32_t> firstTForbit0{
    {295584, 0},
    {529691, 74406400},
    {544490, 26693600},
    {549586, 67400512},
    {559410, 527962048},
    {562849, 179331744},
    {562967, 47243328},
    {562968, 62076224},
    {562969, 74244704},
    {562970, 86353952},
    {562971, 97188864},
    {562972, 109218272},
    {562975, 146579680},
    {562976, 160383072},
    {562979, 193958464}};

  static const std::unordered_map<int, int> nOrbitsPerTF{
    {295584, 1},
    {529691, 128},
    {544490, 32},
    {549586, 32},
    {559410, 32},
    {562849, 32},
    {562967, 32},
    {562968, 32},
    {562969, 32},
    {562970, 32},
    {562971, 32},
    {562972, 32},
    {562975, 32},
    {562976, 32},
    {562979, 32}};

  auto itFirstTForbit0 = firstTForbit0.find(run);
  auto itNOrbitsPerTF = nOrbitsPerTF.find(run);
  if (itFirstTForbit0 == firstTForbit0.end() || itNOrbitsPerTF == nOrbitsPerTF.end()) {
    LOGP(error, "time information not provided for that run");
    exit(1);
  }

  if (ir.orbit < itFirstTForbit0->second) {
    LOGP(error, "IR orbit {} < first orbit of first TF {} !?", ir.orbit, itFirstTForbit0->second);
    exit(2);
  }

  int orbitInTF = ir.orbit - itFirstTForbit0->second;
  if (!isMC) {
    orbitInTF = orbitInTF % itNOrbitsPerTF->second;
  }

  return orbitInTF * o2::constants::lhc::LHCMaxBunches + ir.bc;
}
