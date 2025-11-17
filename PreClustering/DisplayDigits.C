#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>

#include <gsl/span>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "Framework/Logger.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsMID/ROFRecord.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/Digit.h"

using namespace o2;

// const uint32_t nOrbitsPerTF = 128;
const uint32_t nOrbitsPerTF = 32;

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
  {562967, 47243328},
  {562968, 62076224},
  {562970, 86353952},
  {562971, 97188864},
  {562975, 146579680},
  {562976, 160383072}};

uint32_t firstTForbit0(0);

std::tuple<TFile*, TTreeReader*> LoadData(const char* fileName, const char* treeName);
void CreateTimeHistos(std::vector<TH1*>& histos, int BCRange[2], int iTF, const char* extension);
void FillTimeHistos(const std::vector<mch::Digit>& digits, bool selectSignal, bool rejectBackground,
                    gsl::span<TH1*> histos, int deMin = 100, int deMax = 1025);
void FillTimeHistos(const std::vector<mid::ROFRecord>& midROFs, TH1* histo, bool isMC);
void DrawTimeHistos(TH1* hBC, std::vector<TH1*>& histos, std::vector<TH1*>* histosSt);
void CreateCorrelationHistos(std::vector<TH1*>& histos, const char* extension);
void FillCorrelationHistos(const std::vector<mch::Digit>& digits, bool selectSignal, bool rejectBackground,
                           TH1* hist, int deMin = 100, int deMax = 1025);
void FillCorrelationHistos(std::vector<TH1*>& timeHistos, gsl::span<TH1*> corrHistos);
void DrawCorrelationHistos(gsl::span<TH1*> histos, const char* extension);
void CreateChargeHistos(std::vector<TH1*>& histos, const char* extension);
void FillChargeHistos(const std::vector<mch::Digit>& digits, bool selectSignal, bool rejectBackground,
                      gsl::span<TH1*> histos, int deMin = 100, int deMax = 1025);
void DrawChargeHistos(std::vector<TH1*>& histos, const char* extension);
double signalCut(double* x, double* p);
double backgroundCut(double* x, double* p);

int nHistosPerTF = -1; // number of time histos per TF

// parameters defining signal selection and background rejection
uint16_t minNSamplesSignal = 17;
double signalParam[4] = {80., 16., 12., 1.2};
uint16_t minNSamplesBackground = 14;
double backgroundParam[4] = {18., 24., -20., 7.};

int bcIntegrationRange = 6; // time window ([-range, range]) to integrate digits

//_________________________________________________________________________________________________
void DisplayDigits(int runNumber, std::string mchDigitFileName = "mchfdigits.root", std::string mchClusterFileName = "",
                   std::string mchTrackFileName = "mchtracks.root", std::string midTrackFileName = "mid-reco.root",
                   bool selectSignal = false, bool rejectBackground = false, bool isMC = false)
{
  /// show the characteristics of the reconstructed digits

  // find the first orbit of the first TF of this run
  auto itOrbit0 = firstTForbit0perRun.find(runNumber);
  if (itOrbit0 != firstTForbit0perRun.end()) {
    firstTForbit0 = itOrbit0->second;
  } else {
    LOG(warning) << "first orbit not found for this run";
  }

  // load data
  auto [fMCHD, mchReaderD] = LoadData(mchDigitFileName.c_str(), "o2sim");
//  TTreeReaderValue<std::vector<mch::ROFRecord>> mchDigitROFs = {*mchReaderD, "rofs"};
  TTreeReaderValue<std::vector<mch::Digit>> mchDigits = {*mchReaderD, "MCHDigit"};
  TFile* fMCHT = nullptr;
  TTreeReader* mchReaderT = nullptr;
  std::unique_ptr<TTreeReaderValue<std::vector<mch::Digit>>> mchTrackDigits{};
  if (mchClusterFileName.empty()) {
    // read digits from the mchtracks.root
    std::tie(fMCHT, mchReaderT) = LoadData(mchTrackFileName.c_str(), "o2sim");
    if (!mchReaderT->GetTree()->FindBranch("trackdigits")) {
      LOG(error) << "track digits are missing";
      exit(-1);
    }
    mchTrackDigits = std::make_unique<TTreeReaderValue<std::vector<mch::Digit>>>(*mchReaderT, "trackdigits");
  } else {
    // read digits from the mchclusters.root
    LOG(info) << "draw cluster digits instead of track digits";
    std::tie(fMCHT, mchReaderT) = LoadData(mchClusterFileName.c_str(), "o2sim");
    mchTrackDigits = std::make_unique<TTreeReaderValue<std::vector<mch::Digit>>>(*mchReaderT, "clusterdigits");
  }
  auto [fMID, midReader] = LoadData(midTrackFileName.c_str(), "midreco");
  TTreeReaderValue<std::vector<mid::ROFRecord>> midROFs = {*midReader, "MIDTrackROF"};
//  TTreeReaderValue<std::vector<mid::ROFRecord>> midROFs = {*midReader, "MIDTrackClusterROF"};
  int nTF = midReader->GetEntries(false);
  if (mchReaderD->GetEntries(false) != nTF || mchReaderT->GetEntries(false) != nTF) {
    LOG(error) << " input files do not contain the same number of TF";
    exit(-1);
  }
  /*
    // first loop to set the time range
    int64_t BCRange[2] = {0x7FFFFFFFFFFFFFFF, -0x7FFFFFFFFFFFFFFF};
    while (mchReaderD->Next()) {

      // find the reference TF used to defined the digit time
      if (mchDigitROFs->size() == 0) {
        continue;
      }
      const auto& rof = mchDigitROFs->back();
      if (rof.getBCData().orbit < firstTForbit0) {
        LOG(error) << "MCH IR orbit < first orbit of first TF !?";
        exit(-1);
      }
      const auto& digit = mchDigits->at(rof.getLastIdx());
      if (digit.getTime() < 0) {
        LOG(info) << "unable to find the reference TF because of a bug affecting digits with negative time";
        exit(-1);
      }
      uint32_t iTF0 = (rof.getBCData().orbit - firstTForbit0) / nOrbitsPerTF;
      static const uint32_t tfDuration = nOrbitsPerTF * constants::lhc::LHCMaxBunches;
      uint32_t iTFCorr = digit.getTime() / tfDuration;
      InteractionRecord tfRef(0, firstTForbit0 + (iTF0 - iTFCorr) * nOrbitsPerTF);
      int64_t bcRef = tfRef.toLong();

      // set the time range
      for (const auto& digit : *mchDigits) {
        int64_t bc = bcRef + digit.getTime();
        if (bc < BCRange[0]) {
          BCRange[0] = bc;
        }
        if (bc > BCRange[1]) {
          BCRange[1] = bc;
        }
      }
    }

    if (BCRange[0] > BCRange[1]) {
      LOG(error) << "nothing to display !?";
      exit(-1);
    }

    mchReaderD->Restart();
  */
  // create histograms
  std::vector<TH1*> timeHistos{};
  std::vector<TH1*> timeHistosSt[3] = {{}, {}, {}};
//  int BCRange[2] = {-1100000, 1600000};
  int BCRange[2] = {0, nOrbitsPerTF * constants::lhc::LHCMaxBunches};
  TH1* hBC = new TH1F("hBC", "hBC", BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5);
  hBC->SetDirectory(0);
  for (uint32_t i = 0; i < nOrbitsPerTF; ++i) {
    hBC->Fill(i * constants::lhc::LHCMaxBunches + 251, 1000000);
  }
  std::vector<TH1*> chargeHistos{};
  CreateChargeHistos(chargeHistos, "all");
  CreateChargeHistos(chargeHistos, "St1");
  CreateChargeHistos(chargeHistos, "St2");
  CreateChargeHistos(chargeHistos, "St345");
  std::vector<TH1*> chargeHistosCh1{};
  CreateChargeHistos(chargeHistosCh1, "DE100");
  CreateChargeHistos(chargeHistosCh1, "DE101");
  CreateChargeHistos(chargeHistosCh1, "DE102");
  CreateChargeHistos(chargeHistosCh1, "DE103");
  std::vector<TH1*> chargeHistosCh2{};
  CreateChargeHistos(chargeHistosCh2, "DE200");
  CreateChargeHistos(chargeHistosCh2, "DE201");
  CreateChargeHistos(chargeHistosCh2, "DE202");
  CreateChargeHistos(chargeHistosCh2, "DE203");
  std::vector<TH1*> corrHistos{};
  CreateCorrelationHistos(corrHistos, "all");
  CreateCorrelationHistos(corrHistos, "St1");
  CreateCorrelationHistos(corrHistos, "St2");
  CreateCorrelationHistos(corrHistos, "St345");

  int iTF(-1);
  while (mchReaderD->Next() && mchReaderT->Next() && midReader->Next()) {
    ++iTF;

    if (!(*mchTrackDigits)->empty() && iTF == 0) {
      CreateTimeHistos(timeHistos, BCRange, iTF, "all");
      FillTimeHistos(*mchDigits, selectSignal, rejectBackground, {&timeHistos[timeHistos.size() - 9], 4});
      FillTimeHistos(**mchTrackDigits, selectSignal, rejectBackground, {&timeHistos[timeHistos.size() - 5], 4});
      FillTimeHistos(*midROFs, timeHistos[timeHistos.size() - 1], isMC);

      CreateTimeHistos(timeHistosSt[0], BCRange, iTF, "St1");
      FillTimeHistos(*mchDigits, selectSignal, rejectBackground, {&timeHistosSt[0][timeHistosSt[0].size() - 9], 4}, 100, 203);
      FillTimeHistos(**mchTrackDigits, selectSignal, rejectBackground, {&timeHistosSt[0][timeHistosSt[0].size() - 5], 4}, 100, 203);
      FillTimeHistos(*midROFs, timeHistosSt[0][timeHistosSt[0].size() - 1], isMC);

      CreateTimeHistos(timeHistosSt[1], BCRange, iTF, "St2");
      FillTimeHistos(*mchDigits, selectSignal, rejectBackground, {&timeHistosSt[1][timeHistosSt[1].size() - 9], 4}, 300, 403);
      FillTimeHistos(**mchTrackDigits, selectSignal, rejectBackground, {&timeHistosSt[1][timeHistosSt[1].size() - 5], 4}, 300, 403);
      FillTimeHistos(*midROFs, timeHistosSt[1][timeHistosSt[1].size() - 1], isMC);

      CreateTimeHistos(timeHistosSt[2], BCRange, iTF, "St345");
      FillTimeHistos(*mchDigits, selectSignal, rejectBackground, {&timeHistosSt[2][timeHistosSt[2].size() - 9], 4}, 500, 1025);
      FillTimeHistos(**mchTrackDigits, selectSignal, rejectBackground, {&timeHistosSt[2][timeHistosSt[2].size() - 5], 4}, 500, 1025);
      FillTimeHistos(*midROFs, timeHistosSt[2][timeHistosSt[2].size() - 1], isMC);

      FillCorrelationHistos(**mchTrackDigits, selectSignal, rejectBackground, corrHistos[3]);
      FillCorrelationHistos(**mchTrackDigits, selectSignal, rejectBackground, corrHistos[7], 100, 203);
      FillCorrelationHistos(**mchTrackDigits, selectSignal, rejectBackground, corrHistos[11], 300, 403);
      FillCorrelationHistos(**mchTrackDigits, selectSignal, rejectBackground, corrHistos[15], 500, 1025);
    }

    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistos[0], 3});
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistos[3], 3}, 100, 203);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistos[6], 3}, 300, 403);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistos[9], 3}, 500, 1025);

    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh1[0], 3}, 100, 100);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh1[3], 3}, 101, 101);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh1[6], 3}, 102, 102);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh1[9], 3}, 103, 103);

    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh2[0], 3}, 200, 200);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh2[3], 3}, 201, 201);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh2[6], 3}, 202, 202);
    FillChargeHistos(*mchDigits, selectSignal, rejectBackground, {&chargeHistosCh2[9], 3}, 203, 203);
  }

  fMCHD->Close();
  fMCHT->Close();

  FillCorrelationHistos(timeHistos, {&corrHistos[0], 3});
  FillCorrelationHistos(timeHistosSt[0], {&corrHistos[4], 3});
  FillCorrelationHistos(timeHistosSt[1], {&corrHistos[8], 3});
  FillCorrelationHistos(timeHistosSt[2], {&corrHistos[12], 3});

  DrawChargeHistos(chargeHistos, "PerStation");
  DrawChargeHistos(chargeHistosCh1, "Ch1");
  DrawChargeHistos(chargeHistosCh2, "Ch2");
  DrawTimeHistos(hBC, timeHistos, timeHistosSt);
  DrawCorrelationHistos({&corrHistos[0], 4}, "all");
  DrawCorrelationHistos({&corrHistos[4], 4}, "St1");
  DrawCorrelationHistos({&corrHistos[8], 4}, "St2");
  DrawCorrelationHistos({&corrHistos[12], 4}, "St345");
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
void CreateTimeHistos(std::vector<TH1*>& histos, int BCRange[2], int iTF, const char* extension)
{
  /// create digit time histograms

  nHistosPerTF = 9;

  histos.emplace_back(new TH1F(Form("nDigitsTF%d%s", iTF, extension), Form("nDigits - TF %d - %s;BC;nDigits", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("ChargeTF%d%s", iTF, extension), Form("Charge - TF %d - %s;BC;ADC", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("nDigitsIntTF%d%s", iTF, extension), Form("nDigits integral - TF %d - %s;BC;nDigits", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("ChargeIntTF%d%s", iTF, extension), Form("Charge integral - TF %d - %s;BC;ADC", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("nTrackDigitsTF%d%s", iTF, extension), Form("nTrackDigits - TF %d - %s;BC;nTrackDigits", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("TrackChargeTF%d%s", iTF, extension), Form("TrackCharge - TF %d - %s;BC;ADC", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("nTrackDigitsIntTF%d%s", iTF, extension), Form("nTrackDigits integral - TF %d - %s;BC;nTrackDigits", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("TrackChargeIntTF%d%s", iTF, extension), Form("TrackCharge integral - TF %d - %s;BC;ADC", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));
  histos.emplace_back(new TH1F(Form("nMIDObjTF%d%s", iTF, extension), Form("nMIDObj - TF %d - %s;BC;nMIDObj", iTF, extension),
                               BCRange[1] - BCRange[0], BCRange[0] - 0.5, BCRange[1] - 0.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillTimeHistos(const std::vector<mch::Digit>& digits, bool selectSignal, bool rejectBackground,
                    gsl::span<TH1*> histos, int deMin, int deMax)
{
  /// fill digit time histogram

  for (const auto& digit : digits) {
    if (digit.getDetID() < deMin || digit.getDetID() > deMax) {
      continue;
    }
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
    histos[0]->Fill(digit.getTime());
    histos[1]->Fill(digit.getTime(), digit.getADC());
    for (int bc = -bcIntegrationRange; bc <= bcIntegrationRange; ++bc) {
      histos[2]->Fill(digit.getTime() + bc);
      histos[3]->Fill(digit.getTime() + bc, digit.getADC());
    }
  }
}

//_________________________________________________________________________________________________
void FillTimeHistos(const std::vector<mid::ROFRecord>& midROFs, TH1* histo, bool isMC)
{
  /// fill mid time histogram

  for (const auto& rof : midROFs) {
    if (rof.nEntries < 1) {
      continue;
    }
    uint32_t orbitInTF = isMC ? rof.interactionRecord.orbit : (rof.interactionRecord.orbit - firstTForbit0) % nOrbitsPerTF;
    histo->Fill(orbitInTF * constants::lhc::LHCMaxBunches + rof.interactionRecord.bc, /*1000000*/rof.nEntries);
  }
}

//_________________________________________________________________________________________________
void DrawTimeHistos(TH1* hBC, std::vector<TH1*>& histos, std::vector<TH1*>* histosSt)
{
  /// draw digit time histograms

  int nTF = histos.size() / nHistosPerTF;
  for (int iTF = 0; iTF < nTF; ++iTF) {
    static const int nHist = 4;
    TCanvas* cHist[nHist] = {new TCanvas(Form("cTime%d", iTF), Form("cTime%d", iTF), 10, 10, 1200, 1200),
                             new TCanvas(Form("cCharge%d", iTF), Form("cCharge%d", iTF), 10, 10, 1200, 1200),
                             new TCanvas(Form("cTimeIntegral%d", iTF), Form("cTimeIntegral%d", iTF), 10, 10, 1200, 1200),
                             new TCanvas(Form("cChargeIntegral%d", iTF), Form("cChargeIntegral%d", iTF), 10, 10, 1200, 1200)};
    for (int i = 0; i < nHist; ++i) {
      cHist[i]->Divide(1, 4);
      cHist[i]->cd(1);
      histos[nHistosPerTF * iTF + i]->SetStats(false);
      histos[nHistosPerTF * iTF + i]->SetLineColor(4);
      histos[nHistosPerTF * iTF + i]->Draw("hist");
      hBC->SetLineColor(5);
      hBC->SetLineStyle(2);
      hBC->Draw("histsame");
      histos[nHistosPerTF * (iTF + 1) - 1]->SetLineColor(3);
      histos[nHistosPerTF * (iTF + 1) - 1]->SetLineStyle(2);
      histos[nHistosPerTF * (iTF + 1) - 1]->Draw("histsame");
      histos[nHistosPerTF * iTF + i + nHist]->SetLineColor(2);
      histos[nHistosPerTF * iTF + i + nHist]->Draw("histsame");
      for (int iSt = 0; iSt < 3; ++iSt) {
        cHist[i]->cd(iSt + 2);
        histosSt[iSt][nHistosPerTF * iTF + i]->SetStats(false);
        histosSt[iSt][nHistosPerTF * iTF + i]->SetLineColor(4);
        histosSt[iSt][nHistosPerTF * iTF + i]->Draw("hist");
        hBC->SetLineColor(5);
        hBC->SetLineStyle(2);
        hBC->Draw("histsame");
        histosSt[iSt][nHistosPerTF * (iTF + 1) - 1]->SetLineColor(3);
        histosSt[iSt][nHistosPerTF * (iTF + 1) - 1]->SetLineStyle(2);
        histosSt[iSt][nHistosPerTF * (iTF + 1) - 1]->Draw("histsame");
        histosSt[iSt][nHistosPerTF * iTF + i + nHist]->SetLineColor(2);
        histosSt[iSt][nHistosPerTF * iTF + i + nHist]->Draw("histsame");
      }
    }
  }
}

//_________________________________________________________________________________________________
void CreateCorrelationHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create correlation histograms between number of digits and total charge in BC time window

  histos.emplace_back(new TH2F(Form("ChargevsNDigits%s", extension), Form("Charge vs N digits per BC (%s);N digits;ADC", extension),
                               1000, 0, 1000, 10000, 0, 100000));
  histos.emplace_back(new TH2F(Form("ChargevsNDigitsInt%s", extension), Form("Charge vs N digits in BC window (%s);N digits;ADC", extension),
                               1000, 0, 1000, 10000, 0, 100000));
  histos.emplace_back(new TH2F(Form("TrackChargevsNDigits%s", extension), Form("Track charge vs N digits per BC (%s);N digits;ADC", extension),
                               1000, 0, 1000, 10000, 0, 100000));
  histos.emplace_back(new TH2F(Form("TrackChargevsNDigitsInt%s", extension), Form("Integrated track charge vs N digits (%s);N digits;ADC", extension),
                               1000, 0, 1000, 10000, 0, 100000));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillCorrelationHistos(const std::vector<mch::Digit>& digits, bool selectSignal, bool rejectBackground,
                           TH1* hist, int deMin, int deMax)
{
  /// fill correlation histograms between number of digits and total charge associated with tracks

  uint32_t charge(0);
  int nDigits(0);

  for (const auto& digit : digits) {
    if (digit.getDetID() < deMin || digit.getDetID() > deMax) {
      continue;
    }
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
    charge += digit.getADC();
    ++nDigits;
  }

  hist->Fill(nDigits, charge);
}

//_________________________________________________________________________________________________
void FillCorrelationHistos(std::vector<TH1*>& timeHistos, gsl::span<TH1*> corrHistos)
{
  /// fill correlation histograms between number of digits and total charge in BC time window

  if (timeHistos.empty()) {
    return;
  }

  int nBC = timeHistos[0]->GetNbinsX();
  int nTF = timeHistos.size() / nHistosPerTF;
  for (int iTF = 0; iTF < nTF; ++iTF) {
    for (int iBC = 1; iBC <= nBC; ++iBC) {
      for (int i = 0; i < 3; ++i) {
        if (timeHistos[nHistosPerTF * iTF + 2 * i]->GetBinContent(iBC) > 0) {
          corrHistos[i]->Fill(timeHistos[nHistosPerTF * iTF + 2 * i]->GetBinContent(iBC),
                              timeHistos[nHistosPerTF * iTF + 2 * i + 1]->GetBinContent(iBC));
        }
      }
    }
  }
}

//_________________________________________________________________________________________________
void DrawCorrelationHistos(gsl::span<TH1*> histos, const char* extension)
{
  /// draw correlation histograms between number of digits and total charge in BC time window

  TCanvas* cCorr = new TCanvas(Form("cCorr%s", extension), Form("cCorr%s", extension), 10, 10, 800, 800);
  cCorr->Divide(2, 2);
  for (int i = 0; i < 4; ++i) {
    cCorr->cd(i + 1);
    gPad->SetLogz();
    histos[i]->Draw("boxcolz");
  }
}

//_________________________________________________________________________________________________
void CreateChargeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// create digit charge histograms

  histos.emplace_back(new TH1D(Form("ADC%s", extension), Form("ADC (%s);ADC", extension), 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1D(Form("Samples%s", extension), Form("N samples (%s);N samples", extension), 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSample%s", extension), Form("ADC vs N samples (%s);N samples;ADC", extension), 1024, -0.5, 1023.5, 10001, -0.5, 100009.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillChargeHistos(const std::vector<mch::Digit>& digits, bool selectSignal, bool rejectBackground,
                      gsl::span<TH1*> histos, int deMin, int deMax)
{
  /// fill digit charge histograms
  for (const auto digit : digits) {
    if (digit.getDetID() >= deMin && digit.getDetID() <= deMax) {
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
      histos[0]->Fill(digit.getADC());
      histos[1]->Fill(digit.getNofSamples());
      histos[2]->Fill(digit.getNofSamples(), digit.getADC());
    }
  }
}

//_________________________________________________________________________________________________
void DrawChargeHistos(std::vector<TH1*>& histos, const char* extension)
{
  /// draw digit charge histograms

  int nBinsX = histos.size() / 3;
  TCanvas* cHist = new TCanvas(Form("cCharge%s", extension), Form("cCharge%s", extension), 10, 10, nBinsX * 300, 900);
  cHist->Divide(nBinsX, 3);
  for (int iBinX = 0; iBinX < nBinsX; ++iBinX) {
    cHist->cd(iBinX + 1);
    gPad->SetLogy();
    histos[3 * iBinX]->SetStats(false);
    histos[3 * iBinX]->Draw();
    cHist->cd(nBinsX + iBinX + 1);
    gPad->SetLogy();
    histos[3 * iBinX + 1]->SetStats(false);
    histos[3 * iBinX + 1]->Draw();
    cHist->cd(2 * nBinsX + iBinX + 1);
    gPad->SetLogz();
    histos[3 * iBinX + 2]->SetStats(false);
    histos[3 * iBinX + 2]->Draw("colz");
  }

  static TF1* fSignal = new TF1("fSignal", signalCut, 0, 1023, 4);
  fSignal->SetParameters(signalParam);
  fSignal->SetLineColor(2);
  static TF1* fBackground = new TF1("fBackground", backgroundCut, 0, 1023, 4);
  fBackground->SetParameters(backgroundParam);
  fBackground->SetLineColor(4);
  for (int iBinX = 0; iBinX < nBinsX; ++iBinX) {
    cHist->cd(2 * nBinsX + iBinX + 1);
    fSignal->Draw("same");
    fBackground->Draw("same");
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
