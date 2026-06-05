#ifndef DIGITUTILS_H_
#define DIGITUTILS_H_

#include <algorithm>
#include <cmath>
#include <iterator>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <gsl/span>

#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>

#include "DataFormatsMCH/Digit.h"
#include "MCHMappingInterface/Segmentation.h"

using o2::mch::Digit;

/*
 * This file contains utility functions to produce digit control plots.
 * To be used, they require the MCH mapping to be loaded:
 * gSystem->Load("libO2MCHMappingImpl4")
 */

static const std::vector<double> asymmLimits{-0.3, -0.2, -0.1, -0.025, 0.025, 0.1, 0.2, 0.3};
static const std::vector<int> asymmColors{2, 3, 4, 5, 1, 6, 7, 8, 9};

//_________________________________________________________________________________________________
double adcToCharge(uint32_t adc)
{
  /// function to reinterpret digit ADC as calibrated charge in run2
  /// then convert it back in ADC unit using run2 electronic gain

  // 1 ADC channel = 0.61 mV; capa = 0.2 and a0 = 1.25, which is equivalent to gain = 4 mV/fC
  static const double fc2adc = 1. / (1.25 * 0.2 * 0.61);

  float charge(0.);
  std::memcpy(&charge, &adc, sizeof(adc));

  return static_cast<double>(charge) * fc2adc;
}

//_________________________________________________________________________________________________
void CreateDigitTimeInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create digit time histograms

  histos.emplace_back(new TH1F(Form("time%s", extension), "#Deltat vs track time;#Deltat (BC)", 101, -50.5, 50.5));

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillDigitTimeInfo(const Digit& digit, int trackTime, gsl::span<TH1*> histos)
{
  /// fill digit time histograms

  histos[0]->Fill(digit.getTime() + 1.5 - trackTime);
}

//_________________________________________________________________________________________________
TCanvas* DrawDigitTimeInfo(gsl::span<TH1*> histos, const char* extension = "")
{
  /// draw digit time histograms

  TCanvas* c = new TCanvas(Form("digitTimeInfo%s", extension),
                           Form("digit time characteristics %s", extension), 10, 10, 400, 400);
  gPad->SetLogy();
  histos[0]->Draw();

  return c;
}

//_________________________________________________________________________________________________
void CreateDigitChargeInfo(std::vector<TH1*>& histos, const char* extension = "")
{
  /// create digit charge histograms

  histos.emplace_back(new TH1F(Form("ADC%s", extension), "ADC;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("ADCB%s", extension), "ADC Bending;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("ADCNB%s", extension), "ADC Non Bending;ADC", 10001, -0.5, 10000.5));
  histos.emplace_back(new TH1F(Form("Samples%s", extension), "N samples;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH1F(Form("SamplesB%s", extension), "N samples Bending;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH1F(Form("SamplesNB%s", extension), "N samples Non Bending;N samples", 1024, -0.5, 1023.5));
  histos.emplace_back(new TH2F(Form("ADCvsSampleB%s", extension), "ADC vs N samples Bending;N samples;ADC", 401, -0.5, 400.5, 4001, -0.5, 40009.5));
  histos.emplace_back(new TH2F(Form("ADCvsSampleNB%s", extension), "ADC vs N samples Non Bending;N samples;ADC", 401, -0.5, 400.5, 4001, -0.5, 40009.5));

  for (int i = 0; i <= (int)asymmLimits.size(); ++i) {
    histos.emplace_back(new TH2F(Form("ADCvsSampleB%d%s", i, extension), "ADC vs N samples Bending;N samples;ADC", 401, -0.5, 400.5, 2001, -0.5, 40019.5));
    histos.emplace_back(new TH2F(Form("ADCvsSampleNB%d%s", i, extension), "ADC vs N samples Non Bending;N samples;ADC", 401, -0.5, 400.5, 2001, -0.5, 40019.5));
  }

  for (auto h : histos) {
    h->SetDirectory(0);
  }
}

//_________________________________________________________________________________________________
void FillDigitChargeInfo(const Digit& digit, gsl::span<TH1*> histos, double chargeAsymm, bool run2 = false)
{
  /// fill digit charge histograms

  const auto& segmentation = o2::mch::mapping::segmentation(digit.getDetID());

  auto charge = run2 ? adcToCharge(digit.getADC()) : digit.getADC();
  auto i = std::distance(asymmLimits.begin(), std::upper_bound(asymmLimits.begin(), asymmLimits.end(), chargeAsymm));

  histos[0]->Fill(charge);
  histos[3]->Fill(digit.getNofSamples());
  if (segmentation.isBendingPad(digit.getPadID())) {
    histos[1]->Fill(charge);
    histos[4]->Fill(digit.getNofSamples());
    histos[6]->Fill(digit.getNofSamples(), charge);
    histos[2 * i + 8]->Fill(digit.getNofSamples(), charge);
  } else {
    histos[2]->Fill(charge);
    histos[5]->Fill(digit.getNofSamples());
    histos[7]->Fill(digit.getNofSamples(), charge);
    histos[2 * i + 9]->Fill(digit.getNofSamples(), charge);
  }
}

//_________________________________________________________________________________________________
TCanvas* DrawDigitChargeInfo(gsl::span<TH1*> histos, const char* extension = "")
{
  /// draw digit charge histograms

  static int logy[3] = {1, 1, 0};
  static int logz[3] = {0, 0, 1};
  static const char* opt[3] = {"", "", "colz"};

  static std::vector<std::string> asymmLegend{};
  if (asymmLegend.empty()) {
    asymmLegend.emplace_back(fmt::format("< {}", asymmLimits[0]));
    for (int i = 1; i < (int)asymmLimits.size(); ++i) {
      asymmLegend.emplace_back(fmt::format("[{}, {}[", asymmLimits[i - 1], asymmLimits[i]));
    }
    asymmLegend.emplace_back(fmt::format(">= {}", asymmLimits[asymmLimits.size() - 1]));
  }

  TCanvas* c = new TCanvas(Form("digitChargeInfo%s", extension),
                           Form("digit charge characteristics %s", extension), 10, 10, 1200, 800);
  c->Divide(3, 2);

  TLegend* l = new TLegend(0.5, 0.7, 0.9, 0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);

  c->cd(1);
  gPad->SetLogy();
  histos[0]->SetLineColor(1);
  histos[0]->Draw();
  l->AddEntry(histos[0], "all", "l");
  histos[1]->SetStats(false);
  histos[1]->SetLineColor(2);
  histos[1]->Draw("sames");
  l->AddEntry(histos[1], "bending", "l");
  histos[2]->SetStats(false);
  histos[2]->SetLineColor(4);
  histos[2]->Draw("sames");
  l->AddEntry(histos[2], "non bending", "l");
  l->Draw("same");

  c->cd(4);
  gPad->SetLogy();
  histos[3]->SetLineColor(1);
  histos[3]->Draw();
  histos[4]->SetStats(false);
  histos[4]->SetLineColor(2);
  histos[4]->Draw("sames");
  histos[5]->SetStats(false);
  histos[5]->SetLineColor(4);
  histos[5]->Draw("sames");

  c->cd(2);
  gPad->SetLogz();
  gPad->SetGrid();
  histos[6]->Draw("colz");

  c->cd(5);
  gPad->SetLogz();
  gPad->SetGrid();
  histos[7]->Draw("colz");

  TLegend* l2 = new TLegend(0.65, 0.1, 0.9, 0.4);
  l2->SetFillStyle(0);
  l2->SetBorderSize(0);

  for (int i = 0; i <= (int)asymmLimits.size(); ++i) {
    c->cd(3);
    gPad->SetGrid();
    histos[2 * i + 8]->SetStats(false);
    histos[2 * i + 8]->SetLineColor(asymmColors[i]);
    histos[2 * i + 8]->SetLineWidth(2);
    histos[2 * i + 8]->Draw((i == 0) ? "box" : "boxsames0");
    l2->AddEntry(histos[2 * i + 8], asymmLegend[i].c_str(), "l");

    c->cd(6);
    gPad->SetGrid();
    histos[2 * i + 9]->SetStats(false);
    histos[2 * i + 9]->SetLineColor(asymmColors[i]);
    histos[2 * i + 9]->SetLineWidth(2);
    histos[2 * i + 9]->Draw((i == 0) ? "box" : "boxsames0");
  }

  c->cd(3);
  l2->Draw("same");

  return c;
}

#endif // DIGITUTILS_H_
