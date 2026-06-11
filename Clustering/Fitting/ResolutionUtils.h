#ifndef RESOLUTIONUTILS_H_
#define RESOLUTIONUTILS_H_

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <TF1.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TList.h>
#include <TString.h>

#include "DataFormatsMCH/Digit.h"

using o2::mch::Digit;

//_________________________________________________________________________________________________
// create the THnSparse (10 axes) to extract the resolution in the residuals later
THnSparseD* CreatePreClusterInfoMULTI(const char* extension = "")
{
  const Int_t nDim = 10;

  Int_t nbins[nDim] = {
    1000,  // pvalue
    1122,  // residuals
    10001, // ADC_fit
    10001, // ADC_mes
    5001,  // ADC_cluster
    301,   // nSamples
    280,   // Asymm
    3,     // Wire
    2,     // Cathode
    500    // fraction
  };

  Double_t xmin[nDim] = {
    0.,     // pvalue
    -280.5, // residuals
    -0.5,   // ADC_fit
    -0.5,   // ADC_mes
    -5,     // ADC_cluster
    -0.5,   // nSamples
    -0.7,   // Asymm
    -0.5,   // Wire
    -2,     // NonBending
    0.      // fraction
  };

  Double_t xmax[nDim] = {
    1.,      // pvalue
    280.5,   // residuals
    10000.5, // ADC_fit
    10000.5, // ADC_mes
    50005,   // ADC_cluster
    300.5,   // nSamples
    0.7,     // Asymm
    2.5,     // Wire
    2,       // Bending
    1.       // fraction
  };

  TString name = Form("MultiResolutionPreCluster%s", extension);
  TString title = "10D Sparse Histograms for Pre-Cluster ";

  THnSparseD* hSparse = new THnSparseD(name, title, nDim, nbins, xmin, xmax);

  const char* axisTitles[nDim] = {
    "pvalue", "residuals", "ADC_fit", "ADC_mes", "ADC_cluster",
    "nSamples", "Asymm", "Wire", "Cathode", "fraction"};

  for (Int_t i = 0; i < nDim; ++i) {
    hSparse->GetAxis(i)->SetTitle(axisTitles[i]);
  }

  return hSparse;
}
//_________________________________________________________________________________________________
// fill THnSparse (10 axes) with the preclusters characteristics
// plane = 1. (bending) or -1. (non-bending)
// the vector "parameters" is a 6 size vector which is defined as :
// parameters = {sqrt(Qb_tot * Qnb_tot), (NB - B)/(NB + B), distance closest wire, pvalue, chargeB, chargeNB}
void FillResolutionInfo(const Digit& digit, double ADC_fit, double plane,
                        const std::vector<double>& parameters, THnSparseD* h)
{
  Double_t position = -1.;
  if ((std::abs(parameters[2]) < 0.015)) { //"top"
    position = 0.;
  } else if (std::abs(parameters[2]) > 0.075) { //"between"
    position = 2.;
  } else if ((std::abs(parameters[2]) > 0.015) && (std::abs(parameters[2]) < 0.075)) { //"crossover"
    position = 1.;
  }

  Double_t ADC_mes = digit.getADC();
  Double_t residuals = ADC_mes - ADC_fit;
  Double_t nSamples = digit.getNofSamples();
  Double_t fraction = (plane > 0.) ? ADC_mes / parameters[4] : ADC_mes / parameters[5];

  Double_t values[10] = {
    parameters[3], // pvalue
    residuals,     // residuals
    ADC_fit,       // ADC_fit
    ADC_mes,       // ADC_mes
    parameters[0], // ADC_cluster
    nSamples,      // nSamples
    parameters[1], // Asymm
    position,      // Wire
    plane,         // Cathode
    fraction       // fraction
  };

  h->Fill(values);
}
//_________________________________________________________________________________________________
// create the THnSparse (8 axis) for k3 studies
THnSparseD* CreatePreClusterInfoMULTIK3(const char* extension = "")
{
  const Int_t nDim = 8;

  Int_t nbins[nDim] = {
    1000, // pvalue
    40,   // k3x
    40,   // k3y
    5001, // ADC_cluster
    51,   // p
    62,   // phi
    280,  // Asymm
    3     // Wire
  };

  Double_t xmin[nDim] = {
    0.,    // pvalue
    0.,    // k3x
    0.,    // k3y
    -5,    // ADC_cluster
    -0.5,  // p
    -15.5, // phi
    -0.7,  // Asymm
    -0.5   // Wire
  };

  Double_t xmax[nDim] = {
    1.,    // pvalue
    1.,    // k3x
    1.,    // k3y
    50005, // ADC_cluster
    50.5,  // p
    15.5,  // phi
    0.7,   // Asymm
    2.5    // Wire
  };

  TString name = Form("MultiK3PreCluster%s", extension);
  TString title = "9D Sparse Histograms for Pre-Cluster ";

  THnSparseD* hSparse = new THnSparseD(name, title, nDim, nbins, xmin, xmax);

  const char* axisTitles[nDim] = {
    "pvalue", "k3x", "k3y", "ADC_cluster", "p", "phi", "Asymm", "Wire"};

  for (Int_t i = 0; i < nDim; ++i) {
    hSparse->GetAxis(i)->SetTitle(axisTitles[i]);
  }

  return hSparse;
}
//_________________________________________________________________________________________________
// fill THnSparse (8 axis) for k3 studies
// the vector parameters is a 12 size vector which is defined as :
// parameters = {X, Y, K3x, K3y, Qb_tot, Qnb_tot, sqrt(Qb_tot * Qnb_tot), (NB - B)/(NB + B), distance closest wire, pvalue, track angle, track momentum}
void FillK3Info(const std::vector<double>& parameters, THnSparseD* h)
{
  Double_t position = -1.;
  if ((std::abs(parameters[8]) < 0.015)) { //"top"
    position = 0.;
  } else if (std::abs(parameters[8]) > 0.075) { //"between"
    position = 2.;
  } else if ((std::abs(parameters[8]) > 0.015) && (std::abs(parameters[8]) < 0.075)) { //"crossover"
    position = 1.;
  }

  Double_t values[8] = {
    parameters[9],  // pvalue
    parameters[2],  // k3x
    parameters[3],  // k3y
    parameters[6],  // ADC_cluster
    parameters[11], // p
    parameters[10], // phi
    parameters[7],  // Asymm
    position        // Wire
  };

  h->Fill(values);
}

//_________________________________________________________________________________________________
// extract the resolution (std) of the residuals distribution for different ADC which are defined as : residuals = ADC - ADC_fit
// project the TH2D into a corresponding axis ->
// on X : to chose a bin of ADC which as enough statistic to extract a correct resolution
// on Y : to extract the resolution (std) of the residuals distribution of the corresponding ADC binning
// the fit for the std is done 3 times because of the shape of the residuals distribution (see current studies)
// use auto binning (size vary with a define statistic) or harcoded binning (pre defined binning)
// save results in TList with the fit properties (i.e. : chi2, std, mean, ...)
void Resolution(TList*& list, TH2D* hist2D, int statistics, bool auto_bin)
{
  // default digit range value : 20 - 10000 ADC
  int start = 20, end = 10000;
  std::vector<double> wavg_charge;

  TH1D* projX = hist2D->ProjectionX();
  int binStart = projX->FindBin(start);
  int binEnd = projX->FindBin(end);
  int bin_i = binStart;
  double wavg = 0.;

  std::vector<std::pair<int, int>> intervals;

  if (auto_bin) {
    for (int bin_j = binStart; bin_j < binEnd; bin_j++) {

      int integral = projX->Integral(bin_i, bin_j);
      wavg += projX->GetBinCenter(bin_j) * projX->GetBinContent(bin_j);

      if (integral > statistics) {
        wavg_charge.push_back(wavg / integral);
        intervals.push_back(std::make_pair(bin_i, bin_j));
        bin_i = bin_j + 1;
        wavg = 0.; // reset wavg to 0
      }
    }
  } else {
    // hardcoded binning
    std::vector<std::pair<int, int>> charge_bin;

    for (int i = start; i <= 500; i += 1) {
      charge_bin.push_back({i, i});
    }
    for (int i = charge_bin.back().second + 1; i <= 1000; i += 4) {
      charge_bin.push_back({i, i + 3});
    }
    for (int i = charge_bin.back().second + 1; i <= 2500; i += 10) {
      charge_bin.push_back({i, i + 9});
    }
    for (int i = charge_bin.back().second + 1; i <= 5000; i += 16) {
      charge_bin.push_back({i, i + 15});
    }
    for (int i = charge_bin.back().second + 1; i <= end - 25; i += 26) {
      charge_bin.push_back({i, i + 25});
    }

    for (auto interval : charge_bin) {
      int bin_i = projX->FindBin(interval.first);
      int bin_j = projX->FindBin(interval.second);
      int integral = projX->Integral(bin_i, bin_j);

      if (integral > statistics) {
        double wavg = (interval.first + interval.second) / 2.;
        wavg_charge.push_back(wavg);
        intervals.push_back(std::make_pair(bin_i, bin_j));
      }
    }
  }
  // auto rebin after first gaussian fit
  int index = 0;
  for (auto I : intervals) {
    // we save the charge interval into the name of the 1D projection
    TH1D* projY = hist2D->ProjectionY(Form("projY_%s_%d_%d", hist2D->GetName(), static_cast<int>(projX->GetBinCenter(I.first)), static_cast<int>(projX->GetBinCenter(I.second))), I.first, I.second);

    // to guide the gaussian fit
    double initial_sigma = sqrt(wavg_charge[index]);
    double min = -1.8 * initial_sigma;
    double max = 1.8 * initial_sigma;

    TF1* fit = new TF1("fit", "gaus", min, max);
    fit->SetParameter(0, projY->GetMaximum());
    fit->SetParameter(1, 0.0);
    fit->SetParameter(2, 0.5 * initial_sigma);
    fit->SetParLimits(1, min, max);
    projY->Fit(fit, "RQ");

    double sigma1 = std::abs(fit->GetParameter(2));
    double mean1 = fit->GetParameter(1);

    //---------- SECOND FIT ----------
    double range2 = 1.5 * sigma1;
    fit->SetRange(mean1 - range2, mean1 + range2);
    fit->SetParLimits(1, mean1 - range2, mean1 + range2);
    fit->SetParameters(fit->GetParameter(0), mean1, sigma1);
    projY->Fit(fit, "RQ");

    double sigma2 = std::abs(fit->GetParameter(2));
    double mean2 = fit->GetParameter(1);

    //---------- THIRD FIT ----------
    double range3 = 1.5 * sigma2;
    fit->SetRange(mean2 - range3, mean2 + range3);
    fit->SetParLimits(1, mean2 - range3, mean2 + range3);
    fit->SetParameters(fit->GetParameter(0), mean2, sigma2);
    projY->Fit(fit, "RQ");

    double par[3];
    fit->GetParameters(par);
    list->Add(projY);
    delete fit;
    index++;
  }
  delete projX;
}

#endif
