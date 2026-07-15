#ifndef RESOLUTIONUTILS_H_
#define RESOLUTIONUTILS_H_

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <Math/ProbFunc.h>
#include <TAxis.h>
#include <TF1.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
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
    8000,  // residuals
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
    -2000., // residuals
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
    2000.,   // residuals
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
// fit function of the ADC - ADC_fit residuals that convolutes pad charge resolution with ADC threshold
// p[0] = normalization
// p[1] = systematic shift
// p[2] = resolution
// p[3] = mean value of the ADC interval being fitted (should be fixed)
double ADCResolutionWithThreshold(double* x, double* p)
{
  // threshold parameters
  static const double thresholdMean = 22.2;
  static const double thresholdSigma = 2.8;

  // pad charge resolution
  double resolution = TMath::Gaus(*x, p[1], p[2], true);

  // probability that the measured value (= fit value + residual) passes the gaussian threshold
  double threshold = ROOT::Math::gaussian_cdf(p[3] + (*x), thresholdSigma, thresholdMean);

  return p[0] * resolution * threshold;
}

//_________________________________________________________________________________________________
// extract the resolution (std) of the residuals distribution for different ADC which are defined as : residuals = ADC - ADC_fit
// project the TH2D into a corresponding axis ->
// on X : to chose a bin of ADC which has enough statistics to extract a correct resolution
// on Y : to extract the resolution (std) of the residuals distribution of the corresponding ADC binning
// the fit for the std is done several times because of the shape of the residuals distribution (see current studies)
// use auto binning (size vary with a define statistic) or harcoded binning (pre defined binning)
// save results in TList with the fit properties (i.e. : chi2, std, mean, ...)
void Resolution(TList*& list, TH2D* hist2D, int minStat, bool auto_bin)
{
  // default digit range value : 20 - 10000 ADC
  static int start = 20;
  static int end = 10000;

  TH1D* projX = hist2D->ProjectionX();
  int binStart = projX->FindBin(start);
  int binEnd = projX->FindBin(end);

  std::vector<std::pair<int, int>> intervals;

  if (auto_bin) {
    int bin_i = binStart;
    for (int bin_j = binStart; bin_j < binEnd; bin_j++) {

      int integral = projX->Integral(bin_i, bin_j);

      if (integral > minStat) {
        intervals.push_back(std::make_pair(bin_i, bin_j));
        bin_i = bin_j + 1;
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

    for (const auto& [min, max] : charge_bin) {
      int bin_i = projX->FindBin(min);
      int bin_j = projX->FindBin(max);
      int integral = projX->Integral(bin_i, bin_j);

      if (integral > minStat) {
        intervals.push_back(std::make_pair(bin_i, bin_j));
      }
    }
  }

  int index = 0;
  for (const auto& [min, max] : intervals) {

    // we save the charge interval into the name of the 1D projection
    auto minADC = static_cast<int>(std::round(projX->GetBinCenter(min)));
    auto maxADC = static_cast<int>(std::round(projX->GetBinCenter(max)));
    TH1D* projY = hist2D->ProjectionY(Form("projY_%s_%d_%d", hist2D->GetName(), minADC, maxADC), min, max);

    // first gaussian fit to define the binning and set the next fit range around the peak
    static TF1* fit1 = new TF1("fit1", "gausn");
    double norm = projY->GetEntries() * projY->GetBinWidth(1);
    double mean = projY->GetMean();
    double stdDev = projY->GetStdDev();
    fit1->SetParameters(norm, mean, 0.5 * stdDev);
    fit1->SetParLimits(0, 0., 2. * norm);
    fit1->SetParLimits(1, std::max(mean - 2. * stdDev, projY->GetXaxis()->GetXmin()),
                       std::min(mean + 2. * stdDev, projY->GetXaxis()->GetXmax()));
    fit1->SetParLimits(2, 0., 3. * stdDev);
    int status = projY->Fit(fit1, "BQ");
    if (status != 0) {
      printf("%s: first fit first attempt failed with status %d (result = [%f, %f, %f])\n",
             projY->GetName(), status, fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2));
      fit1->SetParameters(norm, 0., stdDev);
      fit1->SetParLimits(2, 0., 2. * stdDev);
      status = projY->Fit(fit1, "BQ", "", -5. * stdDev, 5. * stdDev);
      if (status != 0) {
        printf("%s: first fit second attempt failed with status %d (result = [%f, %f, %f])\n",
               projY->GetName(), status, fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2));
      }
    }

    // rebin
    double mean1 = fit1->GetParameter(1);
    double sigma1 = std::abs(fit1->GetParameter(2));
    if (status != 0 || stdDev < sigma1) {
      sigma1 = stdDev;
      mean1 = mean;
    }
    int rebin = std::round(sigma1 / 10. / projY->GetXaxis()->GetBinWidth(1));
    if (rebin > 0) {
      gErrorIgnoreLevel = kWarning + 1; // silence warning "ngroup=xxx is not an exact divider of nbins=yyy"
      projY->Rebin(rebin);
      gErrorIgnoreLevel = 0;
    }

    // second fit repeated n times to refine the fit range around the peak
    for (int i = 0; i < 5; ++i) {
      static TF1* fit2 = new TF1("fit2", "gausn");
      double xMin = mean1 - 2. * sigma1;
      double xMax = mean1 + 2. * sigma1;
      fit2->SetParameters(norm, mean1, sigma1);
      status = projY->Fit(fit2, "BQ", "", xMin, xMax);
      if (status != 0) {
        printf("%s: second fit first attempt failed with status %d (result = [%f, %f, %f])\n",
               projY->GetName(), status, fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(2));
        fit2->SetParameters(norm, 0., 0.5 * sigma1);
        status = projY->Fit(fit2, "BQ", "", xMin, xMax);
        if (status != 0) {
          printf("%s: second fit second attempt failed with status %d (result = [%f, %f, %f])\n",
                 projY->GetName(), status, fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(2));
          break;
        }
      }
      mean1 = fit2->GetParameter(1);
      sigma1 = std::abs(fit2->GetParameter(2));
    }

    // third fit to extract the resolution
    static TF1* fit3 = new TF1("fit", ADCResolutionWithThreshold, -1000., 1000., 4);
    double xMin = mean1 - 2. * sigma1;
    double xMax = mean1 + 2. * sigma1;
    double meanADC = 0.5 * (minADC + maxADC);
    fit3->SetParameters(norm, mean1, sigma1, meanADC);
    fit3->FixParameter(3, meanADC);
    status = projY->Fit(fit3, "BQ", "", xMin, xMax);
    if (status != 0) {
      printf("%s: third fit first attempt failed with status %d (result = [%f, %f, %f])\n",
             projY->GetName(), status, fit3->GetParameter(0), fit3->GetParameter(1), fit3->GetParameter(2));
      fit3->SetParameters(norm, 0., 0.5 * sigma1, meanADC);
      status = projY->Fit(fit3, "BQ", "", xMin, xMax);
      if (status != 0) {
        printf("%s: third fit second attempt failed with status %d (result = [%f, %f, %f])\n",
               projY->GetName(), status, fit3->GetParameter(0), fit3->GetParameter(1), fit3->GetParameter(2));
      }
    }

    list->Add(projY);
    index++;
  }
  delete projX;
}

#endif
