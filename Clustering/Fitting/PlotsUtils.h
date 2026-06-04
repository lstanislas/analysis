#ifndef PLOTSUTILS_H_
#define PLOTSUTILS_H_

#include <cmath>
#include <vector>
#include <functional>

#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <TH2D.h>
#include <TH1D.h>

// define all used histograms for the analysis
std::vector<TH2D*> h2chi2_ndf;
std::vector<TH1D*> hprob;
std::vector<TH1D*> hk3x;
std::vector<TH1D*> hk3y;
std::vector<TH1D*> hAsymm;

//_________________________________________________________________________________________________
// premade histograms for the analysis
void DelHist();

void LoadHist()
{
  DelHist(); // clear any previous histograms from a prior run in the same ROOT session

  // station index
  int station[3] = {1, 2, 345};

  // chi2 vs ndf
  for (int i = 0; i < 3; ++i) {
    h2chi2_ndf.push_back(new TH2D(Form("h2chi2_ndf_%d", i), Form("#chi^{2} vs ndf (St.%d)", station[i]), 16, -0.5, 15.5, 1000, 0, 50));
    h2chi2_ndf[i]->SetDirectory(0);
    h2chi2_ndf[i]->GetXaxis()->SetTitle("ndf");
    h2chi2_ndf[i]->GetYaxis()->SetTitle("#chi^{2}");
  }
  // p-value distribution
  for (int i = 0; i < 3; ++i) {
    hprob.push_back(new TH1D(Form("hprob_%d", i), Form("P-Value (St.%d)", station[i]), 1000, 0, 1));
    hprob[i]->SetDirectory(0);
    hprob[i]->GetXaxis()->SetTitle("p-value");
  }
  // k3x distribution
  for (int i = 0; i < 3; ++i) {
    hk3x.push_back(new TH1D(Form("hk3x_%d", i), Form("K_{3X} (St.%d)", station[i]), 10000, 0, 1));
    hk3x[i]->SetDirectory(0);
    hk3x[i]->GetXaxis()->SetTitle("K_{3X}");
  }
  // k3y distribution
  for (int i = 0; i < 3; ++i) {
    hk3y.push_back(new TH1D(Form("hk3y_%d", i), Form("K_{3Y} (St.%d)", station[i]), 10000, 0, 1));
    hk3y[i]->SetDirectory(0);
    hk3y[i]->GetXaxis()->SetTitle("K_{3Y}");
  }
  // Asymmetry
  for (int i = 0; i < 6; ++i) {
    hAsymm.push_back(new TH1D(Form("hAsymm_%d", i), Form("Asymmetry bending - nbending (St.%d)", station[i / 2]), 700, -0.7, 0.7));
    hAsymm[i]->SetDirectory(0);
    hAsymm[i]->GetXaxis()->SetTitle("(NB-B)/(NB+B)");
  }
}

//_________________________________________________________________________________________________
// delete properly previous histograms
void DelHist()
{

  for (auto& hist : h2chi2_ndf)
    delete hist;
  h2chi2_ndf.clear();

  for (auto& hist : hprob)
    delete hist;
  hprob.clear();

  for (auto& hist : hk3x)
    delete hist;
  hk3x.clear();

  for (auto& hist : hk3y)
    delete hist;
  hk3y.clear();

  for (auto& hist : hAsymm)
    delete hist;
  hAsymm.clear();
}

//_________________________________________________________________________________________________
// default name, size, title, ...
void FillInfoGraphErrAsymm(TGraphAsymmErrors* graph, const TString& titleX, const TString& titleY, const TString& Station, const TString& Cathode, bool primeFile)
{
  TString file = primeFile ? "(1)" : "(2)";
  TString title = Form("%s vs %s [%s] [%s] %s",
                       titleY.Data(), titleX.Data(),
                       Station.Data(), Cathode.Data(),
                       file.Data());

  graph->SetTitle(title);
  graph->SetName(title);
  graph->GetXaxis()->SetLimits(0, 10001);
  graph->GetXaxis()->SetTitle(titleX);
  graph->GetYaxis()->SetTitle(titleY);
  graph->SetMarkerColor((primeFile) ? kRed : kBlue);
  graph->SetMarkerStyle(47);
  graph->SetMinimum(0);
  graph->SetMaximum(200);
}
//_________________________________________________________________________________________________
// default name, size, title, ...
void FillInfoGraph(TGraph* graph, const TString& titleX, const TString& titleY, const TString& Station, const TString& Cathode, bool primeFile)
{
  TString file = primeFile ? "(1)" : "(2)";
  TString title = Form("%s vs %s [%s] [%s] %s",
                       titleY.Data(), titleX.Data(),
                       Station.Data(), Cathode.Data(),
                       file.Data());

  graph->SetTitle(title);
  graph->SetName(title);
  graph->GetXaxis()->SetLimits(0, 10001);
  graph->GetXaxis()->SetTitle(titleX);
  graph->GetYaxis()->SetTitle(titleY);
  graph->SetMarkerColor((primeFile) ? kRed : kBlue);
  graph->SetMarkerStyle(47);
  if (titleY == "#mu") {
    graph->SetMinimum(-20);
    graph->SetMaximum(20);
  } else {
    graph->SetMinimum(0);
    graph->SetMaximum(200);
  }
}
//_________________________________________________________________________________________________
// default name, size, title, ...
void FillInfoHist(TH1D* h, const TString& titleX, const TString& titleY, const TString& Station, const TString& Cathode, bool primeFile, bool normalize = false)
{
  TString oldName = h->GetName();
  TString file = primeFile ? "(1)" : "(2)";
  TString title = Form("%s vs %s %s [%s] %s %s",
                       titleY.Data(), titleX.Data(),
                       oldName.Data(), Station.Data(), Cathode.Data(), file.Data());

  if (normalize && (h->Integral() != 0)) {
    h->Scale(1.0 / h->Integral());
  }

  h->SetTitle(title);
  h->SetName(title);

  h->GetXaxis()->SetTitle(titleX);
  h->GetYaxis()->SetTitle(titleY);
  h->SetLineColor((primeFile) ? kRed : kBlue);
}
//_________________________________________________________________________________________________
// plot a defined function
TF1* plotNoise(std::string sigma, double alpha, double gamma = 0., double beta = 0., bool Asymm = false)
{
  TF1* Func = nullptr;

  if (sigma == "MC") {
    if (!Asymm) {
      if (!Func) {
        Func = new TF1("func", [](double* x, double* par) {
                              double charge = x[0];
                              double result = (std::pow(charge / par[1], 1. / par[2]) + par[0]);
                              return (0.5 * (sqrt(std::round(result)) + par[4])); }, 0, 10000, 4);
      }
      double signalParam[3] = {14., 13., 1.5};
      Func->SetParameters(signalParam[0], signalParam[1], signalParam[2], alpha);
      Func->SetLineColor(4);
      Func->SetLineWidth(3);
      Func->SetLineStyle(1);
      Func->SetNpx(1000);
    } else {
      if (!Func) {
        Func = new TF1("func", [](double* x, double* par) {
                            double charge = x[0];
                            double result = (std::pow(charge / par[1], 1. / par[2]) + par[0]);
                            return sqrt(0.25 * (std::round(result) + par[3]) + 0.25 * charge * charge * (TMath::Exp(8 * par[4] * par[4]) - TMath::Exp(4 * par[4] * par[4]))); }, 0, 10000, 5);
      }
      double sigmaAsymm = gamma * 0.055;
      double signalParam[3] = {14., 13., 1.5};
      Func->SetParameters(signalParam[0], signalParam[1], signalParam[2], alpha, sigmaAsymm);
      Func->SetLineColor(4);
      Func->SetLineWidth(3);
      Func->SetLineStyle(1);
      Func->SetNpx(1000);
    }
  } else if (sigma == "sADC") {
    if (!Asymm) {
      if (!Func) {
        Func = new TF1("func", [](double* x, double* par) { return par[0] * sqrt(x[0]); }, 0, 10000, 1);
      }

      alpha = (alpha < 0.001) ? 1. : alpha;
      Func->SetParameters(alpha);
      Func->SetLineColor(4);
      Func->SetLineWidth(3);
      Func->SetLineStyle(1);
      Func->SetNpx(1000);
    } else {
      if (!Func) {
        Func = new TF1("func", [](double* x, double* par) {
                            double charge = x[0];
                            double intrinsic = par[0] * sqrt(charge);
                            return sqrt(par[0] * par[0] * sqrt(charge) * sqrt(charge) + 0.25 * charge * charge * (TMath::Exp(8 * par[1] * par[1]) - TMath::Exp(4 * par[1] * par[1]))); }, 0, 10000, 2);
      }
      double sigmaAsymm = gamma * 0.055;
      alpha = (alpha < 0.001) ? 1. : alpha;
      Func->SetParameters(alpha, sigmaAsymm);
      Func->SetLineColor(4);
      Func->SetLineWidth(3);
      Func->SetLineStyle(1);
      Func->SetNpx(1000);
    }
  } else if (sigma == "MULT") {
    if (!Asymm) {
      if (!Func) {
        Func = new TF1("func", [](double* x, double* par) {
          double charge = x[0];
          return par[0] * std::sqrt(charge) + par[1] * charge + par[2] * charge * std::sqrt(charge); }, 0, 10000, 3);
      }
      Func->SetParameters(alpha, beta, gamma);
      Func->SetLineColor(4);
      Func->SetLineWidth(3);
      Func->SetLineStyle(1);
      Func->SetNpx(1000);
    }
  } else {
    std::cerr << "Unknown sigma type!" << std::endl;
    return nullptr;
  }
  return Func;
}

//_________________________________________________________________________________________________
// take a vector of TH1D, divide a canvas into 1xN or 2xN then the i-th and i+1-th element of the vector TH1D
// are on the same plot corresponding to the index i
void plotSAME(std::vector<TH1D*> hist, const char* name, const char* title, bool all = false)
{
  TCanvas* c = new TCanvas(name, title, 800, 600);
  int N = hist.size();
  int COL = 1;
  if (all) {
    COL = 2;
  }

  c->Divide(COL, (N / 2 + COL - 1) / COL);
  for (int i = 1; i < (N / 2 + 1); ++i) {
    TVirtualPad* pad = c->cd(i);
    if (!pad) {
      std::cerr << "plotSAME: cd(" << i << ") returned null\n";
      continue;
    }
    pad->SetLogy();
    hist[2 * (i - 1)]->SetLineColor(2);
    hist[2 * (i - 1)]->Draw("HIST");
    hist[2 * (i - 1) + 1]->Draw("HIST SAME");
  }
  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________
// take a vector of TH1D, divide a canvas into 1xN or 2xN and plot them
void plot1D(std::vector<TH1D*> hist, const char* name, const char* title, bool all = false)
{
  TCanvas* c = new TCanvas(name, title, 800, 600);
  int N = hist.size();
  int COL = 1;
  if (all) {
    COL = 2;
  }

  c->Divide(COL, (N + COL - 1) / COL);
  for (int i = 1; i < N + 1; ++i) {
    TVirtualPad* pad = c->cd(i);
    if (!pad) {
      std::cerr << "plot1D: cd(" << i << ") returned null\n";
      continue;
    }
    pad->SetLogy();
    hist[i - 1]->Draw("HIST");
  }
  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________
// take a vector of TH2D, divide a canvas into 1xN or 2xN and plot them
void plot2D(std::vector<TH2D*> hist, const char* name, const char* title, bool all = false)
{
  TCanvas* c = new TCanvas(name, title, 800, 600);
  int N = hist.size();
  int COL = 1;
  if (all) {
    COL = 2;
  }

  c->Divide(COL, (N + COL - 1) / COL);
  for (int i = 1; i < N + 1; ++i) {
    TVirtualPad* pad = c->cd(i);
    if (!pad) {
      std::cerr << "plot2D: cd(" << i << ") returned null\n";
      continue;
    }
    pad->SetLogz();
    gStyle->SetPalette(kRainBow);
    hist[i - 1]->Draw("COLZ");
  }
  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________
// build the theoretical noise TF1 from a noise model string: "MULT_XpX_XpX_XpX", "sADC_XpX", "MC_XpX"
// ('p' encodes the decimal point, same convention as BuildToyMC)
TF1* noiseFromString(const std::string& noise)
{
  auto extractOne = [](const std::string& s, const std::string& pat) -> double {
    size_t pos = s.find(pat);
    if (pos == std::string::npos)
      return 1.0;
    std::string tok = s.substr(pos + pat.size());
    size_t next = tok.find('_');
    if (next != std::string::npos)
      tok = tok.substr(0, next);
    std::replace(tok.begin(), tok.end(), 'p', '.');
    try {
      return std::stod(tok);
    } catch (...) {
      return 1.0;
    }
  };

  if (noise.find("MULT_") != std::string::npos) {
    std::string s = noise.substr(noise.find("MULT_") + 5);
    std::replace(s.begin(), s.end(), 'p', '.');
    size_t p1 = s.find('_'),
           p2 = (p1 != std::string::npos) ? s.find('_', p1 + 1) : std::string::npos;
    double alpha = 1., beta = 0., gamma = 0.;
    if (p1 != std::string::npos && p2 != std::string::npos) {
      try {
        alpha = std::stod(s.substr(0, p1));
        beta = std::stod(s.substr(p1 + 1, p2 - p1 - 1));
        gamma = std::stod(s.substr(p2 + 1)); // rest of string, no trailing '_' needed
      } catch (...) {
      }
    }
    return plotNoise("MULT", alpha, gamma, beta);
  } else if (noise.find("sADC_") != std::string::npos) {
    return plotNoise("sADC", extractOne(noise, "sADC_"));
  } else if (noise.find("MC_") != std::string::npos) {
    return plotNoise("MC", extractOne(noise, "MC_"));
  }
  return plotNoise("sADC", 1.0);
}

//_________________________________________________________________________________________________
// take two vectors of TGraphAsymmErrors object corresponding to a certain station, divide a canvas in two (B / NB),
// plots the graphs along with a defined function from plotNoise
void tGraphErrAsymm(std::vector<TGraphAsymmErrors*>& graphs1, std::vector<TGraphAsymmErrors*>& graphs2, const std::string& name, int station, const std::string& noise = "")
{

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1200, 800);
  c->Divide(1, 2);

  auto func = noiseFromString(noise);

  for (int i = 0; i < 2; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();

    auto* g1 = graphs1[2 * station + i];
    auto* g2 = graphs2[2 * station + i];

    g1->Draw("APE");
    g2->Draw("PE SAME");
    func->Draw("L SAME");

    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(g1, g1->GetTitle(), "p");
    leg->AddEntry(g2, g2->GetTitle(), "p");
    leg->AddEntry(func, "ToyMC th. noise", "l");
    leg->Draw();
  }

  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________
// take two vectors of TGraph object corresponding to a certain station, divide a canvas in four (B / NB) x (DATA / TMC), plots the graphs
void tGraph(std::vector<TGraph*>& graphs1, std::vector<TGraph*>& graphs2, const std::string& name, int station)
{

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1200, 800);
  c->Divide(2, 2);

  TGraph* plots[4] = {
    graphs1[2 * station], graphs2[2 * station],
    graphs1[2 * station + 1], graphs2[2 * station + 1]};

  for (int i = 0; i < 4; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();

    plots[i]->Draw("AP");

    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(plots[i], plots[i]->GetTitle(), "p");
    leg->Draw();
  }

  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________
// take two vectors of TH1D object corresponding to a certain station and a certain cathode, divide a canvas in
// eight (corresponding to a choice of ADC binning), plots the graphs
void tHist(std::vector<TH1D*>& hists1, std::vector<TH1D*>& hists2, const std::string& name, int station, int cathode)
{
  std::vector<TH1D*> subh1(hists1.begin() + 8 * (2 * station + cathode), hists1.begin() + 8 * (2 * station + cathode + 1));
  std::vector<TH1D*> subh2(hists2.begin() + 8 * (2 * station + cathode), hists2.begin() + 8 * (2 * station + cathode + 1));

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1200, 800);
  c->Divide(2, 4);

  for (int i = 0; i < 8; ++i) {
    c->cd(i + 1);
    gPad->SetGrid();
    gPad->SetLogy();

    double max1 = subh1[i]->GetMaximum();
    double max2 = subh2[i]->GetMaximum();
    double ymax = std::max(max1, max2) * 1.5;

    subh1[i]->GetXaxis()->SetRangeUser(-50, 50);
    subh1[i]->SetMaximum(ymax);
    subh1[i]->Draw("HIST");
    subh2[i]->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    TString title1 = Form("%s (N=%d)", subh1[i]->GetTitle(), (int)subh1[i]->GetEntries());
    TString title2 = Form("%s (N=%d)", subh2[i]->GetTitle(), (int)subh2[i]->GetEntries());
    leg->AddEntry(subh1[i], title1.Data(), "l");
    leg->AddEntry(subh2[i], title2.Data(), "l");
    leg->Draw();
  }

  c->Update();
  c->Draw();
}

//_________________________________________________________________________________________________
// take two vectors of TGraphAsymmErrors object corresponding to a station with a parameter alpha (to define a function with plotFunction),
// divide a canvas in two (B / NB), plots the ratio of the y-data points between the two TGraphAsymmErrors object along with the corresponding errors
void tRatio(std::vector<TGraphAsymmErrors*>& graphs1, std::vector<TGraphAsymmErrors*>& graphs2, const std::string& name, int station, const std::string& noise)
{

  // parse noise model string (same format as BuildToyMC: "MULT_XpX_XpX_XpX", "sADC_XpX", "MC_XpX")
  auto extractOne = [](const std::string& s, const std::string& pat) -> double {
    size_t pos = s.find(pat);
    if (pos == std::string::npos)
      return 1.0;
    std::string tok = s.substr(pos + pat.size());
    size_t next = tok.find('_');
    if (next != std::string::npos)
      tok = tok.substr(0, next);
    std::replace(tok.begin(), tok.end(), 'p', '.');
    try {
      return std::stod(tok);
    } catch (...) {
      return 1.0;
    }
  };

  TF1* func = noiseFromString(noise);
  double alpha = 1.0, beta = 0.0, gamma = 0.0;
  std::function<double(double)> lambda;

  if (noise.find("MULT_") != std::string::npos) {
    std::string s = noise.substr(noise.find("MULT_") + 5);
    std::replace(s.begin(), s.end(), 'p', '.');
    size_t p1 = s.find('_'),
           p2 = (p1 != std::string::npos) ? s.find('_', p1 + 1) : std::string::npos;
    if (p1 != std::string::npos && p2 != std::string::npos) {
      try {
        alpha = std::stod(s.substr(0, p1));
        beta = std::stod(s.substr(p1 + 1, p2 - p1 - 1));
        gamma = std::stod(s.substr(p2 + 1)); // rest of string, no trailing '_' needed
      } catch (...) {
      }
    }
    lambda = [alpha, beta, gamma](double charge) {
      return alpha * std::sqrt(charge) + beta * charge + gamma * charge * std::sqrt(charge);
    };
    std::cout << "Using MULT_ equation with alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << std::endl;
  } else if (noise.find("sADC_") != std::string::npos) {
    alpha = extractOne(noise, "sADC_");
    lambda = [alpha](double charge) { return alpha * std::sqrt(charge); };
    std::cout << "Using sADC_ equation with alpha=" << alpha << std::endl;
  } else if (noise.find("MC_") != std::string::npos) {
    alpha = extractOne(noise, "MC_");
    lambda = [alpha](double charge) {
      double result = std::pow(charge / 13., 1. / 1.5) + 14.;
      return 0.5 * (std::sqrt(std::round(result)) + alpha);
    };
    std::cout << "Using MC_ equation with alpha=" << alpha << std::endl;
  } else {
    lambda = [](double charge) { return std::sqrt(charge); };
    std::cout << "DEFAULT: Using sADC_ equation with alpha=1" << std::endl;
  }

  TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 1000, 800);
  c->Divide(1, 2);
  for (int i = 0; i < 2; ++i) {

    c->cd(i + 1);
    gPad->SetGrid();

    TGraphAsymmErrors* gi = graphs1[2 * station + i];
    TGraphAsymmErrors* gj = graphs2[2 * station + i];
    TGraphAsymmErrors* ratio = new TGraphAsymmErrors();        // draw (#sigma_{th}(#sigma_{DATA} / #sigma_{TMC})
    TGraphAsymmErrors* ratio_simple = new TGraphAsymmErrors(); // draw (#sigma_{DATA} / #sigma_{TMC})

    TString Cathode = (i == 0) ? "Bending" : "NonBending";
    TString title = Form("#sigma_{th}(#sigma_{DATA} / #sigma_{TMC}) - %s", Cathode.Data());
    TString title_simple = Form("(#sigma_{DATA} / #sigma_{TMC}) - %s", Cathode.Data());

    ratio->SetTitle(title);
    ratio_simple->SetTitle(title_simple);

    ratio->GetXaxis()->SetTitle("ADC");
    ratio->GetYaxis()->SetTitle("#sigma");
    ratio->SetMarkerStyle(2);
    ratio_simple->SetMarkerStyle(2);
    ratio_simple->SetMarkerColor(kGreen);

    // since it is not said that the two TGraphAsymmErrors object should have the same number of elements with the same x values
    // we take as reference the first TGraphAsymmErrors object and loop over its elements (i index)
    // (in the i-loop) we loop over the elements of second TGraphAsymmErrors object (j index)
    // if the difference between | x_j - x_i | is lower than 10^-6 we do the ratio of y_i / y_j
    // futhermore, because of how we construct the two TGraphAsymmErrors object
    // we are sure that if a couple (x_j, x_i) respect the previous condition then
    // x_j (resp. x_i) cannot respect the condtion with the elements other than x_i (resp. x_j)
    int Ni = gi->GetN();
    for (int k = 0; k < Ni; ++k) {
      double xi, yi;
      gi->GetPoint(k, xi, yi);
      if (yi == 0) // check if the y value is not zero
        continue;
      double yi_err = 0.5 * (gi->GetErrorYlow(k) + gi->GetErrorYhigh(k)); // low and high errors are the same by construction

      int Nj = gj->GetN();
      for (int l = 0; l < Nj; ++l) {
        double xj, yj;
        gj->GetPoint(l, xj, yj);
        if (std::abs(xi - xj) < 1e-6 && yj != 0) {
          double yj_err = 0.5 * (gj->GetErrorYlow(l) + gj->GetErrorYhigh(l));
          double r = yi / yj;
          double rel_err = std::sqrt(std::pow(yi_err / yi, 2) + std::pow(yj_err / yj, 2)); // propagate errors
          double r_err = r * rel_err * lambda(xi);
          double r_err_simple = r * rel_err;

          ratio->AddPoint(xi, (r * lambda(xi))); // ratio #sigma_{th}(#sigma_{DATA} / #sigma_{TMC})
          ratio_simple->AddPoint(xi, r);         // second ratio (#sigma_{DATA} / #sigma_{TMC})
          ratio->SetPointError(ratio->GetN() - 1, 0, 0, r_err, r_err);
          ratio_simple->SetPointError(ratio_simple->GetN() - 1, 0, 0, r_err_simple, r_err_simple);
          break;
        }
      }
    }

    ratio->Draw("AP");
    ratio_simple->Draw("P SAME");
    func->Draw("L SAME");
    TLegend* leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(ratio, ratio->GetTitle(), "p");
    leg->AddEntry(ratio_simple, ratio_simple->GetTitle(), "p");
    leg->AddEntry(func, "ToyMC th. noise", "l");
    leg->Draw();
  }
  c->Update();
  c->Draw();
}

#endif
