// This file contains the implementation of a 2D mass fit analysis using RooFit and ROOT.
// The analysis consists of three parts:
// 1. Fitting the D0 mass (1D fit) to extract signal and background parameters.
// 2. Fitting the ΔM (DeltaM) mass (1D fit) to extract its parameters.
// 3. Performing a 2D fit by combining the two distributions, with various BDT cut values.
// 
// The code uses a TTree to read events, creates RooDataSets from histograms, defines signal and background models,
// performs the fits, and prints the resulting yields. The function `Automated_2D_Classifier` can be called
// from a ROOT macro to run the entire analysis.

// Include necessary RooFit, ROOT, and standard library headers.
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooBifurGauss.h"
#include "RooPlot.h"
#include "RooJohnson.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TLine.h"
#include "RooHist.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

// Structure to store parameters from the 1D D0 mass fit
struct D0FitParams {
    double mean;
    double lambda;
    double gamma;
    double delta;
    double sigma_L;
    double sigma_R;
    double f_bifur;
    double a1;
    double a2;
    double nsig;
    double nbkg;
};

// Structure to store parameters from the 1D ΔM mass fit
struct DeltaMFitParams {
    double mean;
    double lambda;
    double gamma;
    double delta;
    double alpha;
    double beta;
    double nsig;
    double nbkg;
};

// Structure to store the yields from the 2D mass fit
struct TwoDYields {
    double signal;
    double combinatorial;
    double random_slowpion;
    double misreco_D0;
};

//
// Helper function to plot the 2D mass fits (D0 and ΔM) and their pull distributions.
// The plots are produced on two separate canvases (one for D0 and one for ΔM) and can be saved
// using a suffix that identifies the cut values applied.
//
void Plot2DMass(RooRealVar &D0_refit, RooRealVar &deltam_refit,
                RooDataSet &data, RooAddPdf &total, const std::string &suffix) {
  // Create a frame for the D0 mass distribution
  RooPlot* d0_frame = D0_refit.frame(RooFit::Title("D0 Mass Fit"));
  data.plotOn(d0_frame, RooFit::Name("d0_data"));
  total.plotOn(d0_frame, RooFit::LineColor(kBlue), RooFit::LineWidth(8), RooFit::Name("total_D0_fit"));

  // Create a frame for the ΔM mass distribution
  RooPlot* deltam_frame = deltam_refit.frame(RooFit::Title("ΔM Mass Fit"));
  data.plotOn(deltam_frame, RooFit::Name("deltam_data"));
  total.plotOn(deltam_frame, RooFit::LineColor(kBlue), RooFit::LineWidth(8), RooFit::Name("total_deltam_fit"));

  // Create pull distribution frames for D0 and ΔM
  RooPlot* d0_pull_frame = D0_refit.frame(RooFit::Title("Pull Distribution D0"));
  RooHist* d0_pullHist = d0_frame->pullHist();
  d0_pull_frame->addPlotable(d0_pullHist, "P");

  RooPlot* deltam_pull_frame = deltam_refit.frame(RooFit::Title("Pull Distribution ΔM"));
  RooHist* deltam_pullHist = deltam_frame->pullHist();
  deltam_pull_frame->addPlotable(deltam_pullHist, "P");

  // Overlay individual signal and background components on the D0 frame
  total.plotOn(d0_frame, RooFit::Components("signal"), RooFit::LineColor(kRed),
               RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("D0_signal_fit"));
  total.plotOn(d0_frame, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2),
               RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::Name("D0_combinatorial_fit"));
  total.plotOn(d0_frame, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3),
               RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("D0_random_slowpion_fit"));
  total.plotOn(d0_frame, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7),
               RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::Name("D0_misreco_fit"));
  data.plotOn(d0_frame, RooFit::Name("d0_data"));

  // Overlay individual signal and background components on the ΔM frame
  total.plotOn(deltam_frame, RooFit::Components("signal"), RooFit::LineColor(kRed),
               RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("deltam_signal_fit"));
  total.plotOn(deltam_frame, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2),
               RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::Name("deltam_combinatorial_fit"));
  total.plotOn(deltam_frame, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3),
               RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("deltam_random_slowpion_fit"));
  total.plotOn(deltam_frame, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7),
               RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::Name("deltam_misreco_fit"));
  data.plotOn(deltam_frame, RooFit::Name("deltam_data"));
  total.paramOn(d0_frame);
  total.paramOn(deltam_frame);

  // Draw the D0 mass plot on a canvas with a pull distribution below
  double max_data_point = d0_frame->GetMaximum();
  d0_frame->SetMaximum(1.1 * max_data_point);
  d0_frame->SetMinimum(1);
  TCanvas* c1 = new TCanvas(("c1_" + suffix).c_str(), "D0 Fit Canvas", 900, 900);
  c1->SetTitle("");
  c1->Divide(1, 2, 0.0, 0.0);
  c1->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
  d0_frame->GetYaxis()->SetTitleFont(62);
  d0_frame->GetYaxis()->SetLabelFont(62);
  d0_frame->GetYaxis()->SetTitleSize(0.06);
  d0_frame->GetYaxis()->SetLabelSize(0.06);
  d0_frame->GetYaxis()->SetTitleOffset(1.4);
  d0_frame->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
  d0_frame->Draw();
  TLegend* legend1 = new TLegend(0.15, 0.7, 0.4, 0.9);
  legend1->SetTextFont(62);
  legend1->AddEntry(d0_frame->findObject("d0_data"), "Data", "p");
  legend1->AddEntry(d0_frame->findObject("total_D0_fit"), "Total", "l");
  legend1->AddEntry(d0_frame->findObject("D0_signal_fit"), "Signal", "l");
  legend1->AddEntry(d0_frame->findObject("D0_random_slowpion_fit"), "Random Slow Pion", "l");
  legend1->AddEntry(d0_frame->findObject("D0_combinatorial_fit"), "Combinatorial", "l");
  legend1->AddEntry(d0_frame->findObject("D0_misreco_fit"), "Misreconstructed D^{0}", "l");
  legend1->SetBorderSize(0);
  legend1->Draw();
  c1->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
  d0_pull_frame->GetYaxis()->SetTitle("Pull");
  d0_pull_frame->GetYaxis()->SetTitleSize(0.13);
  d0_pull_frame->GetYaxis()->SetTitleOffset(0.25);
  d0_pull_frame->GetYaxis()->SetLabelSize(0.10);
  d0_pull_frame->GetXaxis()->SetLabelSize(0.13);
  d0_pull_frame->GetXaxis()->SetTitleSize(0.13);
  d0_pull_frame->GetXaxis()->SetTitleFont(62);
  d0_pull_frame->GetYaxis()->SetTitleFont(62);
  d0_pull_frame->GetYaxis()->SetLabelFont(62);
  d0_pull_frame->GetXaxis()->SetLabelFont(62);
  d0_pull_frame->GetYaxis()->SetNdivisions(6);
  d0_pull_frame->GetYaxis()->SetRangeUser(-5, 5);
  d0_pull_frame->GetXaxis()->SetNdivisions(6);
  d0_pull_frame->Draw();
  TLine* zeroLine = new TLine(1820., 0.0, 1910., 0.0);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(9);
  zeroLine->SetLineWidth(2);
  zeroLine->Draw();
  c1->Update();

  // The ΔM canvas plotting block is currently commented out.
  // Uncomment and adjust if separate plots for ΔM are required.
  /*
  double max_data_pointdm = deltam_frame->GetMaximum();
  deltam_frame->SetMaximum(1.1 * max_data_pointdm);
  deltam_frame->SetMinimum(1);
  TCanvas* c2 = new TCanvas(("c2_" + suffix).c_str(), "ΔM Fit Canvas", 900, 900);
  c2->SetTitle("");
  c2->Divide(1, 2, 0.0, 0.0);
  c2->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
  deltam_frame->GetYaxis()->SetTitleFont(62);
  deltam_frame->GetYaxis()->SetLabelFont(62);
  deltam_frame->GetYaxis()->SetTitleSize(0.06);
  deltam_frame->GetYaxis()->SetLabelSize(0.06);
  deltam_frame->GetYaxis()->SetTitleOffset(1.3);
  deltam_frame->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
  deltam_frame->Draw();
  TLegend* legend2 = new TLegend(0.15, 0.7, 0.4, 0.9);
  legend2->SetTextFont(62);
  legend2->AddEntry(deltam_frame->findObject("deltam_data"), "Data", "p");
  legend2->AddEntry(deltam_frame->findObject("total_deltam_fit"), "Total", "l");
  legend2->AddEntry(deltam_frame->findObject("deltam_signal_fit"), "Signal", "l");
  legend2->AddEntry(deltam_frame->findObject("deltam_random_slowpion_fit"), "Random Slow Pion", "l");
  legend2->AddEntry(deltam_frame->findObject("deltam_combinatorial_fit"), "Combinatorial", "l");
  legend2->AddEntry(deltam_frame->findObject("deltam_misreco_fit"), "Mis-reconstructed D^{0}", "l");
  legend2->SetBorderSize(0);
  legend2->Draw();
  c2->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
  deltam_pull_frame->GetYaxis()->SetTitle("Pull");
  deltam_pull_frame->GetYaxis()->SetTitleSize(0.13);
  deltam_pull_frame->GetYaxis()->SetTitleOffset(0.25);
  deltam_pull_frame->GetYaxis()->SetLabelSize(0.10);
  deltam_pull_frame->GetXaxis()->SetLabelSize(0.13);
  deltam_pull_frame->GetXaxis()->SetTitleSize(0.13);
  deltam_pull_frame->GetXaxis()->SetTitleFont(62);
  deltam_pull_frame->GetYaxis()->SetTitleFont(62);
  deltam_pull_frame->GetYaxis()->SetLabelFont(62);
  deltam_pull_frame->GetXaxis()->SetLabelFont(62);
  deltam_pull_frame->GetYaxis()->SetRangeUser(-5, 5);
  deltam_pull_frame->GetYaxis()->SetNdivisions(6);
  deltam_pull_frame->GetXaxis()->SetNdivisions(6);
  deltam_pull_frame->Draw();
  TLine* zeroLine2 = new TLine(142.0, 0.0, 152.0, 0.0);
  zeroLine2->SetLineColor(kBlack);
  zeroLine2->SetLineStyle(9);
  zeroLine2->SetLineWidth(2);
  zeroLine2->Draw();
  c2->Update();
  */
}

//
// Function to perform the 1D D0 mass fit.
// This function opens the ROOT file, applies cuts based on BDTG responses,
// fills a histogram with the selected D0 mass values, converts the histogram to a RooDataSet,
// defines the signal (Johnson + bifurcated Gaussian) and background (Chebychev) models,
// performs the fit, prints the fit results, and returns the fitted parameters.
//
D0FitParams FitD0Mass(const char* file_path, float comb_cut, float rsp_cut) {
    // Open the ROOT file and retrieve the TTree
    TFile* file = TFile::Open(file_path);
    TTree* tree = (TTree*)file->Get("outtree");

    // Variables read from the tree branches
    float BDTG_response_comb;
    float BDTG_response_rsp;
    float D0_refit_mass;
    tree->SetBranchAddress("BDTG_response_comb", &BDTG_response_comb);
    tree->SetBranchAddress("BDTG_response_rsp", &BDTG_response_rsp);
    tree->SetBranchAddress("Dst_ReFit_D0_M_best", &D0_refit_mass);

    // Create a histogram to store the D0 mass distribution (binning may be adjusted as needed)
    TH1F* h_mass = new TH1F("h_mass", "D0 Refit Mass;Mass (MeV/c^{2});Events", 100, 1820, 1910);
    Long64_t nentries = tree->GetEntries();
    int selectedEvents = 0;
    // Loop over events and fill the histogram if the event passes the BDTG cuts
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (BDTG_response_comb > comb_cut && BDTG_response_rsp > rsp_cut) {
            h_mass->Fill(D0_refit_mass);
            selectedEvents++;
        }
    }
    std::cout << "D0Fit: Total entries = " << nentries << ", Selected = " << selectedEvents << std::endl;

    // Define the RooRealVar for D0 mass and build a RooDataSet from the histogram
    RooRealVar D0_refit("Dst_ReFit_D0_M_best", "m(D^{0}) [MeV/c^{2}]", 1820., 1910.);
    RooDataSet data("data", "dataset", RooArgSet(D0_refit));
    for (int bin = 1; bin <= h_mass->GetNbinsX(); bin++) {
        double content = h_mass->GetBinContent(bin);
        double center = h_mass->GetBinCenter(bin);
        for (int j = 0; j < content; j++) {
            D0_refit.setVal(center);
            data.add(RooArgSet(D0_refit));
        }
    }
    std::cout << "D0Fit: Entries in RooDataSet = " << data.numEntries() << std::endl;

    // Define the signal model: a combination of a Johnson function and a bifurcated Gaussian
    RooRealVar mean("mean", "Signal mean", 1865., 1863., 1867);
    RooRealVar lambda("lambda", "deltam_lambda", 9, 0, 15);
    RooRealVar gamma("gamma", "deltam_gamma", 0.31, -2, 2);
    RooRealVar delta("delta", "deltam_delta", 1.38, 0, 5);
    RooJohnson d0_johnson("d0_johnson", "Johnson signal", D0_refit, mean, lambda, gamma, delta);

    RooRealVar sigma_L("sigma_L", "Left width", 6, 0, 100);
    RooRealVar sigma_R("sigma_R", "Right width", 7, 0, 20);
    RooBifurGauss bifurGauss("bifurGauss", "Bifurcated Gaussian", D0_refit, mean, sigma_L, sigma_R);

    RooRealVar f_bifur("f_bifur", "Fraction", 0.4, 0.0, 1.0);
    RooAddPdf signal("signal", "Signal Model", RooArgList(d0_johnson, bifurGauss), RooArgList(f_bifur));

    // Define the background model using a Chebychev polynomial
    RooRealVar a1("a1", "Coefficient a1", -0.16, -1.0, 2.0);
    RooRealVar a2("a2", "Coefficient a2", 0.11, -1.0, 2.0);
    RooChebychev background("background", "Chebychev Bkg", D0_refit, RooArgList(a1, a2));

    // Define the total PDF as the sum of signal and background components
    RooRealVar nsig("nsig", "Number of signal events", 4090000, 0, 1e8);
    RooRealVar nbkg("nbkg", "Number of background events", 35000000, 0, 1e9);
    RooAddPdf total("total", "Signal+Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // Perform the fit and print the results
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print("v");

    // Store the fitted parameters in the structure to return
    D0FitParams params;
    params.mean   = mean.getVal();
    params.lambda = lambda.getVal();
    params.gamma  = gamma.getVal();
    params.delta  = delta.getVal();
    params.sigma_L = sigma_L.getVal();
    params.sigma_R = sigma_R.getVal();
    params.f_bifur = f_bifur.getVal();
    params.a1     = a1.getVal();
    params.a2     = a2.getVal();
    params.nsig   = nsig.getVal();
    params.nbkg   = nbkg.getVal();

    file->Close();
    return params;
}

//
// Function to perform the 1D ΔM mass fit.
// This function follows a similar procedure as FitD0Mass: it opens the ROOT file,
// applies the selection cuts, creates a histogram and RooDataSet for ΔM, defines the signal
// and background PDFs (here the background is modeled with a polynomial via RooGenericPdf),
// performs the fit, prints the results, and returns the fitted parameters.
//
DeltaMFitParams FitDeltaMMass(const char* file_path, float comb_cut, float rsp_cut) {
    // Open the ROOT file and retrieve the TTree
    TFile* file = TFile::Open(file_path);
    TTree* tree = (TTree*)file->Get("outtree");

    // Variables read from the tree branches
    float BDTG_response_comb;
    float BDTG_response_rsp;
    float deltam_refit_mass;
    tree->SetBranchAddress("BDTG_response_comb", &BDTG_response_comb);
    tree->SetBranchAddress("BDTG_response_rsp", &BDTG_response_rsp);
    tree->SetBranchAddress("deltam_ReFit", &deltam_refit_mass);

    // Create a histogram for the ΔM mass distribution (range and binning can be adjusted)
    TH1F* h_mass = new TH1F("h_mass", "DeltaM Refit Mass;Mass (MeV/c^{2});Events", 100, 142., 152.);
    Long64_t nentries = tree->GetEntries();
    int selectedEvents = 0;
    // Loop over events and fill the histogram if the event passes the BDTG cuts
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (BDTG_response_comb > comb_cut && BDTG_response_rsp > rsp_cut) {
            h_mass->Fill(deltam_refit_mass);
            selectedEvents++;
        }
    }
    std::cout << "DeltaMFit: Total entries = " << nentries << ", Selected = " << selectedEvents << std::endl;

    // Define the RooRealVar for ΔM and build a RooDataSet from the histogram
    RooRealVar deltam_refit("deltam_ReFit", "#DeltaM [MeV/c^{2}]", 142., 152.);
    RooDataSet data("data", "dataset", RooArgSet(deltam_refit));
    for (int bin = 1; bin <= h_mass->GetNbinsX(); bin++) {
        double content = h_mass->GetBinContent(bin);
        double center = h_mass->GetBinCenter(bin);
        for (int j = 0; j < content; j++) {
            deltam_refit.setVal(center);
            data.add(RooArgSet(deltam_refit));
        }
    }
    std::cout << "DeltaMFit: Entries in RooDataSet = " << data.numEntries() << std::endl;

    // Define the signal model using a Johnson PDF
    RooRealVar mean("mean", "DeltaM Mean", 145, 144, 146);
    RooRealVar lambda("lambda", "DeltaM Lambda", 0.5, 0, 1);
    RooRealVar gamma("gamma", "DeltaM Gamma", -0.7, -1, 1);
    RooRealVar delta("delta", "DeltaM Delta", 1, 0, 10);
    RooJohnson signal("signal", "DeltaM Johnson", deltam_refit, mean, lambda, gamma, delta);

    // Define the background model as a polynomial using RooGenericPdf
    RooRealVar alpha("alpha", "alpha", 1.0, 0., 10.);
    RooRealVar beta("beta", "beta", 0.1, -10., 10.);
    RooGenericPdf background("background", "Poly background",
                             "pow((deltam_ReFit - 139.5), alpha) * exp(-1*beta*(deltam_ReFit - 139.5))",
                             RooArgSet(deltam_refit, alpha, beta));

    // Define the total PDF as the sum of signal and background components
    RooRealVar nsig("nsig", "Number of signal events", 13023, 0, 1e7);
    RooRealVar nbkg("nbkg", "Number of background events", 2074774, 0, 1e7);
    RooAddPdf total("total", "Signal+Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // Perform the fit and print the results
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print("v");

    // Store the fitted parameters in the structure to return
    DeltaMFitParams params;
    params.mean   = mean.getVal();
    params.lambda = lambda.getVal();
    params.gamma  = gamma.getVal();
    params.delta  = delta.getVal();
    params.alpha  = alpha.getVal();
    params.beta   = beta.getVal();
    params.nsig   = nsig.getVal();
    params.nbkg   = nbkg.getVal();

    file->Close();
    return params;
}

//
// Function to perform the 2D mass fit.
// It reads the ROOT file and applies 2D cuts on BDTG responses,
// builds a 2D RooDataSet using both D0 and ΔM mass variables,
// constructs signal PDFs for each mass dimension using parameters fixed from the 1D fits,
// builds background PDFs, and finally constructs the total 2D PDF as a sum of several components.
// After performing the fit, it returns the yields of the different components.
//
TwoDYields Fit2DMass(const char* file_path, float comb_cut, float rsp_cut,
                      const D0FitParams &d0Params,
                      const DeltaMFitParams &deltaMParams,
                      const std::string &suffix) {
  // Open the ROOT file and retrieve the TTree
  TFile* file = TFile::Open(file_path);
  TTree* tree = (TTree*)file->Get("outtree");

  // Define the RooRealVars for D0 and ΔM masses
  RooRealVar D0_refit("Dst_ReFit_D0_M_best", "m(D^{0}) [MeV/c^{2}]", 1820, 1910);
  RooRealVar deltam_refit("deltam_ReFit", "#DeltaM [MeV/c^{2}]", 142, 152);

  // Build the 2D RooDataSet using the cuts on both BDTG responses
  RooDataSet data("data", "dataset", RooArgSet(D0_refit, deltam_refit));
  Long64_t nentries = tree->GetEntries();
  float D0_refit_mass, deltam_refit_mass;
  float BDTG_response_comb, BDTG_response_rsp;
  tree->SetBranchAddress("BDTG_response_comb", &BDTG_response_comb);
  tree->SetBranchAddress("BDTG_response_rsp",  &BDTG_response_rsp);
  tree->SetBranchAddress("Dst_ReFit_D0_M_best", &D0_refit_mass);
  tree->SetBranchAddress("deltam_ReFit", &deltam_refit_mass);
  
  // Loop over events and add to the dataset if the event passes both cuts
  for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    if (BDTG_response_comb > comb_cut && BDTG_response_rsp > rsp_cut) {
      D0_refit.setVal(D0_refit_mass);
      deltam_refit.setVal(deltam_refit_mass);
      data.add(RooArgSet(D0_refit, deltam_refit));
    }
  }
  std::cout << "2DFit: For suffix " << suffix << " => Entries in RooDataSet = " << data.numEntries() << std::endl;

  // Build the D0 signal PDF using fixed parameters from the 1D D0 mass fit
  RooRealVar d0_mean("d0_mean", "D0 Mean", d0Params.mean);
  RooRealVar d0_lambda("d0_lambda", "D0 Lambda", d0Params.lambda);
  RooRealVar d0_gamma("d0_gamma", "D0 Gamma", d0Params.gamma);
  RooRealVar d0_delta("d0_delta", "D0 Delta", d0Params.delta);
  RooJohnson d0_johnson("d0_johnson", "D0 Johnson", D0_refit, d0_mean, d0_lambda, d0_gamma, d0_delta);
  RooRealVar d0_sigma_L("d0_sigma_L", "Left sigma", d0Params.sigma_L);
  RooRealVar d0_sigma_R("d0_sigma_R", "Right sigma", d0Params.sigma_R);
  RooBifurGauss bifurGauss("bifurGauss", "Bifur Gauss", D0_refit, d0_mean, d0_sigma_L, d0_sigma_R);
  RooRealVar f_bifur("f_bifur", "Fraction", d0Params.f_bifur);
  RooAddPdf d0_signal_pdf("d0_signal_pdf", "D0 Signal PDF",
                          RooArgList(d0_johnson, bifurGauss), RooArgList(f_bifur));

  // Build the ΔM signal PDF using fixed parameters from the 1D ΔM mass fit
  RooRealVar deltam_mean("deltam_mean", "ΔM Mean", deltaMParams.mean);
  RooRealVar deltam_lambda("deltam_lambda", "ΔM Lambda", deltaMParams.lambda);
  RooRealVar deltam_gamma("deltam_gamma", "ΔM Gamma", deltaMParams.gamma);
  RooRealVar deltam_delta("deltam_delta", "ΔM Delta", deltaMParams.delta);
  RooJohnson deltam_signal_pdf("deltam_signal_pdf", "ΔM Signal", deltam_refit,
                               deltam_mean, deltam_lambda, deltam_gamma, deltam_delta);

  // Build background PDFs for each dimension using the fixed parameters
  RooRealVar bkg_alpha("bkg_alpha", "alpha", deltaMParams.alpha);
  RooRealVar bkg_beta("bkg_beta", "beta", deltaMParams.beta);
  RooGenericPdf deltam_background_pdf("deltam_background_pdf", "Poly bkg",
      "pow((deltam_ReFit-139.5),bkg_alpha)*exp(-1*bkg_beta*(deltam_ReFit-139.5))",
      RooArgSet(deltam_refit, bkg_alpha, bkg_beta));
  RooRealVar a1("a1", "a1", d0Params.a1);
  RooRealVar a2("a2", "a2", d0Params.a2);
  RooChebychev d0_background_pdf("d0_background_pdf", "Chebychev bkg", D0_refit, RooArgList(a1, a2));

  // Construct the 2D PDFs as products of the 1D PDFs
  RooProdPdf signal("signal", "Signal 2D", RooArgList(d0_signal_pdf, deltam_signal_pdf));
  RooProdPdf combinatorial("combinatorial", "Combinatorial bkg", RooArgList(d0_background_pdf, deltam_background_pdf));
  RooProdPdf random_slowpion("random_slowpion", "Random slow pion",
                              RooArgList(d0_signal_pdf, deltam_background_pdf));
  RooProdPdf misreco_D0("misreco_D0", "Misreconstructed D0",
                        RooArgList(d0_background_pdf, deltam_signal_pdf));

  // Define the yields for each component of the 2D fit
  RooRealVar Signal_yield("Signal_yield", "Signal yield", 10000, 0, 1e7);
  RooRealVar Combinatorial_yield("Combinatorial_yield", "Comb. yield", 10000000, 0, 1e8);
  RooRealVar Random_slowpion_yield("Random_slowpion_yield", "Random slow yield", 10000000, 0, 1e8);
  RooRealVar Misreco_yield("Misreco_yield", "Misreco yield", 100000, 200, 1e7);

  // Define the total 2D PDF as a sum of the four components
  RooAddPdf total("total", "Total 2D PDF",
                  RooArgList(signal, combinatorial, random_slowpion, misreco_D0),
                  RooArgList(Signal_yield, Combinatorial_yield, Random_slowpion_yield, Misreco_yield));

  // Perform the 2D fit and print the results
  RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
  fit_result->Print("v");

  // Retrieve the yields from the fit results
  TwoDYields yields;
  yields.signal = Signal_yield.getVal();
  yields.combinatorial = Combinatorial_yield.getVal();
  yields.random_slowpion = Random_slowpion_yield.getVal();
  yields.misreco_D0 = Misreco_yield.getVal();

  // Optionally produce the plots (currently commented out)
  // Plot2DMass(D0_refit, deltam_refit, data, total, suffix);

  file->Close();
  return yields;
}

//
// Main function: loops over arrays of BDT cut values, performs 1D fits to extract parameters,
// runs the 2D mass fit for each combination of cuts, and prints a summary of the yields and significance.
//
int main() {
    // Define arrays of BDT cut values for the combinatorial (comb) and random slow pion (rsp) selections
    std::vector<float> BDTcuts_comb;
    std::vector<float> BDTcuts_rsp;

    for (float cut = -1; cut <= 0.4; cut += 0.05) {
        BDTcuts_comb.push_back(cut);
    }
    for (float cut = -1; cut <= -0.8; cut += 0.05) {
        BDTcuts_rsp.push_back(cut);
    }

    // Vectors to store the yields and corresponding cut values for later summary output
    std::vector<float> combCuts, rspCuts;
    std::vector<double> signalYields, combYields, randSlowYields, misrecoYields;

    // Use a generic file name (change as necessary)
    const char* filePath1 = "input.root";
    
    // Loop over each combination of BDT cut values
    for (size_t i = 0; i < BDTcuts_comb.size(); i++) {
        for (size_t j = 0; j < BDTcuts_rsp.size(); j++) {
            float comb_cut = BDTcuts_comb[i];
            float rsp_cut  = BDTcuts_rsp[j];
            
            // Run the 1D fits to get fixed parameters for the 2D fit
            D0FitParams d0Params = FitD0Mass(filePath1, comb_cut, rsp_cut);
            DeltaMFitParams deltaMParams = FitDeltaMMass(filePath1, comb_cut, rsp_cut);

            // Create a suffix string based on the current cut values for plotting/identification purposes
            std::ostringstream oss;
            oss << "comb_" << comb_cut << "_rsp_" << rsp_cut;
            std::string suffix = oss.str();

            std::cout << "=============================================" << std::endl;
            std::cout << "Running 2D fit for:" << std::endl;
            std::cout << "  BDTG_response_comb > " << comb_cut << std::endl;
            std::cout << "  BDTG_response_rsp  > " << rsp_cut << std::endl;

            // Perform the 2D mass fit using the fixed parameters from the 1D fits
            TwoDYields yields = Fit2DMass(filePath1, comb_cut, rsp_cut, d0Params, deltaMParams, suffix);

            // Store the cut values and yields
            combCuts.push_back(comb_cut);
            rspCuts.push_back(rsp_cut);
            signalYields.push_back(yields.signal);
            combYields.push_back(yields.combinatorial);
            randSlowYields.push_back(yields.random_slowpion);
            misrecoYields.push_back(yields.misreco_D0);
        }
    }

    // Print a summary table of the yields and calculated significance for each cut combination
    std::cout << "\nSummary of all 2D fit results:" << std::endl;
    std::cout << "Comb Cut, RSP Cut, Signal Yield, Combinatorial Yield, Random Slow Pion Yield, Misreco D0 Yield, Significance" << std::endl;
    for (size_t i = 0; i < signalYields.size(); i++) {
        double significance = signalYields[i] / sqrt(signalYields[i] + combYields[i] + randSlowYields[i] + misrecoYields[i]);
        std::cout << combCuts[i] << ", "
                  << rspCuts[i] << ", "
                  << signalYields[i] << ", "
                  << combYields[i] << ", "
                  << randSlowYields[i] << ", "
                  << misrecoYields[i] << ", "
                  << significance << std::endl;
    }

    return 0;
}

// This function can be invoked from a ROOT macro to run the complete analysis.
void Automated_2D_Classifier() {
    main();
}
