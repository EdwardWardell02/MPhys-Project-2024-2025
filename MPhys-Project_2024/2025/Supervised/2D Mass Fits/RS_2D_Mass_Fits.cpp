// This script performs a two-dimensional mass fit for RS (Right-Sign) events
// by fitting the Δm and D0 refitted mass (Dst_ReFit_D0_M_best) distributions.
// It creates separate RooDataSets for the two mass variables, defines fixed signal
// and background PDFs for both variables, and then combines them into a total PDF
// for a 2D fit. The final plots include the mass fits and their pull distributions.
// To run the analysis, call the function twoD_Mass_Fit_RS_DeltaM_with_Dst_ReFit_D0_M_best()
// from a ROOT macro or interactive ROOT session.

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooBifurGauss.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TLatex.h"

// Use the RooFit namespace for cleaner code
using namespace RooFit;

void twoD_Mass_Fit_RS_DeltaM_with_Dst_ReFit_D0_M_best() {
    // Replace the file path with a generic one for portability.
    // Uncomment the desired file path if you have multiple options.
    // const char* file_path = "input_log.root";
    const char* file_path = "input_sw_check.root";
    
    // Open the ROOT file and get the TTree named "DecayTree"
    TFile* root_file = TFile::Open(file_path);
    TTree* tree = (TTree*)root_file->Get("DecayTree");

    // ------------------------------
    // Define Variables and Create Datasets
    // ------------------------------
    // Create RooRealVars for the Δm and D0 refitted mass with their respective ranges.
    RooRealVar deltam_refit("deltam_ReFit", "#Delta m [MeV/c^{2}]", 142.0, 152.0);
    RooRealVar D0_refit("Dst_ReFit_D0_M_best", "m(D^{0}) [MeV/c^{2}]", 1820, 1910);
    
    // Create a combined dataset containing both mass variables.
    RooDataSet data("data", "Data", tree, RooArgSet(D0_refit, deltam_refit));
    // Also create individual datasets for each variable (useful for 1D projections).
    RooDataSet d0_data("d0_data", "Data", tree, RooArgSet(D0_refit));
    RooDataSet deltam_data("deltam_data", "Data", tree, RooArgSet(deltam_refit));

    //---------------------------------
    // Define Fixed PDFs for Δm and D0
    //---------------------------------
    // --- Δm Fixed PDFs ---
    // Set fixed parameters for the Δm signal PDF.
    RooRealVar deltam_mean("deltam_mean", "DeltaM Mean", 145.44907);
    RooRealVar deltam_lambda("deltam_lambda", "DeltaM Lambda", 0.38952);
    RooRealVar deltam_gamma("deltam_gamma", "DeltaM Gamma", 0.1130);
    RooRealVar deltam_delta("deltam_delta", "DeltaM Delta", 1.4235);
    // Define additional parameters for the Δm PDF components.
    RooRealVar deltam_sigma("deltam_sigma", "DeltaM Sigma", 1.097);
    RooRealVar deltam_sigma_R("deltam_sigma_R", "DeltaM Right Sigma", 0.459, 0., 500);
    RooRealVar deltam_sigma_L("deltam_sigma_L", "DeltaM left Sigma", 1.241, 0., 500);
    RooRealVar deltam_f_bifur("deltam_f_bifur", "DeltaM Fraction of bifur", 0.7815);
    RooRealVar deltam_f_gauss("deltam_f_gauss", "DeltaM Fraction of Gaussian", 0.1218);
    
    // Construct Δm signal PDF components:
    RooGaussian deltam_Gauss("deltam_Gauss", "DeltaM Gaussian", deltam_refit, deltam_mean, deltam_sigma);
    RooJohnson deltam_Johnson("deltam_Johnson", "DeltaM Johnson", deltam_refit, deltam_mean, deltam_lambda, deltam_gamma, deltam_delta);
    RooBifurGauss deltam_bifurGauss("deltam_bifurGauss", "DeltaM Bifurcated Gaussian", deltam_refit, deltam_mean, deltam_sigma_L, deltam_sigma_R);
    // Combine Δm signal components into a single PDF.
    RooAddPdf deltam_signal_pdf("deltam_signal_pdf", "DeltaM Signal PDF", 
                                RooArgList(deltam_bifurGauss, deltam_Johnson, deltam_Gauss),
                                RooArgList(deltam_f_gauss, deltam_f_bifur));

    // --- D0 Fixed PDFs ---
    // Set fixed parameters for the D0 signal PDF.
    RooRealVar d0_mean("d0_mean", "D0 Mean", 1864.0927);
    RooRealVar d0_lambda("d0_lambda", "D0 lambda", 17.00009);
    RooRealVar d0_gamma("d0_gamma", "D0 gamma", -0.209282);
    RooRealVar d0_delta("d0_delta", "D0 delta", 2.35568);
    RooRealVar d0_sigma_R("d0_sigma_R", "D0 Right Sigma", 10.87);
    RooRealVar d0_sigma_L("d0_sigma_L", "D0 Left Sigma", 15.01);
    RooRealVar d0_f_DG("d0_f_DG", "Fraction of Double Gauss", 0.8560);
    // Build the D0 signal PDF from a Johnson and a bifurcated Gaussian.
    RooJohnson d0_Johnson("d0_Johnso", "D0 Johnson", D0_refit, d0_mean, d0_lambda, d0_gamma, d0_delta);
    RooBifurGauss d0_refit_bifurGauss("d0_refit_bifurGauss", "Dst_ReFit_D0_M_best Bifurcated Gaussian", D0_refit, d0_mean, d0_sigma_L, d0_sigma_R);
    RooAddPdf d0_refit_signal_pdf("d0_refit_signal_pdf", "Dst_ReFit_D0_M_best Signal PDF", 
                                  RooArgList(d0_Johnson, d0_refit_bifurGauss),
                                  RooArgList(d0_f_DG));

    // Define the D0 background using an exponential with a negative slope.
    RooRealVar d0_lambda_bkg("#lambda_{bkg}^{D^{0}}", "lambda", -0.0010990, -1, 0.0);
    RooExponential d0_refit_background_pdf("d0_refit_background_pdf", "Exponential Background", D0_refit, d0_lambda_bkg);

    // -------- Δm Background PDF --------
    // Use a generic polynomial background model for Δm.
    RooRealVar alpha("alpha", "alpha", 0.5832);
    RooRealVar beta("beta", "beta", 0.06237);
    RooGenericPdf deltam_background_pdf(
        "background", "Polynomial background", 
        "pow((deltam_ReFit - 139.5), alpha) * exp(-1 * beta * (deltam_ReFit - 139.5))",
        RooArgSet(deltam_refit,  alpha, beta));

    //-------------------------------
    // Combine PDFs into 2D Components
    //-------------------------------
    // Combine D0 and Δm fixed signal PDFs into a product PDF.
    RooProdPdf signal("signal", "signal", RooArgList(d0_refit_signal_pdf, deltam_signal_pdf));
    // Similarly, combine the background PDFs.
    RooProdPdf combinatorial("combinatorial", "Combinatorial background", RooArgList(d0_refit_background_pdf, deltam_background_pdf));
    RooProdPdf random_slowpion("random_slowpion", "true D0 with incorrect slow pion", RooArgList(d0_refit_signal_pdf, deltam_background_pdf));
    RooProdPdf misreco_D0("misreco_D0", "incorrect D0 with true slow pion", RooArgList(d0_refit_background_pdf, deltam_signal_pdf));
    
    //-----------------------
    // Define Yields for Each Component
    //-----------------------
    RooRealVar Signal_yield("Signal_yield", "Number of signal events for DeltaM", 4000000, 0, 1000000000);
    RooRealVar Combinatorial_yield("Combinatorial_yield", "Number of background events for DeltaM", 1000000, 0, 1000000000);
    RooRealVar Random_slowpion_yield("Random_slowpion_yield", "Number of signal events for DeltaM", 1000000, 0, 1000000000);
    RooRealVar Misreco_yield("Misreco_yield", "Number of background events for DeltaM", 100000, 0, 1000000000);

    //--------------------------------------------------
    // Combine All Components into the Total PDF for the 2D Fit
    //--------------------------------------------------
    RooAddPdf total("total", "Signal + Background",
        RooArgList(signal, combinatorial, random_slowpion, misreco_D0),
        RooArgList(Signal_yield, Combinatorial_yield, Random_slowpion_yield, Misreco_yield));

    // ------------------------------
    // Perform the 2D Fit
    // ------------------------------
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print();

    // ------------------------------
    // Create Projections for D0 and Δm and Plot Results
    // ------------------------------
    // D0 Mass Projection
    RooPlot* d0_frame = D0_refit.frame(RooFit::Title("Mass Fit"));
    data.plotOn(d0_frame);
    total.plotOn(d0_frame, RooFit::LineColor(kBlue), RooFit::LineWidth(8), RooFit::Name("total_D0_fit"));

    // Δm Mass Projection
    RooPlot* deltam_frame = deltam_refit.frame(RooFit::Title("Mass Fit"));
    data.plotOn(deltam_frame);
    total.plotOn(deltam_frame, RooFit::LineColor(kBlue), RooFit::LineWidth(8), RooFit::Name("total_deltam_fit"));

    // ------------------------------
    // Create Pull Distributions
    // ------------------------------
    RooPlot* d0_pull_frame = D0_refit.frame(RooFit::Title("Pull Distribution D0"));
    RooHist* d0_pullHist = d0_frame->pullHist();
    d0_pull_frame->addPlotable(d0_pullHist, "P");
    
    RooPlot* deltam_pull_frame = deltam_refit.frame(RooFit::Title("Pull Distribution DeltaM"));
    RooHist* deltam_pullHist = deltam_frame->pullHist();
    deltam_pull_frame->addPlotable(deltam_pullHist, "P");

    // Overlay individual components on the D0 projection.
    total.plotOn(d0_frame, RooFit::Components("signal"), RooFit::LineColor(kRed),
        RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("D0_signal_fit"));
    total.plotOn(d0_frame, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen + 2),
        RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::Name("D0_combinatorial_fit"));
    total.plotOn(d0_frame, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3),
        RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("D0_random_slowpion_fit"));
    total.plotOn(d0_frame, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7),
        RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::Name("D0_misreco_fit"));
    data.plotOn(d0_frame, RooFit::Name("d0_data"));

    // Overlay individual components on the Δm projection.
    total.plotOn(deltam_frame, RooFit::Components("signal"), RooFit::LineColor(kRed),
        RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("deltam_signal_fit"));
    total.plotOn(deltam_frame, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen + 2),
        RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::Name("deltam_combinatorial_fit"));
    total.plotOn(deltam_frame, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3),
        RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("deltam_random_slowpion_fit"));
    total.plotOn(deltam_frame, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7),
        RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::Name("deltam_misreco_fit"));
    data.plotOn(deltam_frame, RooFit::Name("deltam_data"));
    total.paramOn(d0_frame);
    total.paramOn(deltam_frame);

    //-----------------
    // D0 Mass Plotting
    //-----------------
    double max_data_point = d0_frame->GetMaximum();
    d0_frame->SetMaximum(1.2 * max_data_point);
    d0_frame->SetMinimum(1);
    TCanvas* c1 = new TCanvas("c1", "Fit Canvas", 900, 900);
    c1->SetTitle("");  // Disable canvas title
    c1->Divide(1, 2, 0.0, 0.0);
    c1->cd(1);
    c1->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
    d0_frame->GetYaxis()->SetTitleFont(62);
    d0_frame->GetYaxis()->SetLabelFont(62);
    d0_frame->GetYaxis()->SetTitleSize(0.06);
    d0_frame->GetYaxis()->SetLabelSize(0.06);
    d0_frame->GetYaxis()->SetTitleOffset(1);
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
    c1->cd(2);
    c1->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
    d0_pull_frame->GetYaxis()->SetTitle("Pull");
    d0_pull_frame->GetYaxis()->SetTitleSize(0.13); //Title size
    d0_pull_frame->GetYaxis()->SetTitleOffset(0.25); //Title Offset
    d0_pull_frame->GetYaxis()->SetLabelSize(0.10); //y-axis labels size
    d0_pull_frame->GetXaxis()->SetLabelSize(0.13); //x-axis labels size
    d0_pull_frame->GetXaxis()->SetTitleSize(0.13); //x-axis title
    d0_pull_frame->GetXaxis()->SetTitleFont(62);
    d0_pull_frame->GetYaxis()->SetTitleFont(62);
    d0_pull_frame->GetYaxis()->SetLabelFont(62);
    d0_pull_frame->GetXaxis()->SetLabelFont(62);
    d0_pull_frame->GetYaxis()->SetNdivisions(4);
    d0_pull_frame->GetYaxis()->SetRangeUser(-15, 15);
    d0_pull_frame->GetXaxis()->SetNdivisions(6);
    d0_pull_frame->Draw();
    TLine* zeroLine = new TLine(1820., 0.0, 1910., 0.0);
    zeroLine->SetLineColor(kBlack);
    zeroLine->SetLineStyle(9);
    zeroLine->SetLineWidth(2);
    zeroLine->Draw();
    c1->Draw();

    //----------------
    // DeltaM Plotting
    //----------------
    double max_data_pointdm = deltam_frame->GetMaximum();
    deltam_frame->SetMaximum(1.2 * max_data_pointdm);
    deltam_frame->SetMinimum(1);
    TCanvas* c2 = new TCanvas("c2", "Fit Canvas", 900, 900);
    c2->SetTitle("");  // Disable canvas title
    c2->Divide(1, 2, 0.0, 0.0);
    c2->cd(1);
    c2->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
    deltam_frame->GetYaxis()->SetTitleFont(62);
    deltam_frame->GetYaxis()->SetLabelFont(62);
    deltam_frame->GetYaxis()->SetTitleSize(0.06);
    deltam_frame->GetYaxis()->SetLabelSize(0.06);
    deltam_frame->GetYaxis()->SetTitleOffset(1.1);
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
    c2->cd(2);
    c2->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
    deltam_pull_frame->GetYaxis()->SetTitle("Pull");
    deltam_pull_frame->GetYaxis()->SetTitleSize(0.13); //Title size
    deltam_pull_frame->GetYaxis()->SetTitleOffset(0.25); //Title Offset
    deltam_pull_frame->GetYaxis()->SetLabelSize(0.10); //y-axis labels size
    deltam_pull_frame->GetXaxis()->SetLabelSize(0.13); //x-axis labels size
    deltam_pull_frame->GetXaxis()->SetTitleSize(0.13); //x-axis title
    deltam_pull_frame->GetXaxis()->SetTitleFont(62);
    deltam_pull_frame->GetYaxis()->SetTitleFont(62);
    deltam_pull_frame->GetYaxis()->SetLabelFont(62);
    deltam_pull_frame->GetXaxis()->SetLabelFont(62);
    deltam_pull_frame->GetYaxis()->SetRangeUser(-31, 31);
    deltam_pull_frame->GetYaxis()->SetNdivisions(4);
    deltam_pull_frame->GetXaxis()->SetNdivisions(6);
    deltam_pull_frame->Draw();
    TLine* zeroLine2 = new TLine(142.0, 0.0, 152.0, 0.0);
    zeroLine2->SetLineColor(kBlack);
    zeroLine2->SetLineStyle(9);
    zeroLine2->SetLineWidth(2);
    zeroLine2->Draw();
    c2->Draw();
}
