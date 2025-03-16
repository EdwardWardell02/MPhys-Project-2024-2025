// This script performs a 1D mass fit for the RS D0 refitted mass variable
// (Dst_ReFit_D0_M_best) using a signal model composed of a Johnson PDF combined
// with a bifurcated Gaussian, and an exponential background model.
// The code reads data from a ROOT TTree, performs the fit, plots the results
// (including a pull distribution), and prints key fit statistics.
//
// To run the analysis, call the function 
// RS_Dst_ReFit_D0_M_best_Bifur_sig_Negative_Exponential_bkg_1D_fit()
// from a ROOT macro or an interactive ROOT session.

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

void RS_Dst_ReFit_D0_M_best_Bifur_sig_powerlaw_bkg_1D_fit() {
    // Replace the absolute file path with a generic placeholder for portability.
    const char* file_path = "input.root";
    
    // Open the input ROOT file and retrieve the TTree named "DecayTree".
    TFile* root_file = TFile::Open(file_path);
    TTree* tree = (TTree*)root_file->Get("DecayTree");

    // ------------------------------
    // Define the mass variable and dataset
    // ------------------------------
    // Create a RooRealVar for the D0 refitted mass in the range [1820, 1910] MeV/c^2.
    RooRealVar D0_refit("Dst_ReFit_D0_M_best", "m(D^{0}) [MeV/c^{2}]", 1820., 1910.);
    // Build a RooDataSet from the TTree using the defined mass variable.
    RooDataSet data("data", "dataset", tree, RooArgSet(D0_refit));

    // ------------------------------
    // Signal Model: Johnson + Bifurcated Gaussian
    // ------------------------------
    // Define the common signal mean with allowed variation.
    RooRealVar mean("#mu", "Signal mean", 1865., 1863., 1867.);
    // Define Johnson PDF parameters.
    RooRealVar lambda("#lambda", "deltam_lambda", 19, 17, 21);
    RooRealVar gamma("#gamma", "deltam_gamma", 0.001, -1, 1);
    RooRealVar delta("#delta", "deltam_delta", 2.7, 0, 5);
    // Create the Johnson PDF using the D0_refit variable.
    RooJohnson d0_johnson("deltam_johnson", "signal", D0_refit, mean, lambda, gamma, delta);

    // Define parameters for the bifurcated Gaussian.
    RooRealVar sigma_L("#sigma_{L}", "Left width of Bifurcated Gaussian", 18, 15, 20);
    RooRealVar sigma_R("#sigma_{R}", "Right width of Bifurcated Gaussian", 12, 10, 15);
    RooBifurGauss bifurGauss("bifurGauss", "Bifurcated Gaussian", D0_refit, mean, sigma_L, sigma_R);

    // Define the fraction that weights the bifurcated Gaussian contribution.
    RooRealVar f_bifur("f_{bifur}", "Fraction of Bifurcated Gaussian", 0.9, 0.0, 1.0);

    // Combine the Johnson PDF and the bifurcated Gaussian into a final signal model.
    RooAddPdf signal("signal", "Double Gaussian + Bifurcated Gaussian Signal", 
                     RooArgList(d0_johnson, bifurGauss), RooArgList(f_bifur));
    
    // ------------------------------
    // Background Model: Exponential Background
    // ------------------------------
    // Define the exponential background parameter with a negative slope.
    RooRealVar lambda_bkg("#lambda_{bkg}", "Exponential slope", -0.29, -1.0, 0.0);
    // Create the exponential background model using the D0_refit variable.
    RooExponential background("background", "Exponential Background", D0_refit, lambda_bkg);

    // ------------------------------
    // Total PDF: Signal + Background
    // ------------------------------
    // Define the number of signal and background events.
    RooRealVar nsig("nsig", "Number of signal events", 4090000, 0, 100000000);
    RooRealVar nbkg("nbkg", "Number of background events", 35000000, 0, 1000000000);
    // Combine the signal and background PDFs into the total model.
    RooAddPdf total("total", "Signal + Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // ------------------------------
    // Perform the Fit
    // ------------------------------
    // Fit the total PDF to the data, saving the fit result and using Minos for error estimation.
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print();

    // ------------------------------
    // Plot the Fit Results
    // ------------------------------
    // Create a frame for the mass variable and plot the data.
    RooPlot* frame = D0_refit.frame(RooFit::Title("Mass Fit"));
    data.plotOn(frame);
    total.plotOn(frame, RooFit::LineColor(kBlue), RooFit::LineWidth(10), RooFit::Name("total_fit"));

    // ------------------------------
    // Add Chi2/NDOF Text to the Plot
    // ------------------------------
    // Calculate the chi-square per degree of freedom.
    double chi2 = frame->chiSquare();
    int ndata = data.numEntries();
    int ndf = fit_result->floatParsFinal().getSize();
    int dof = ndata - ndf;
    // Display the chi-square per NDOF on the plot.
    TText* chi2_text = new TText(0.6, 0.8, Form("#chi^{2}/NDOF = %.2f/%d", chi2, dof));
    chi2_text->SetTextSize(0.04);
    chi2_text->SetTextFont(42);
    chi2_text->SetTextAlign(12);
    chi2_text->Draw();

    // ------------------------------
    // Create the Pull Plot
    // ------------------------------
    // Create a separate frame for the pull distribution.
    RooPlot* pull_frame = D0_refit.frame(RooFit::Title("Pull Distribution"));
    // Extract the pull histogram from the main frame.
    RooHist* pullHist = frame->pullHist();
    pull_frame->addPlotable(pullHist, "P");

    // Overlay the individual signal and background components on the main frame.
    total.plotOn(frame, RooFit::Components("signal"), RooFit::LineColor(kRed), 
                 RooFit::LineStyle(9), RooFit::LineWidth(10), RooFit::Name("signal_fit"));
    total.plotOn(frame, RooFit::Components("background"), RooFit::LineColor(kGreen + 2), 
                 RooFit::LineStyle(2), RooFit::LineWidth(10), RooFit::Name("background_fit"));
    // Plot the data on top so it remains visible.
    data.plotOn(frame, RooFit::Name("data"));
    
    // Adjust the vertical range of the frame.
    double max_data_point = frame->GetMaximum();
    frame->SetMaximum(1.1 * max_data_point);
    frame->SetMinimum(1);
    // Display the fitted parameter values on the plot.
    total.paramOn(frame);

    // ------------------------------
    // Set Up the Canvas and Draw the Plots
    // ------------------------------
    // Create a canvas divided into two pads: top pad for the main fit and bottom pad for the pull distribution.
    TCanvas* c = new TCanvas("c", "Fit Canvas", 900, 900);
    c->SetTitle(""); // Disable canvas title
    c->Divide(1, 2, 0.0, 0.0);

    // Draw the main fit plot in the top pad.
    c->cd(1);
    c->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
    frame->GetYaxis()->SetTitleFont(62);
    frame->GetYaxis()->SetLabelFont(62);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetLabelSize(0.06);
    frame->GetYaxis()->SetTitleOffset(1);
    frame->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame->Draw();

    // Add a legend in the top pad.
    TLegend* legend = new TLegend(0.15, 0.7, 0.4, 0.9);
    legend->SetTextFont(62);
    legend->AddEntry(frame->findObject("data"), "Data", "p");
    legend->AddEntry(frame->findObject("total_fit"), "Total", "l");
    legend->AddEntry(frame->findObject("signal_fit"), "Signal", "l");
    legend->AddEntry(frame->findObject("background_fit"), "Background", "l");
    legend->SetBorderSize(0);
    legend->Draw();

    // Draw the pull plot in the bottom pad.
    c->cd(2);
    c->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
    pull_frame->GetYaxis()->SetTitle("Pull");
    pull_frame->GetYaxis()->SetTitleSize(0.13);
    pull_frame->GetYaxis()->SetTitleOffset(0.3);
    pull_frame->GetYaxis()->SetLabelSize(0.10);
    pull_frame->GetXaxis()->SetLabelSize(0.13);
    pull_frame->GetXaxis()->SetTitleSize(0.13);
    pull_frame->GetXaxis()->SetTitleFont(62);
    pull_frame->GetYaxis()->SetTitleFont(62);
    pull_frame->GetYaxis()->SetLabelFont(62);
    pull_frame->GetXaxis()->SetLabelFont(62);
    pull_frame->GetYaxis()->SetNdivisions(5);
    pull_frame->GetYaxis()->SetRangeUser(-11, 11);
    pull_frame->GetXaxis()->SetNdivisions(6);
    pull_frame->Draw();

    // Add a horizontal zero-line to the pull plot for reference.
    TLine* zeroLine = new TLine(1820., 0.0, 1910., 0.0);
    zeroLine->SetLineColor(kBlack);
    zeroLine->SetLineStyle(9);
    zeroLine->SetLineWidth(2);
    zeroLine->Draw();

    // Print the minimum and maximum values of the pull distribution.
    double minPull = pull_frame->GetMinimum();
    double maxPull = pull_frame->GetMaximum();
    std::cout << "Minimum Pull: " << minPull << std::endl;
    std::cout << "Maximum Pull: " << maxPull << std::endl;
    
    // Draw the complete canvas.
    c->Draw();
}

