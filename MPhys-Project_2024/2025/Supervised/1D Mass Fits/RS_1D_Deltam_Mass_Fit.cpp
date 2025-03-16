// This script performs a 1D mass fit for the RS Δm distribution using a combination of
// a Johnson PDF, a bifurcated Gaussian, and two Gaussian components to model the signal,
// along with a power-law background modeled by a RooGenericPdf.
// The code reads data from a ROOT TTree, performs the fit, plots the results (including
// a pull distribution), and prints key fit statistics.
// To run the analysis, call the function RS_deltam_massfit_Johnson_and_Bifur_sig_powerlaw_bkg_1D_fit()
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

// Use the RooFit namespace for cleaner code
using namespace RooFit;

void RS_deltam_massfit_Johnson_and_Bifur_sig_powerlaw_bkg_1D_fit() {
    // Replace the absolute file path with a generic file name for GitHub distribution.
    const char* file_path = "input.root";
    
    // Open the input ROOT file and retrieve the TTree named "DecayTree"
    TFile* root_file = TFile::Open(file_path);
    TTree* tree = (TTree*)root_file->Get("DecayTree");

    // Define the Δm variable over the range [142, 152] MeV/c^2.
    RooRealVar deltam_refit("deltam_ReFit", "#Delta m [MeV/c^{2}]", 142.0, 152.0);
    
    // Create a RooDataSet from the TTree using the Δm variable.
    RooDataSet data("data", "dataset", tree, RooArgSet(deltam_refit));

    // ------------------------------
    // Signal Model: Johnson + Bifurcated Gaussian + Two Gaussian Components
    // ------------------------------
    // Define the common mean for the signal components.
    RooRealVar mean("#mu", "Signal mean", 143., 142.0, 152.0);
    
    // Bifurcated Gaussian parameters (modeling asymmetric tails).
    RooRealVar sigma_L("#sigma_{L}", "Left width of Bifurcated Gaussian", 0.23, 0., 3);
    RooRealVar sigma_R("#sigma_{R}", "Right width of Bifurcated Gaussian", 0.26, 0., 3);
    RooBifurGauss bifurGauss("bifurGauss", "Bifurcated Gaussian", deltam_refit, mean, sigma_L, sigma_R);
    
    // Gaussian component parameters.
    RooRealVar sigma("#sigma", "Signal width", 0.1, 0.0, 3);
    RooGaussian gauss("gauss", "Gaussian", deltam_refit, mean, sigma);
    
    // A second Gaussian component, with its own fraction.
    RooGaussian gauss1("gauss1", "Gaussian1", deltam_refit, mean, sigma);
    
    // Fraction parameters to weight the contributions of the Gaussian components.
    RooRealVar f_gauss("f_{gauss}", "Fraction of Gaussian", 0.3, 0.0, 1.0);
    RooRealVar f_gauss1("f_{gauss1}", "Fraction of Gaussian1", 0.3, 0.0, 1.0);
    
    // Fraction between the Johnson and Bifurcated Gaussian components.
    RooRealVar f_bifur("f_{bifur}", "Fraction of Bifurcated Gaussian", 0.7, 0.0, 1.0);
    
    // Define the Johnson PDF parameters for the signal component.
    RooRealVar lambda("#lambda", "deltam_lambda", 0.28, 0, 2);
    RooRealVar gamma("#gamma", "deltam_gamma", 0.18, -1, 2);
    RooRealVar delta("#delta", "deltam_delta", 0.9, 0, 2);
    RooJohnson deltam_johnson("deltam_johnson", "signal", deltam_refit, mean, lambda, gamma, delta);
    
    // Combine the signal components into one final signal PDF.
    // The ordering in RooArgList corresponds to the fractions provided in the RooArgList for yields.
    RooAddPdf signal("signal", "Johnson + Bifurcated Gaussian Signal", 
                     RooArgList(bifurGauss, deltam_johnson, gauss, gauss1),
                     RooArgList(f_gauss, f_gauss1, f_bifur));
    
    // ------------------------------
    // Background Model: Power-Law (Polynomial) Background using RooGenericPdf
    // ------------------------------
    // Define coefficients for the background model.
    RooRealVar alpha1("alpha1", "First coefficient", 0.1, -2.0, 2.0);
    RooRealVar alpha2("alpha2", "Second coefficient", 0.1, -2.0, 2.0);
    RooRealVar beta("beta", "Exponent coefficient", 0.75, 0.0, 1.0);
    
    // Construct the background PDF as a power-law polynomial function of (deltam_ReFit - 139.5).
    RooGenericPdf background(
        "background", "Polynomial background", 
        "alpha1 * pow((deltam_ReFit - 139.5), beta) + alpha2 * (deltam_ReFit - 139.5)", 
        RooArgSet(deltam_refit, alpha1, alpha2, beta)
    );
    
    // ------------------------------
    // Total PDF: Signal + Background
    // ------------------------------
    // Define the expected number of signal and background events.
    RooRealVar nsig("nsig", "Number of signal events", 4090000, 0, 1000000000);
    RooRealVar nbkg("nbkg", "Number of background events", 35000000, 0, 1000000000);
    
    // Combine the signal and background PDFs into a total PDF.
    RooAddPdf total("total", "Signal + Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // ------------------------------
    // Perform the Fit
    // ------------------------------
    // Fit the total PDF to the data. Use Minos for more accurate error estimation.
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print();

    // ------------------------------
    // Plot the Fit Results
    // ------------------------------
    // Create a frame for the Δm variable and plot the data.
    RooPlot* frame = deltam_refit.frame(RooFit::Title("Mass Fit"));
    data.plotOn(frame);
    total.plotOn(frame, RooFit::LineColor(kBlue), RooFit::LineWidth(8), RooFit::Name("total_fit"), Format("NEALU"));

    // Add chi-square per degree-of-freedom information to the plot.
    double chi2 = frame->chiSquare();
    int ndata = data.numEntries();
    int ndf = fit_result->floatParsFinal().getSize();
    int dof = ndata - ndf;  // Degrees of freedom
    TText* chi2_text = new TText(0.6, 0.8, Form("#chi^{2}/NDOF = %.2f/%d", chi2, dof));
    chi2_text->SetTextSize(0.04);
    chi2_text->SetTextFont(42);
    chi2_text->SetTextAlign(12); // Center alignment
    chi2_text->Draw();

    // ------------------------------
    // Create a Pull Plot
    // ------------------------------
    // Create a separate frame for the pull distribution.
    RooPlot* pull_frame = deltam_refit.frame(RooFit::Title("Pull Distribution"));
    // Get the pull histogram from the main frame and add it to the pull frame.
    RooHist* pullHist = frame->pullHist();
    pull_frame->addPlotable(pullHist, "P");

    // Overlay the individual signal and background components on the main frame.
    total.plotOn(frame, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::Name("signal_fit"));
    total.plotOn(frame, RooFit::Components("background"), RooFit::LineColor(kGreen + 2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::Name("background_fit"));
    // Plot the data on top to ensure it is visible.
    data.plotOn(frame, RooFit::Name("data"));
    
    // Adjust the vertical scale of the frame.
    double max_data_point = frame->GetMaximum();
    frame->SetMaximum(1.2 * max_data_point);
    frame->SetMinimum(1);

    // Display the fitted parameters on the plot.
    total.paramOn(frame);

    // ------------------------------
    // Set Up the Canvas and Draw the Plots
    // ------------------------------
    // Create a canvas divided into two pads: the top pad for the main fit plot and the bottom pad for the pull plot.
    TCanvas* c = new TCanvas("c", "Fit Canvas", 900, 900);
    c->SetTitle("");  // Disable canvas title
    c->Divide(1, 2, 0.0, 0.0);

    // Draw the main fit plot on the top pad.
    c->cd(1);
    c->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
    frame->GetYaxis()->SetTitleFont(62);
    frame->GetYaxis()->SetLabelFont(62);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetLabelSize(0.06);
    frame->GetYaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
    frame->Draw();

    // Add a legend in the top pad.
    TLegend* legend = new TLegend(0.15, 0.7, 0.4, 0.9);
    legend->SetTextFont(62);
    legend->AddEntry(frame->findObject("data"), "Data", "p");
    legend->AddEntry(frame->findObject("total_fit"), "Total", "l");
    legend->AddEntry(frame->findObject("signal_fit"), "Signal", "l");
    legend->AddEntry(frame->findObject("background_fit"), "Background", "l");
    legend->AddEntry(frame->findObject("gauss_fit"), "Added Gaussian", "l");
    legend->SetBorderSize(0);
    legend->Draw();

    // Draw the pull plot on the bottom pad.
    c->cd(2);
    c->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
    c->cd(2)->SetRightMargin(0.15);
    pull_frame->GetYaxis()->SetTitle("Pull");
    pull_frame->GetYaxis()->SetTitleSize(0.13); // Title size
    pull_frame->GetYaxis()->SetTitleOffset(0.3); // Title offset
    pull_frame->GetYaxis()->SetLabelSize(0.10);   // Y-axis label size
    pull_frame->GetXaxis()->SetLabelSize(0.13);   // X-axis label size
    pull_frame->GetXaxis()->SetTitleSize(0.13);     // X-axis title size
    pull_frame->GetXaxis()->SetTitleFont(62);
    pull_frame->GetYaxis()->SetTitleFont(62);
    pull_frame->GetYaxis()->SetLabelFont(62);
    pull_frame->GetXaxis()->SetLabelFont(62);
    pull_frame->GetYaxis()->SetRangeUser(-19, 19);
    pull_frame->GetYaxis()->SetNdivisions(4);
    pull_frame->GetXaxis()->SetNdivisions(6);
    pull_frame->Draw();

    // Add a horizontal zero-line to the pull plot for reference.
    TLine* zeroLine = new TLine(142., 0.0, 152., 0.0);
    zeroLine->SetLineColor(kBlack);
    zeroLine->SetLineStyle(9);
    zeroLine->SetLineWidth(2);
    zeroLine->Draw();

    // Print the minimum and maximum pull values.
    double minPull = pull_frame->GetMinimum();
    double maxPull = pull_frame->GetMaximum();
    std::cout << "Minimum Pull: " << minPull << std::endl;
    std::cout << "Maximum Pull: " << maxPull << std::endl;

    // Draw the complete canvas.
    c->Draw();
}

