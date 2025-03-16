// This script performs a 1D mass fit for ΔM using RooFit.
// It uses a Johnson signal model and a polynomial background model (via RooGenericPdf)
// to describe the mass distribution from data in a ROOT file.
// The script also produces fit plots, pull distributions, and calculates the integrated yields
// in a specified mass range to estimate purity.
//
// To run the analysis, call the function 
// deltam_massfit_johnson_sig_polynomial_bkg_mass_fit_1D_WS_fit()
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

void deltam_massfit_johnson_sig_polynomial_bkg_mass_fit_1D_WS_fit() {
    // Set the file path to a generic input ROOT file (update as needed)
    const char* file_path = "input.root";
    
    // Open the ROOT file
    TFile* root_file = TFile::Open(file_path);
    if (!root_file || root_file->IsZombie()) {
        std::cerr << "Error: Could not open file " << file_path << std::endl;
        return;
    }
    
    // Retrieve the TTree named "outtree" from the file
    TTree* tree = (TTree*)root_file->Get("outtree");
    if (!tree) {
        std::cerr << "Error: Could not find tree in file " << file_path << std::endl;
        return;
    }
    
    // Get and print the total number of events in the tree
    int number_of_events = tree->GetEntries();
    std::cout << "Number of events in the tree: " << number_of_events << std::endl;

    // Define the RooRealVar for the ΔM mass variable in the range [142, 152] MeV/c^2
    RooRealVar deltam_refit("deltam_ReFit", "#Delta m [MeV/c^{2}]", 142., 152.);

    // Create a RooDataSet from the TTree using the defined mass variable
    RooDataSet data("data", "dataset", tree, RooArgSet(deltam_refit));

    // ------------------------------
    // Define the Signal Model: Johnson PDF
    // ------------------------------
    // Set up the parameters for the Johnson function
    RooRealVar mean("#mu", "Signal mean", 145, 144, 146);
    RooRealVar lambda("#lambda", "deltam_lambda", 0.23, 0, 10);
    RooRealVar gamma("#gamma", "deltam_gamma", 0.99, -1, 2);
    RooRealVar delta("#delta", "deltam_delta", 0.0266, 0, 3);
    
    // Create the Johnson signal model using the mass variable and parameters
    RooJohnson signal("signal", "signal", deltam_refit, mean, lambda, gamma, delta);

    // ------------------------------
    // Define the Background Model: Polynomial using RooGenericPdf
    // ------------------------------
    // Set up the parameters for the background model
    RooRealVar alpha("alpha", "alpha", 1.0, 0., 10.);
    RooRealVar beta("beta", "beta", 0.1, -10., 10.);
    
    // Define the background model as a polynomial function of (deltam_ReFit - 139.5)
    RooGenericPdf background(
        "background", "Polynomial background", 
        "pow((deltam_ReFit - 139.5), alpha) * exp(-1 * beta * (deltam_ReFit - 139.5))",
        RooArgSet(deltam_refit, alpha, beta)
    );

    // ------------------------------
    // Define the Total PDF: Signal + Background
    // ------------------------------
    // Define yield variables for signal and background components
    RooRealVar nsig("nsig", "Number of signal events", 13023, 0, 100000000);
    RooRealVar nbkg("nbkg", "Number of background events", 274774, 0, 100000000);

    // Combine the signal and background models into a total model
    RooAddPdf total("total", "Signal + Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // ------------------------------
    // Perform the Fit
    // ------------------------------
    // Fit the total model to the data and save the fit result.
    // RooFit::Minos(true) enables a more accurate error estimation.
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print();

    // ------------------------------
    // Plot the Fit Results
    // ------------------------------
    // Create a frame for the mass variable and plot the data and the total fit model.
    RooPlot* frame = deltam_refit.frame(RooFit::Title("Mass Fit"));
    data.plotOn(frame);
    total.plotOn(frame, RooFit::LineColor(kBlue), RooFit::LineWidth(8), RooFit::Name("total_fit"));

    // ------------------------------
    // Add Chi2/NDOF Text to the Plot
    // ------------------------------
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
    // Create a separate frame for the pull distribution and add the pull histogram
    RooPlot* pull_frame = deltam_refit.frame(RooFit::Title("Pull Distribution"));
    RooHist* pullHist = frame->pullHist();
    pull_frame->addPlotable(pullHist, "P");

    // Overlay the individual signal and background components on the main frame
    total.plotOn(frame, RooFit::Components("signal"), RooFit::LineColor(kRed),
                 RooFit::LineWidth(8), RooFit::LineStyle(9), RooFit::Name("signal_fit"));
    total.plotOn(frame, RooFit::Components("background"), RooFit::LineColor(kGreen + 2),
                 RooFit::LineWidth(8), RooFit::LineStyle(2), RooFit::Name("background_fit"));
    total.paramOn(frame); // Automatically print parameter values on the frame

    // ------------------------------
    // Set Up the Canvas for Plotting
    // ------------------------------
    // Create a canvas divided into two pads: one for the main fit and one for the pull plot.
    TCanvas* c = new TCanvas("c", "Fit Canvas", 900, 900);
    c->SetTitle("");  // Disable canvas title
    c->Divide(1, 2, 0.0, 0.0);

    // Draw the main fit frame in the top pad.
    c->cd(1);
    c->cd(1)->SetPad(0, 0.35, 0.95, 1.0);
    frame->GetYaxis()->SetTitleFont(62);
    frame->GetYaxis()->SetLabelFont(62);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetLabelSize(0.06);
    frame->GetYaxis()->SetTitleOffset(1.3);
    frame->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
    double max_data_point = frame->GetMaximum();
    frame->SetMaximum(1.2 * max_data_point);
    frame->Draw();

    // Add a legend to the top pad.
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
    pull_frame->GetYaxis()->SetTitleOffset(0.25);
    pull_frame->GetYaxis()->SetLabelSize(0.10);
    pull_frame->GetXaxis()->SetLabelSize(0.13);
    pull_frame->GetXaxis()->SetTitleSize(0.13);
    pull_frame->GetXaxis()->SetTitleFont(62);
    pull_frame->GetYaxis()->SetTitleFont(62);
    pull_frame->GetYaxis()->SetLabelFont(62);
    pull_frame->GetXaxis()->SetLabelFont(62);
    pull_frame->GetYaxis()->SetRangeUser(-5, 5);
    pull_frame->GetYaxis()->SetNdivisions(6);
    pull_frame->GetXaxis()->SetNdivisions(6);
    pull_frame->Draw();

    // Add a horizontal zero-line to the pull plot for reference.
    TLine* zeroLine = new TLine(142.0, 0.0, 152.0, 0.0);
    zeroLine->SetLineColor(kBlack);
    zeroLine->SetLineStyle(9);
    zeroLine->SetLineWidth(2);
    zeroLine->Draw();

    // ------------------------------
    // Calculate Integrated Signal and Background Yields
    // ------------------------------
    // Define the integration range for the signal (e.g., [144.5, 146] MeV/c^2)
    double rangeMin = 144.5;
    double rangeMax = 146.0;
    deltam_refit.setRange("signalRange", rangeMin, rangeMax);

    // Integrate the signal PDF over the defined range and scale by the signal yield.
    RooAbsReal* signalIntegral = signal.createIntegral(deltam_refit, RooFit::NormSet(deltam_refit), RooFit::Range("signalRange"));
    double signalIntegralValue = signalIntegral->getVal() * nsig.getVal();
    std::cout << "Signal events in the interval [" << rangeMin << ", " << rangeMax << "] MeV/c^2: " << signalIntegralValue << std::endl;

    // Integrate the background PDF over the defined range and scale by the background yield.
    RooAbsReal* backgroundIntegral = background.createIntegral(deltam_refit, RooFit::NormSet(deltam_refit), RooFit::Range("signalRange"));
    double backgroundIntegralValue = backgroundIntegral->getVal() * nbkg.getVal();
    std::cout << "Background events in the interval [" << rangeMin << ", " << rangeMax << "] MeV/c^2: " << backgroundIntegralValue << std::endl;

    // Calculate and print the purity in the specified mass range.
    double purity = signalIntegralValue / (backgroundIntegralValue);
    std::cout << "Purity in the interval [" << rangeMin << ", " << rangeMax << "] MeV/c^2: " << purity << std::endl;

    // Update and display the canvas.
    c->Draw();
}
