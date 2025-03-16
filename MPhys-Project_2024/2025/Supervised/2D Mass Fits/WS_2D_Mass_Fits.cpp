// This script performs a two-dimensional mass fit for WS (Wrong-Sign) events
// by fitting the Δm and D0 refitted mass (Dst_ReFit_D0_M_best) distributions.
// It creates separate RooDataSets for the two mass variables, defines fixed signal
// and background PDFs for both variables, and then combines them into a total PDF
// for a 2D fit. The final plots include the mass fits and their pull distributions.
// To run the analysis, call the function twoD_Mass_Fit_WS_DeltaM_mass_projections()
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

// Function to perform a 2D mass fit with projections for Δm and D0 mass distributions
void twoD_Mass_Fit_WS_DeltaM_mass_projections() {
    // Open the ROOT file with the data (replace "your_data_file.root" with your actual file name)
    const char* file_path = "your_data_file.root";
    TFile* root_file = TFile::Open(file_path);
    // Retrieve the TTree from the file (replace "DecayTree" with your tree name if different)
    TTree* tree = (TTree*)root_file->Get("DecayTree");

    // Define the two observables (Δm and D0 mass) with their ranges
    RooRealVar deltam_refit("deltam_ReFit", "#Δm [MeV/c^{2}]", 142.0, 152.0);
    RooRealVar D0_refit("Dst_ReFit_D0_M_best", "m(D^{0}) [MeV/c^{2}]", 1820, 1910);

    // Create datasets from the tree using the defined variables
    RooDataSet data("data", "Data", tree, RooArgSet(D0_refit, deltam_refit));
    RooDataSet d0_data("d0_data", "Data", tree, RooArgSet(D0_refit));
    RooDataSet deltam_data("deltam_data", "Data", tree, RooArgSet(deltam_refit));

    //---------------------------------
    // Define Fixed PDFs for Signal and Background
    //---------------------------------
    // Fixed parameters for the Δm distribution
    RooRealVar deltam_mean("deltam_mean", "Δm Mean", 145.4252);
    RooRealVar deltam_lambda("deltam_lambda", "Δm Lambda", 0.2376);
    RooRealVar deltam_gamma("deltam_gamma", "Δm Gamma", 0.99);
    RooRealVar deltam_delta("deltam_delta", "Δm Delta", 0.0226);

    // Fixed parameters for the D0 mass distribution
    RooRealVar d0_mean("d0_mean", "D0 Mean", 1866.145);
    RooRealVar d0_lambda("d0_lambda", "D0 Lambda", 20.927);
    RooRealVar d0_gamma("d0_gamma", "D0 Gamma", 0.0898);
    RooRealVar d0_delta("d0_delta", "D0 Delta", 2.8323);
    
    // Define a Johnson PDF for the D0 signal shape
    RooJohnson d0_johnson("deltam_johnson", "signal", D0_refit, d0_mean, d0_lambda, d0_gamma, d0_delta);

    // Define a bifurcated Gaussian for the D0 distribution
    RooRealVar d0_sigma_L("#sigma_{L}", "Left width of Bifurcated Gaussian", 16.6);
    RooRealVar d0_sigma_R("#sigma_{R}", "Right width of Bifurcated Gaussian", 11.3);
    RooBifurGauss bifurGauss("bifurGauss", "Bifurcated Gaussian", D0_refit, d0_mean, d0_sigma_L, d0_sigma_R);
    // Fraction of the bifurcated Gaussian contribution
    RooRealVar f_bifur("f_{bifur}", "Fraction of Bifurcated Gaussian", 0.9966);

    // Define the Δm signal PDF using a Johnson function
    RooJohnson deltam_signal_pdf("deltam_signal_pdf", "Δm Signal PDF", deltam_refit, deltam_mean, deltam_lambda, deltam_gamma, deltam_delta);
    // Combine the D0 PDFs into a single signal PDF
    RooAddPdf d0_refit_signal_pdf("signal", "Double Gaussian + Bifurcated Gaussian Signal", 
                     RooArgList(d0_johnson, bifurGauss), RooArgList(f_bifur));

    // Parameters for the Δm background PDF (polynomial with exponential damping)
    RooRealVar alpha("alpha", "alpha", 0.46334);
    RooRealVar beta("beta", "beta", 0.01080);

    // Generic PDF for the Δm background shape
    RooGenericPdf deltam_background_pdf(
        "background", "Polynomial background", 
        "pow((deltam_ReFit - 139.5), alpha) * exp(-1 * beta * (deltam_ReFit - 139.5))",
        RooArgSet(deltam_refit,  alpha, beta)
    );
    
    // Define background for the D0 mass distribution using a second-order Chebyshev polynomial
    RooRealVar a1("a1", "First coefficient", 0.00692);
    RooRealVar a2("a2", "Second coefficient", -0.04276);
    RooChebychev d0_refit_background_pdf("d0_refit_background_pdf", "Second-order Chebyshev Background",
                        D0_refit, RooArgList(a1, a2));

    //-------------------------------
    // Combine Signal and Background PDFs
    //-------------------------------
    // Create PDFs for various components by taking products of the D0 and Δm PDFs
    RooProdPdf signal("signal", "signal", RooArgList(d0_refit_signal_pdf, deltam_signal_pdf));
    RooProdPdf combinatorial("combinatorial", "Combinatorial background", RooArgList(d0_refit_background_pdf, deltam_background_pdf));
    RooProdPdf random_slowpion("random_slowpion", "true D0 with incorrect slow pion", RooArgList(d0_refit_signal_pdf, deltam_background_pdf));
    RooProdPdf misreco_D0("misreco_D0", "incorrect D0 with true slow pion", RooArgList(d0_refit_background_pdf, deltam_signal_pdf));
    
    //-----------------------
    // Define Yields for Each PDF Component
    //-----------------------
    RooRealVar Signal_yield("Signal_yield", "Number of signal events for Δm", 100000, 0, 276000);
    RooRealVar Combinatorial_yield("Combinatorial_yield", "Number of background events for Δm", 100000, 0, 276000);
    RooRealVar Random_slowpion_yield("Random_slowpion_yield", "Number of signal events for Δm", 100000, 0, 276000);
    RooRealVar Misreco_yield("Misreco_yield", "Number of background events for Δm", 100000, 0, 276000);

    //--------------------------------------------------
    // Total PDF: Sum of all components weighted by their yields
    //--------------------------------------------------
    RooAddPdf total("total", "Signal + Background",
        RooArgList(signal, combinatorial, random_slowpion, misreco_D0),
        RooArgList(Signal_yield, Combinatorial_yield, Random_slowpion_yield, Misreco_yield)
    );

    // Fit the total PDF to the data and save the fit result
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1));
    fit_result->Print();

    // Define mass ranges for the D0 variable to create separate projections
    D0_refit.setRange("range1", 1820, 1860);
    D0_refit.setRange("range2", 1860, 1870);
    D0_refit.setRange("range3", 1870, 1910);

    // Create frames for Δm projections in different D0 mass ranges
    RooPlot* frame1 = deltam_refit.frame(RooFit::Title("Δm Fit Projection (1820 < m(D^{0}) < 1860)"));
    RooPlot* frame2 = deltam_refit.frame(RooFit::Title("Δm Fit Projection (1860 < m(D^{0}) < 1870)"));
    RooPlot* frame3 = deltam_refit.frame(RooFit::Title("Δm Fit Projection (1870 < m(D^{0}) < 1910)"));

    // Plot data and total fit on the frames (using projection to the specific D0 ranges)
    data.plotOn(frame1, RooFit::CutRange("range1"), RooFit::Name("data1"));
    total.plotOn(frame1, RooFit::ProjectionRange("range1"), RooFit::Name("total_fit1"));
    // Create a pull distribution for the first range
    RooPlot* pull_frame1 = deltam_refit.frame();
    RooHist* pullHist1 = frame1->pullHist("data1", "total_fit1");
    pull_frame1->addPlotable(pullHist1, "P");
    // Plot individual components for range1 with different styles and colors
    total.plotOn(frame1, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("signal_fit1"));
    total.plotOn(frame1, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("combinatorial_fit1"));
    total.plotOn(frame1, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("random_slowpion_fit1"));
    total.plotOn(frame1, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("misreco_fit1"));

    // Similar plotting for the other D0 mass ranges (range2 and range3)
    data.plotOn(frame2, RooFit::CutRange("range2"), RooFit::Name("data2"));
    total.plotOn(frame2, RooFit::ProjectionRange("range2"), RooFit::Name("total_fit2"));
    RooPlot* pull_frame2 = deltam_refit.frame();
    RooHist* pullHist2 = frame2->pullHist("data2", "total_fit2");
    pull_frame2->addPlotable(pullHist2, "P");
    total.plotOn(frame2, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("signal_fit2"));
    total.plotOn(frame2, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("combinatorial_fit2"));
    total.plotOn(frame2, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("random_slowpion_fit2"));
    total.plotOn(frame2, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("misreco_fit2"));

    data.plotOn(frame3, RooFit::CutRange("range3"), RooFit::Name("data3"));
    total.plotOn(frame3, RooFit::ProjectionRange("range3"), RooFit::Name("total_fit3"));
    RooPlot* pull_frame3 = deltam_refit.frame();
    RooHist* pullHist3 = frame3->pullHist("data3", "total_fit3");
    pull_frame3->addPlotable(pullHist3, "P");
    total.plotOn(frame3, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("signal_fit3"));
    total.plotOn(frame3, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("combinatorial_fit3"));
    total.plotOn(frame3, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("random_slowpion_fit3"));
    total.plotOn(frame3, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1),
                RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("misreco_fit3"));

    ////////////////////////////////////////////////////////////////////////////////
    // Create a canvas and pads to display the Δm projections with pull distributions
    ////////////////////////////////////////////////////////////////////////////////
    TCanvas* c = new TCanvas("c", "Δm Projections", 1800, 600);

    // Create pads for each D0 mass range projection
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.0, 0.33, 1.0);
    TPad* pad2 = new TPad("pad2", "pad2", 0.33, 0.0, 0.66, 1.0);
    TPad* pad3 = new TPad("pad3", "pad3", 0.66, 0.0, 1.0, 1.0);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();

    // Set logarithmic y-axis for better visualization of distributions
    pad1->SetLogy();
    pad2->SetLogy();
    pad3->SetLogy();

    /////////////////
    // Range1 plots (first D0 mass window)
    /////////////////
    pad1->cd();
    // Create sub-pads: one for the main frame and one for the pull distribution
    TPad* pad1_1 = new TPad("pad1_1", "pad1_1", 0.0, 0.25, 1.0, 1.0);
    TPad* pad1_2 = new TPad("pad1_2", "pad1_2", 0.0, 0.0, 1.0, 0.25);
    pad1_1->SetBottomMargin(0.02);
    pad1_2->SetTopMargin(0.02);
    pad1_1->Draw();
    pad1_2->Draw();

    pad1_1->cd();
    frame1->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
    frame1->GetYaxis()->SetTitleOffset(1.5);
    frame1->GetYaxis()->SetTitleFont(62);
    frame1->GetYaxis()->SetLabelFont(62);
    frame1->GetYaxis()->SetTitleSize(0.07);
    frame1->GetYaxis()->SetLabelSize(0.07);
    frame1->Draw();
    // Add a label indicating the D0 mass range for this plot
    TLatex* D0_left = new TLatex(0.18, 0.95, "m(D^{0}) #in [1820, 1860]MeV/c^{2}");
    D0_left->SetNDC();
    D0_left->SetTextFont(62);
    D0_left->SetTextSize(0.06);
    D0_left->SetTextColor(kBlack);
    D0_left->Draw();

    pad1_2->cd();
    pull_frame1->GetYaxis()->SetTitle("Pull");
    pull_frame1->GetYaxis()->SetTitleSize(0.15);
    pull_frame1->GetYaxis()->SetTitleOffset(0.4);
    pull_frame1->GetYaxis()->SetLabelSize(0.15);
    pull_frame1->GetXaxis()->SetLabelSize(0.20);
    pull_frame1->GetXaxis()->SetTitleSize(0.17);
    pull_frame1->GetXaxis()->SetTitleFont(62);
    pull_frame1->GetYaxis()->SetTitleFont(62);
    pull_frame1->GetYaxis()->SetLabelFont(62);
    pull_frame1->GetXaxis()->SetLabelFont(62);
    pull_frame1->GetYaxis()->SetRangeUser(-5, 5);
    pull_frame1->GetYaxis()->SetNdivisions(5);
    pull_frame1->GetXaxis()->SetNdivisions(8);
    pull_frame1->Draw();
    // Draw a horizontal zero line for the pull distribution
    TLine* zeroLine1 = new TLine(142.0, 0.0, 152.0, 0.0);
    zeroLine1->SetLineColor(kBlack);
    zeroLine1->SetLineStyle(9);
    zeroLine1->SetLineWidth(2);
    zeroLine1->Draw();

    // ... (Additional plotting code for other mass ranges and canvases would follow)
    
    // Finally, draw the canvas
    c->Draw();
}
