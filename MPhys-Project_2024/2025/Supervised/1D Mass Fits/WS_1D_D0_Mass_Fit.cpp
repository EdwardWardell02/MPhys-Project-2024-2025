// This script performs a 1D mass fit for the D0 refitted mass variable (Dst_ReFit_D0_M_best)
// using a combination of a Johnson signal model and a bifurcated Gaussian to model the signal,
// and a second-order Chebyshev polynomial to model the background.
// The script also produces fit plots and a pull distribution for quality assessment.

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

void Dst_ReFit_D0_M_best_Johnson_Bifur_sig_powerlaw_bkg_1D_WS_fit() {
    // Use a generic file path (update "input.root" to your actual file if needed)
    const char* file_path = "input.root";
    
    // Open the input ROOT file and retrieve the TTree named "outtree"
    TFile* root_file = TFile::Open(file_path);
    TTree* tree = (TTree*)root_file->Get("outtree");

    // ------------------------------
    // Define the mass variable and dataset
    // ------------------------------
    // Create a RooRealVar for the D0 refitted mass in the range [1820, 1910] MeV/c^2.
    RooRealVar D0_refit("Dst_ReFit_D0_M_best", "m(D^{0}) [MeV/c^{2}]", 1820., 1910.);
    // Create a RooDataSet from the TTree using the defined mass variable.
    RooDataSet data("data", "dataset", tree, RooArgSet(D0_refit));

    // ------------------------------
    // Define the Signal Model
    // ------------------------------
    // Signal is modeled using a Johnson PDF combined with a bifurcated Gaussian.
    // Set up the parameters for the Johnson function.
    RooRealVar mean("#mu", "Signal mean", 1865., 1863., 1867.);
    RooRealVar lambda("#lambda", "deltam_lambda", 19, 0, 30);
    RooRealVar gamma("#gamma", "deltam_gamma", 0.001, -1, 1);
    RooRealVar delta("#delta", "deltam_delta", 2.7, 0, 5);
    // Create the Johnson signal model using the mass variable and parameters.
    RooJohnson d0_johnson("deltam_johnson", "signal", D0_refit, mean, lambda, gamma, delta);

    // Define the bifurcated Gaussian component.
    RooRealVar sigma_L("#sigma_{L}", "Left width of Bifurcated Gaussian", 18, 10, 100);
    RooRealVar sigma_R("#sigma_{R}", "Right width of Bifurcated Gaussian", 12, 0, 20);
    RooBifurGauss bifurGauss("bifurGauss", "Bifurcated Gaussian", D0_refit, mean, sigma_L, sigma_R);

    // Define the fraction between the Johnson and bifurcated Gaussian components.
    RooRealVar f_bifur("f_{bifur}", "Fraction of Bifurcated Gaussian", 0.9, 0.0, 1.0);

    // Combine the two components into the final signal model.
    RooAddPdf signal("signal", "Double Gaussian + Bifurcated Gaussian Signal", 
                     RooArgList(d0_johnson, bifurGauss), RooArgList(f_bifur));

    // ------------------------------
    // Define the Background Model
    // ------------------------------
    // Here we use a second-order Chebyshev polynomial to model the background.
    RooRealVar a1("a1", "First coefficient", 0.33, -1.0, 1.0);
    RooRealVar a2("a2", "Second coefficient", 0.17, -1.0, 1.0);
    RooChebychev background("background", "Second-order Chebyshev Background",
                            D0_refit, RooArgList(a1, a2));

    // ------------------------------
    // Combine Signal and Background into the Total PDF
    // ------------------------------
    // Define the yield parameters for signal and background.
    RooRealVar nsig("nsig", "Number of signal events", 4090000, 0, 100000000);
    RooRealVar nbkg("nbkg", "Number of background events", 35000000, 0, 1000000000);
    // Combine the signal and background PDFs using their yields.
    RooAddPdf total("total", "Signal + Background", RooArgList(signal, background), RooArgList(nsig, nbkg));

    // ------------------------------
    // Perform the Fit
    // ------------------------------
    // Fit the total model to the data, saving the result and using Minos for error estimation.
    RooFitResult* fit_result = total.fitTo(data, RooFit::Save(), RooFit::PrintLevel(-1), RooFit::Minos(true));
    fit_result->Print();

    // ------------------------------
    // Plot the Fit Results
    // ------------------------------
    // Create a frame for the mass variable and plot the data and total fit.
    RooPlot* frame = D0_refit.frame(RooFit::Title("Mass Fit"));
    data.plotOn(frame);
    total.plotOn(frame, RooFit::LineColor(kBlue), RooFit::LineWidth(10), RooFit::Name("total_fit"));

    // ------------------------------
    // Add Chi2/NDOF Text to the Plot
    // ------------------------------
    // Calculate the chi-square and degrees of freedom.
    double chi2 = frame->chiSquare();
    int ndata = data.numEntries();
    int ndf = fit_result->floatParsFinal().getSize();
    int dof = ndata - ndf;  // Degrees of freedom
    // Create a TText object to display the chi-square per degree of freedom.
    TText* chi2_text = new TText(0.6, 0.8, Form("#chi^{2}/NDOF = %.2f/%d", chi2, dof));
    chi2_text->SetTextSize(0.04);
    chi2_text->SetTextFont(42);
    chi2_text->SetTextAlign(12); // Center alignment
    chi2_text->Draw();

    // ------------------------------
    // Create the Pull Plot
    // ------------------------------
    // Create a frame for the pull distribution.
    RooPlot* pull_frame = D0_refit.frame(RooFit::Title("Pull Distribution"));
    // Obtain the pull histogram from the main frame.
    RooHist* pullHist = frame->pullHist();
    pull_frame->addPlotable(pullHist, "P");

    // Overlay the signal and background components on the main frame.
    total.plotOn(frame, RooFit::Components("signal"), RooFit::LineColor(kRed),
                 RooFit::LineStyle(9), RooFit::LineWidth(10), RooFit::Name("signal_fit"));
    total.plotOn(frame, RooFit::Components("background"), RooFit::LineColor(kGreen + 2),
                 RooFit::LineStyle(2), RooFit::LineWidth(10), RooFit::Name("background_fit"));
    // Plot data on top to ensure it remains visible.
    data.plotOn(frame, RooFit::Name("data"));
    
    // Adjust the frame's vertical range.
    double max_data_point = frame->GetMaximum();
    frame->SetMaximum(1.1 * max_data_point);
    frame->SetMinimum(1);
    // Display the parameter values on the frame.
    total.paramOn(frame);

    // ------------------------------
    // Set Up the Canvas and Draw the Plots
    // ------------------------------
    // Create a canvas divided into two pads: the top for the main fit and the bottom for the pull plot.
    TCanvas* c = new TCanvas("c", "Fit Canvas", 900, 900);
    c->SetTitle("");  // Disable canvas title
    c->Divide(1, 2, 0.0, 0.0);

    // Top pad for the main fit plot.
    c->cd(1);
    c->cd(1)->SetPad(0.05, 0.35, 0.95, 1.0);
    frame->GetYaxis()->SetTitleFont(62);
    frame->GetYaxis()->SetLabelFont(62);
    frame->GetYaxis()->SetTitleSize(0.06);
    frame->GetYaxis()->SetLabelSize(0.06);
    frame->GetYaxis()->SetTitleOffset(1.4);
    frame->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
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

    // Bottom pad for the pull plot.
    c->cd(2);
    c->cd(2)->SetPad(0.05, 0.0, 0.95, 0.35);
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
    pull_frame->GetYaxis()->SetRangeUser(-5, 5);
    pull_frame->GetYaxis()->SetNdivisions(6);
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
    
    // Draw the canvas.
    c->Draw();
}
