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

    D0_refit.setRange("range1", 1820, 1860);
    D0_refit.setRange("range2", 1860, 1870);
    D0_refit.setRange("range3", 1870, 1910);

    RooPlot* frame1 = deltam_refit.frame(RooFit::Title("DeltaM Fit Projection (1820 < m(D^{0}) < 1850)"));
    RooPlot* frame2 = deltam_refit.frame(RooFit::Title("DeltaM Fit Projection (1850 < m(D^{0}) < 1880)"));
    RooPlot* frame3 = deltam_refit.frame(RooFit::Title("DeltaM Fit Projection (1880 < m(D^{0}) < 1910)"));


    data.plotOn(frame1, RooFit::CutRange("range1"), RooFit::Name("data1"));
    total.plotOn(frame1, RooFit::ProjectionRange("range1"), RooFit::Name("total_fit1"));
    RooPlot* pull_frame1 = deltam_refit.frame();
    RooHist* pullHist1 = frame1->pullHist("data1", "total_fit1");
    pull_frame1->addPlotable(pullHist1, "P");
    total.plotOn(frame1, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("signal_fit1"));
    total.plotOn(frame1, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("combinatorial_fit1"));
    total.plotOn(frame1, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("random_slowpion_fit1"));
    total.plotOn(frame1, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range1"), RooFit::Name("misreco_fit1"));

    data.plotOn(frame2, RooFit::CutRange("range2"), RooFit::Name("data2"));
    total.plotOn(frame2, RooFit::ProjectionRange("range2"), RooFit::Name("total_fit2"));
    RooPlot* pull_frame2 = deltam_refit.frame();
    RooHist* pullHist2 = frame2->pullHist("data2", "total_fit2");
    pull_frame2->addPlotable(pullHist2, "P");
    total.plotOn(frame2, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("signal_fit2"));
    total.plotOn(frame2, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("combinatorial_fit2"));
    total.plotOn(frame2, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("random_slowpion_fit2"));
    total.plotOn(frame2, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range2"), RooFit::Name("misreco_fit2"));

    data.plotOn(frame3, RooFit::CutRange("range3"), RooFit::Name("data3"));
    total.plotOn(frame3, RooFit::ProjectionRange("range3"), RooFit::Name("total_fit3"));
    RooPlot* pull_frame3 = deltam_refit.frame();
    RooHist* pullHist3 = frame3->pullHist("data3", "total_fit3");
    pull_frame3->addPlotable(pullHist3, "P");
    total.plotOn(frame3, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("signal_fit3"));
    total.plotOn(frame3, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("combinatorial_fit3"));
    total.plotOn(frame3, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("random_slowpion_fit3"));
    total.plotOn(frame3, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range3"), RooFit::Name("misreco_fit3"));


    // Create a canvas to display the plots
    TCanvas* c = new TCanvas("c", "DeltaM Projections", 1800, 600);

    // Create pads for each range
    TPad* pad1 = new TPad("pad1", "pad1", 0.0, 0.0, 0.33, 1.0);
    TPad* pad2 = new TPad("pad2", "pad2", 0.33, 0.0, 0.66, 1.0);
    TPad* pad3 = new TPad("pad3", "pad3", 0.66, 0.0, 1.0, 1.0);

    pad1->Draw();
    pad2->Draw();
    pad3->Draw();

    pad1->SetLogy();
    pad2->SetLogy();
    pad3->SetLogy();
    // Plot for range1
    pad1->cd();
    TPad* pad1_1 = new TPad("pad1_1", "pad1_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad1_2 = new TPad("pad1_2", "pad1_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad1_1->SetBottomMargin(0.02);
    pad1_2->SetTopMargin(0.02);
    pad1_1->Draw();
    pad1_2->Draw();

    pad1_1->cd();
    frame1->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
    frame1->GetYaxis()->SetTitleOffset(1.3); // Adjusted
    frame1->GetYaxis()->SetTitleFont(62);
    frame1->GetYaxis()->SetLabelFont(62);
    frame1->GetYaxis()->SetTitleSize(0.07);
    frame1->GetYaxis()->SetLabelSize(0.07);
    frame1->Draw();
    TLatex* D0_left = new TLatex(0.18, 0.95, "m(D^{0}) #in [1820, 1860]MeV/c^{2}");
    D0_left->SetNDC();
    D0_left->SetTextFont(62);
    D0_left->SetTextSize(0.06);
    D0_left->SetTextColor(kBlack);
    D0_left->Draw();

    pad1_2->cd();
    pull_frame1->GetYaxis()->SetTitle("Pull");
    pull_frame1->GetYaxis()->SetTitleSize(0.15); // Adjusted
    pull_frame1->GetYaxis()->SetTitleOffset(0.4); // Adjusted
    pull_frame1->GetYaxis()->SetLabelSize(0.15); // Adjusted
    pull_frame1->GetXaxis()->SetLabelSize(0.20); // Adjusted
    pull_frame1->GetXaxis()->SetTitleSize(0.17); // Adjusted
    pull_frame1->GetXaxis()->SetTitleFont(62);
    pull_frame1->GetYaxis()->SetTitleFont(62);
    pull_frame1->GetYaxis()->SetLabelFont(62);
    pull_frame1->GetXaxis()->SetLabelFont(62);
    pull_frame1->GetYaxis()->SetRangeUser(-25, 25);
    pull_frame1->GetYaxis()->SetNdivisions(4);
    pull_frame1->GetXaxis()->SetNdivisions(8);
    pull_frame1->Draw();
    TLine* zeroLine1 = new TLine(142.0, 0.0, 152.0, 0.0);
    zeroLine1->SetLineColor(kBlack);
    zeroLine1->SetLineStyle(9);
    zeroLine1->SetLineWidth(2);
    zeroLine1->Draw();

    // Plot for range2
    pad2->cd();
    TPad* pad2_1 = new TPad("pad2_1", "pad2_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad2_2 = new TPad("pad2_2", "pad2_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad2_1->SetBottomMargin(0.02);
    pad2_2->SetTopMargin(0.02);
    pad2_1->Draw();
    pad2_2->Draw();

    pad2_1->cd();
    frame2->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
    frame2->GetYaxis()->SetTitleOffset(1.3); // Unchanged
    frame2->GetYaxis()->SetTitleFont(62);
    frame2->GetYaxis()->SetLabelFont(62);
    frame2->GetYaxis()->SetTitleSize(0.07);
    frame2->GetYaxis()->SetLabelSize(0.07);
    frame2->Draw();
    TLatex* D0_peak = new TLatex(0.17, 0.95, "m(D^{0}) #in [1860, 1870]MeV/c^{2}");
    D0_peak->SetNDC();
    D0_peak->SetTextFont(62);
    D0_peak->SetTextSize(0.06);
    D0_peak->SetTextColor(kBlack);
    D0_peak->Draw();

    pad2_2->cd();
    pull_frame2->GetYaxis()->SetTitle("Pull");
    pull_frame2->GetYaxis()->SetTitleSize(0.15); // Adjusted
    pull_frame2->GetYaxis()->SetTitleOffset(0.4); // Adjusted
    pull_frame2->GetYaxis()->SetLabelSize(0.15); // Adjusted
    pull_frame2->GetXaxis()->SetLabelSize(0.20); // Adjusted
    pull_frame2->GetXaxis()->SetTitleSize(0.17); // Adjusted
    pull_frame2->GetXaxis()->SetTitleFont(62);
    pull_frame2->GetYaxis()->SetTitleFont(62);
    pull_frame2->GetYaxis()->SetLabelFont(62);
    pull_frame2->GetXaxis()->SetLabelFont(62);
    pull_frame2->GetYaxis()->SetRangeUser(-25, 25);
    pull_frame2->GetYaxis()->SetNdivisions(4);
    pull_frame2->GetXaxis()->SetNdivisions(8);
    pull_frame2->Draw();
    TLine* zeroLine2 = new TLine(142.0, 0.0, 152.0, 0.0);
    zeroLine2->SetLineColor(kBlack);
    zeroLine2->SetLineStyle(9);
    zeroLine2->SetLineWidth(2);
    zeroLine2->Draw();

    // Plot for range3
    pad3->cd();
    TPad* pad3_1 = new TPad("pad3_1", "pad3_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad3_2 = new TPad("pad3_2", "pad3_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad3_1->SetBottomMargin(0.02);
    pad3_2->SetTopMargin(0.02);
    pad3_1->Draw();
    pad3_2->Draw();

    pad3_1->cd();
    frame3->GetYaxis()->SetTitle("Events / (0.1 MeV/c^{2})");
    frame3->GetYaxis()->SetTitleOffset(1.3); // Adjusted
    frame3->GetYaxis()->SetTitleFont(62);
    frame3->GetYaxis()->SetLabelFont(62);
    frame3->GetYaxis()->SetTitleSize(0.07);
    frame3->GetYaxis()->SetLabelSize(0.07);
    frame3->Draw();
    TLatex* D0_right = new TLatex(0.17, 0.95, "m(D^{0}) #in [1870, 1910]MeV/c^{2}");
    D0_right->SetNDC();
    D0_right->SetTextFont(62);
    D0_right->SetTextSize(0.06);
    D0_right->SetTextColor(kBlack);
    D0_right->Draw();

    pad3_2->cd();
    pull_frame3->GetYaxis()->SetTitle("Pull");
    pull_frame3->GetYaxis()->SetTitleSize(0.15); // Adjusted
    pull_frame3->GetYaxis()->SetTitleOffset(0.4); // Adjusted
    pull_frame3->GetYaxis()->SetLabelSize(0.15); // Adjusted
    pull_frame3->GetXaxis()->SetLabelSize(0.20); // Adjusted
    pull_frame3->GetXaxis()->SetTitleSize(0.17); // Adjusted
    pull_frame3->GetXaxis()->SetTitleFont(62);
    pull_frame3->GetYaxis()->SetTitleFont(62);
    pull_frame3->GetYaxis()->SetLabelFont(62);
    pull_frame3->GetXaxis()->SetLabelFont(62);
    pull_frame3->GetYaxis()->SetRangeUser(-25, 25);
    pull_frame3->GetYaxis()->SetNdivisions(4);
    pull_frame3->GetXaxis()->SetNdivisions(8);
    pull_frame3->Draw();
    TLine* zeroLine3 = new TLine(142.0, 0.0, 152.0, 0.0);
    zeroLine3->SetLineColor(kBlack);
    zeroLine3->SetLineStyle(9);
    zeroLine3->SetLineWidth(2);
    zeroLine3->Draw();

    c->Draw();

    // Define ranges for deltam
    deltam_refit.setRange("range4", 142.0, 144.5);
    deltam_refit.setRange("range5", 144.5, 146.5);
    deltam_refit.setRange("range6", 146.5, 152.0);

    RooPlot* frame4 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (142 < #Delta m < 144.5)"));
    RooPlot* frame5 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (144.5 < #Delta m < 146.5)"));
    RooPlot* frame6 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (146.5 < #Delta m < 152)"));

    data.plotOn(frame4, RooFit::CutRange("range4"), RooFit::Name("data4"));
    total.plotOn(frame4, RooFit::ProjectionRange("range4"), RooFit::Name("total_fit4"));
    RooPlot* pull_frame4 = D0_refit.frame(RooFit::Title("Pull Distribution (142 < #Delta m < 144.5)"));
    RooHist* pullHist4 = frame4->pullHist("data4", "total_fit4");
    pull_frame4->addPlotable(pullHist4, "P");
    total.plotOn(frame4, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range4"), RooFit::Name("signal_fit4"));
    total.plotOn(frame4, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range4"), RooFit::Name("combinatorial_fit4"));
    total.plotOn(frame4, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range4"), RooFit::Name("random_slowpion_fit4"));
    total.plotOn(frame4, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range4"), RooFit::Name("misreco_fit4"));

    data.plotOn(frame5, RooFit::CutRange("range5"), RooFit::Name("data5"));
    total.plotOn(frame5, RooFit::ProjectionRange("range5"), RooFit::Name("total_fit5"));
    RooPlot* pull_frame5 = D0_refit.frame(RooFit::Title("Pull Distribution (144.5 < #Delta m < 146.5)"));
    RooHist* pullHist5 = frame5->pullHist("data5", "total_fit5");
    pull_frame5->addPlotable(pullHist5, "P");
    total.plotOn(frame5, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range5"), RooFit::Name("signal_fit5"));
    total.plotOn(frame5, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range5"), RooFit::Name("combinatorial_fit5"));
    total.plotOn(frame5, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range5"), RooFit::Name("random_slowpion_fit5"));
    total.plotOn(frame5, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range5"), RooFit::Name("misreco_fit5"));

    data.plotOn(frame6, RooFit::CutRange("range6"), RooFit::Name("data6"));
    total.plotOn(frame6, RooFit::ProjectionRange("range6"), RooFit::Name("total_fit6"));
    RooPlot* pull_frame6 = D0_refit.frame(RooFit::Title("Pull Distribution (146.5 < #Delta m < 152)"));
    RooHist* pullHist6 = frame6->pullHist("data6", "total_fit6");
    pull_frame6->addPlotable(pullHist6, "P");
    total.plotOn(frame6, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range6"), RooFit::Name("signal_fit6"));
    total.plotOn(frame6, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range6"), RooFit::Name("combinatorial_fit6"));
    total.plotOn(frame6, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range6"), RooFit::Name("random_slowpion_fit6"));
    total.plotOn(frame6, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range6"), RooFit::Name("misreco_fit6"));

    // Create a canvas to display the plots for D0 mass projections
    TCanvas* c2 = new TCanvas("c2", "D0 Mass Projections", 1800, 600);

    // Create pads for each range
    TPad* pad4 = new TPad("pad4", "pad4", 0.0, 0.0, 0.33, 1.0);
    TPad* pad5 = new TPad("pad5", "pad5", 0.33, 0.0, 0.66, 1.0);
    TPad* pad6 = new TPad("pad6", "pad6", 0.66, 0.0, 1.0, 1.0);

    pad4->Draw();
    pad5->Draw();
    pad6->Draw();

    pad4->SetLogy();
    pad5->SetLogy();
    pad6->SetLogy();
    // Plot for range4
    pad4->cd();
    TPad* pad4_1 = new TPad("pad4_1", "pad4_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad4_2 = new TPad("pad4_2", "pad4_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad4_1->SetBottomMargin(0.02);
    pad4_2->SetTopMargin(0.02);
    pad4_1->Draw();
    pad4_2->Draw();

    pad4_1->cd();
    frame4->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame4->GetYaxis()->SetTitleOffset(1.5); // Adjusted
    frame4->GetYaxis()->SetTitleFont(62);
    frame4->GetYaxis()->SetLabelFont(62);
    frame4->GetYaxis()->SetTitleSize(0.06);
    frame4->GetYaxis()->SetLabelSize(0.05);
    frame4->Draw();
    TLatex* deltam_left = new TLatex(0.18, 0.95, "#Delta m #in [142, 144.5]MeV/c^{2}");
    deltam_left->SetNDC();
    deltam_left->SetTextFont(62);
    deltam_left->SetTextSize(0.06);
    deltam_left->SetTextColor(kBlack);
    deltam_left->Draw();

    pad4_2->cd();
    pull_frame4->GetYaxis()->SetTitle("Pull");
    pull_frame4->GetYaxis()->SetTitleSize(0.15); // Adjusted
    pull_frame4->GetYaxis()->SetTitleOffset(0.4); // Adjusted
    pull_frame4->GetYaxis()->SetLabelSize(0.15); // Adjusted
    pull_frame4->GetXaxis()->SetLabelSize(0.20); // Adjusted
    pull_frame4->GetXaxis()->SetTitleSize(0.17); // Adjusted
    pull_frame4->GetXaxis()->SetTitleFont(62);
    pull_frame4->GetYaxis()->SetTitleFont(62);
    pull_frame4->GetYaxis()->SetLabelFont(62);
    pull_frame4->GetXaxis()->SetLabelFont(62);
    pull_frame4->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame4->GetYaxis()->SetNdivisions(4);
    pull_frame4->GetXaxis()->SetNdivisions(8);
    pull_frame4->Draw();
    TLine* zeroLine4 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine4->SetLineColor(kBlack);
    zeroLine4->SetLineStyle(9);
    zeroLine4->SetLineWidth(2);
    zeroLine4->Draw();

    // Plot for range5
    pad5->cd();
    TPad* pad5_1 = new TPad("pad5_1", "pad5_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad5_2 = new TPad("pad5_2", "pad5_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad5_1->SetBottomMargin(0.02);
    pad5_2->SetTopMargin(0.02);
    pad5_1->Draw();
    pad5_2->Draw();

    pad5_1->cd();
    frame5->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame5->GetYaxis()->SetTitleOffset(1.4); // Unchanged
    frame5->GetYaxis()->SetTitleFont(62);
    frame5->GetYaxis()->SetLabelFont(62);
    frame5->GetYaxis()->SetTitleSize(0.06);
    frame5->GetYaxis()->SetLabelSize(0.07);
    frame5->Draw();
    TLatex* deltam_peak = new TLatex(0.17, 0.95, "#Delta m #in [144.5, 146.5]MeV/c^{2}");
    deltam_peak->SetNDC();
    deltam_peak->SetTextFont(62);
    deltam_peak->SetTextSize(0.06);
    deltam_peak->SetTextColor(kBlack);
    deltam_peak->Draw();

    pad5_2->cd();
    pull_frame5->GetYaxis()->SetTitle("Pull");
    pull_frame5->GetYaxis()->SetTitle("Pull");
    pull_frame5->GetYaxis()->SetTitleSize(0.15); // Adjusted
    pull_frame5->GetYaxis()->SetTitleOffset(0.4); // Adjusted
    pull_frame5->GetYaxis()->SetLabelSize(0.15); // Adjusted
    pull_frame5->GetXaxis()->SetLabelSize(0.20); // Adjusted
    pull_frame5->GetXaxis()->SetTitleSize(0.17); // Adjusted
    pull_frame5->GetXaxis()->SetTitleFont(62);
    pull_frame5->GetYaxis()->SetTitleFont(62);
    pull_frame5->GetYaxis()->SetLabelFont(62);
    pull_frame5->GetXaxis()->SetLabelFont(62);
    pull_frame5->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame5->GetYaxis()->SetNdivisions(4);
    pull_frame5->GetXaxis()->SetNdivisions(8);
    pull_frame5->Draw();
    TLine* zeroLine5 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine5->SetLineColor(kBlack);
    zeroLine5->SetLineStyle(9);
    zeroLine5->SetLineWidth(2);
    zeroLine5->Draw();

    // Plot for range6
    pad6->cd();
    TPad* pad6_1 = new TPad("pad6_1", "pad6_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad6_2 = new TPad("pad6_2", "pad6_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad6_1->SetBottomMargin(0.02);
    pad6_2->SetTopMargin(0.02);
    pad6_1->Draw();
    pad6_2->Draw();

    pad6_1->cd();
    frame6->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame6->GetYaxis()->SetTitleOffset(1.5); // Adjusted
    frame6->GetYaxis()->SetTitleFont(62);
    frame6->GetYaxis()->SetLabelFont(62);
    frame6->GetYaxis()->SetTitleSize(0.06);
    frame6->GetYaxis()->SetLabelSize(0.05);
    frame6->Draw();
    TLatex* deltam_right = new TLatex(0.17, 0.95, "#Delta m #in [146.5, 152]MeV/c^{2}");
    deltam_right->SetNDC();
    deltam_right->SetTextFont(62);
    deltam_right->SetTextSize(0.06);
    deltam_right->SetTextColor(kBlack);
    deltam_right->Draw();

    pad6_2->cd();
    pull_frame6->GetYaxis()->SetTitle("Pull");
    pull_frame6->GetYaxis()->SetTitleSize(0.15); // Adjusted
    pull_frame6->GetYaxis()->SetTitleOffset(0.4); // Adjusted
    pull_frame6->GetYaxis()->SetLabelSize(0.15); // Adjusted
    pull_frame6->GetXaxis()->SetLabelSize(0.20); // Adjusted
    pull_frame6->GetXaxis()->SetTitleSize(0.17); // Adjusted
    pull_frame6->GetXaxis()->SetTitleFont(62);
    pull_frame6->GetYaxis()->SetTitleFont(62);
    pull_frame6->GetYaxis()->SetLabelFont(62);
    pull_frame6->GetXaxis()->SetLabelFont(62);
    pull_frame6->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame6->GetYaxis()->SetNdivisions(4);
    pull_frame6->GetXaxis()->SetNdivisions(8);
    pull_frame6->Draw();
    TLine* zeroLine6 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine6->SetLineColor(kBlack);
    zeroLine6->SetLineStyle(9);
    zeroLine6->SetLineWidth(2);
    zeroLine6->Draw();

    c2->Draw();


























    deltam_refit.setRange("range_1", 142.0, 143.5);
    deltam_refit.setRange("range_2", 143.5, 145.5);
    deltam_refit.setRange("range_3", 145.5, 149);
    deltam_refit.setRange("range_4", 149, 152.0);

    RooPlot* frame_1 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (142 < #Delta m < 143.5)"));
    RooPlot* frame_2 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (143.5 < #Delta m < 145.5)"));
    RooPlot* frame_3 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (145.5 < #Delta m < 149)"));
    RooPlot* frame_4 = D0_refit.frame(RooFit::Title("D0 Mass Fit Projection (149 < #Delta m < 152)"));


    data.plotOn(frame_1, RooFit::CutRange("range_1"), RooFit::Name("data_1"));
    total.plotOn(frame_1, RooFit::ProjectionRange("range_1"), RooFit::Name("total_fit_1"));
    RooPlot* pull_frame_1 = D0_refit.frame(RooFit::Title("Pull Distribution (142 < #Delta m < 143.5)"));
    RooHist* pullHist_1 = frame_1->pullHist("data_1", "total_fit_1");
    pull_frame_1->addPlotable(pullHist_1, "P");
    total.plotOn(frame_1, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_1"), RooFit::Name("signal_fit_1"));
    total.plotOn(frame_1, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range_1"), RooFit::Name("combinatorial_fit_1"));
    total.plotOn(frame_1, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_1"), RooFit::Name("random_slowpion_fit_1"));
    total.plotOn(frame_1, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range_1"), RooFit::Name("misreco_fit_1"));

    data.plotOn(frame_2, RooFit::CutRange("range_2"), RooFit::Name("data_2"));
    total.plotOn(frame_2, RooFit::ProjectionRange("range_2"), RooFit::Name("total_fit_2"));
    RooPlot* pull_frame_2 = D0_refit.frame(RooFit::Title("Pull Distribution (143.5 < #Delta m < 145.5)"));
    RooHist* pullHist_2 = frame_2->pullHist("data_2", "total_fit_2");
    pull_frame_2->addPlotable(pullHist_2, "P");
    total.plotOn(frame_2, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_2"), RooFit::Name("signal_fit_2"));
    total.plotOn(frame_2, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range_2"), RooFit::Name("combinatorial_fit_2"));
    total.plotOn(frame_2, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_2"), RooFit::Name("random_slowpion_fit_2"));
    total.plotOn(frame_2, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range_2"), RooFit::Name("misreco_fit_2"));

    data.plotOn(frame_3, RooFit::CutRange("range_3"), RooFit::Name("data_3"));
    total.plotOn(frame_3, RooFit::ProjectionRange("range_3"), RooFit::Name("total_fit_3"));
    RooPlot* pull_frame_3 = D0_refit.frame(RooFit::Title("Pull Distribution (145.5 < #Delta m < 149)"));
    RooHist* pullHist_3 = frame_3->pullHist("data_3", "total_fit_3");
    pull_frame_3->addPlotable(pullHist_3, "P");
    total.plotOn(frame_3, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_3"), RooFit::Name("signal_fit_3"));
    total.plotOn(frame_3, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range_3"), RooFit::Name("combinatorial_fit_3"));
    total.plotOn(frame_3, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_3"), RooFit::Name("random_slowpion_fit_3"));
    total.plotOn(frame_3, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range_3"), RooFit::Name("misreco_fit_3"));

    data.plotOn(frame_4, RooFit::CutRange("range_4"), RooFit::Name("data_4"));
    total.plotOn(frame_4, RooFit::ProjectionRange("range_4"), RooFit::Name("total_fit_4"));
    RooPlot* pull_frame_4 = D0_refit.frame(RooFit::Title("Pull Distribution (149 < #Delta m < 152)"));
    RooHist* pullHist_4 = frame_4->pullHist("data_4", "total_fit_4");
    pull_frame_4->addPlotable(pullHist_4, "P");
    total.plotOn(frame_4, RooFit::Components("signal"), RooFit::LineColor(kRed), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_4"), RooFit::Name("signal_fit_4"));
    total.plotOn(frame_4, RooFit::Components("combinatorial"), RooFit::LineColor(kGreen+2), RooFit::LineStyle(2), RooFit::LineWidth(8), RooFit::ProjectionRange("range_4"), RooFit::Name("combinatorial_fit_4"));
    total.plotOn(frame_4, RooFit::Components("random_slowpion"), RooFit::LineColor(kMagenta+3), RooFit::LineStyle(9), RooFit::LineWidth(8), RooFit::ProjectionRange("range_4"), RooFit::Name("random_slowpion_fit_4"));
    total.plotOn(frame_4, RooFit::Components("misreco_D0"), RooFit::LineColor(kOrange+7), RooFit::LineStyle(1), RooFit::LineWidth(8), RooFit::ProjectionRange("range_4"), RooFit::Name("misreco_fit_4"));

    // Create a canvas to display the plots for D0 mass projections
    TCanvas* c3 = new TCanvas("c3", "D0 Mass Projections", 1800, 600);

    // Create pads for each range
    TPad* pad_1 = new TPad("pad_1", "pad_1", 0.0, 0.0, 0.25, 1.0);
    TPad* pad_2 = new TPad("pad_2", "pad_2", 0.25, 0.0, 0.50, 1.0);
    TPad* pad_3 = new TPad("pad_3", "pad_3", 0.50, 0.0, 0.75, 1.0);
    TPad* pad_4 = new TPad("pad_4", "pad_4", 0.75, 0.0, 1.0, 1.0);

    pad_1->Draw();
    pad_2->Draw();
    pad_3->Draw();
    pad_4->Draw();

    pad_1->SetLogy();
    pad_2->SetLogy();
    pad_3->SetLogy();
    pad_4->SetLogy();
    // Plot for range4
    pad_1->cd();
    TPad* pad_1_1 = new TPad("pad_1_1", "pad_1_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad_1_2 = new TPad("pad_1_2", "pad_1_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad_1_1->SetBottomMargin(0.02);
    pad_1_2->SetTopMargin(0.02);
    pad_1_1->Draw();
    pad_1_2->Draw();

    pad_1_1->cd();
    frame_1->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame_1->GetYaxis()->SetTitleOffset(1.5); // Adjusted
    frame_1->GetYaxis()->SetTitleFont(62);
    frame_1->GetYaxis()->SetLabelFont(62);
    frame_1->GetYaxis()->SetTitleSize(0.06);
    frame_1->GetYaxis()->SetLabelSize(0.05);
    frame_1->Draw();
    TLatex* deltam_left2 = new TLatex(0.18, 0.95, "#Delta m #in [142, 143.5]MeV/c^{2}");
    deltam_left2->SetNDC();
    deltam_left2->SetTextFont(62);
    deltam_left2->SetTextSize(0.06);
    deltam_left2->SetTextColor(kBlack);
    deltam_left2->Draw();

    pad_1_2->cd();
    pull_frame_1->GetYaxis()->SetTitle("Pull");
    pull_frame_1->GetYaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_1->GetYaxis()->SetTitleOffset(0.45); // Adjusted
    pull_frame_1->GetYaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_1->GetXaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_1->GetXaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_1->GetXaxis()->SetTitleFont(62);
    pull_frame_1->GetYaxis()->SetTitleFont(62);
    pull_frame_1->GetYaxis()->SetLabelFont(62);
    pull_frame_1->GetXaxis()->SetLabelFont(62);
    pull_frame_1->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame_1->GetYaxis()->SetNdivisions(4);
    pull_frame_1->GetXaxis()->SetNdivisions(5);
    pull_frame_1->Draw();
    TLine* zeroLine_1 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine_1->SetLineColor(kBlack);
    zeroLine_1->SetLineStyle(9);
    zeroLine_1->SetLineWidth(2);
    zeroLine_1->Draw();

    // Plot for range5
    pad_2->cd();
    TPad* pad_2_1 = new TPad("pad_2_1", "pad_2_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad_2_2 = new TPad("pad_2_2", "pad_2_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad_2_1->SetBottomMargin(0.02);
    pad_2_2->SetTopMargin(0.02);
    pad_2_1->Draw();
    pad_2_2->Draw();

    pad_2_1->cd();
    frame_2->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame_2->GetYaxis()->SetTitleOffset(1.4); // Unchanged
    frame_2->GetYaxis()->SetTitleFont(62);
    frame_2->GetYaxis()->SetLabelFont(62);
    frame_2->GetYaxis()->SetTitleSize(0.06);
    frame_2->GetYaxis()->SetLabelSize(0.07);
    frame_2->Draw();
    TLatex* deltam_peak2 = new TLatex(0.17, 0.95, "#Delta m #in [143.5, 145.5]MeV/c^{2}");
    deltam_peak2->SetNDC();
    deltam_peak2->SetTextFont(62);
    deltam_peak2->SetTextSize(0.06);
    deltam_peak2->SetTextColor(kBlack);
    deltam_peak2->Draw();

    pad_2_2->cd();
    pull_frame_2->GetYaxis()->SetTitle("Pull");
    pull_frame_2->GetYaxis()->SetTitle("Pull");
    pull_frame_2->GetYaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_2->GetYaxis()->SetTitleOffset(0.45); // Adjusted
    pull_frame_2->GetYaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_2->GetXaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_2->GetXaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_2->GetXaxis()->SetTitleFont(62);
    pull_frame_2->GetYaxis()->SetTitleFont(62);
    pull_frame_2->GetYaxis()->SetLabelFont(62);
    pull_frame_2->GetXaxis()->SetLabelFont(62);
    pull_frame_2->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame_2->GetYaxis()->SetNdivisions(4);
    pull_frame_2->GetXaxis()->SetNdivisions(5);
    pull_frame_2->Draw();
    TLine* zeroLine_2 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine_2->SetLineColor(kBlack);
    zeroLine_2->SetLineStyle(9);
    zeroLine_2->SetLineWidth(2);
    zeroLine_2->Draw();

    // Plot for range6
    pad_3->cd();
    TPad* pad_3_1 = new TPad("pad_3_1", "pad_3_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad_3_2 = new TPad("pad_3_2", "pad_3_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad_3_1->SetBottomMargin(0.02);
    pad_3_2->SetTopMargin(0.02);
    pad_3_1->Draw();
    pad_3_2->Draw();

    pad_3_1->cd();
    frame_3->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame_3->GetYaxis()->SetTitleOffset(1.5); // Adjusted
    frame_3->GetYaxis()->SetTitleFont(62);
    frame_3->GetYaxis()->SetLabelFont(62);
    frame_3->GetYaxis()->SetTitleSize(0.06);
    frame_3->GetYaxis()->SetLabelSize(0.05);
    frame_3->Draw();
    TLatex* deltam_right1 = new TLatex(0.17, 0.95, "#Delta m #in [145.5, 149]MeV/c^{2}");
    deltam_right1->SetNDC();
    deltam_right1->SetTextFont(62);
    deltam_right1->SetTextSize(0.06);
    deltam_right1->SetTextColor(kBlack);
    deltam_right1->Draw();

    pad_3_2->cd();
    pull_frame_3->GetYaxis()->SetTitle("Pull");
    pull_frame_3->GetYaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_3->GetYaxis()->SetTitleOffset(0.45); // Adjusted
    pull_frame_3->GetYaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_3->GetXaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_3->GetXaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_3->GetXaxis()->SetTitleFont(62);
    pull_frame_3->GetYaxis()->SetTitleFont(62);
    pull_frame_3->GetYaxis()->SetLabelFont(62);
    pull_frame_3->GetXaxis()->SetLabelFont(62);
    pull_frame_3->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame_3->GetYaxis()->SetNdivisions(4);
    pull_frame_3->GetXaxis()->SetNdivisions(5);
    pull_frame_3->Draw();
    TLine* zeroLine_3 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine_3->SetLineColor(kBlack);
    zeroLine_3->SetLineStyle(9);
    zeroLine_3->SetLineWidth(2);
    zeroLine_3->Draw();

    pad_4->cd();
    TPad* pad_4_1 = new TPad("pad_4_1", "pad_4_1", 0.0, 0.25, 1.0, 1.0); // Adjusted
    TPad* pad_4_2 = new TPad("pad_4_2", "pad_4_2", 0.0, 0.0, 1.0, 0.25); // Adjusted
    pad_4_1->SetBottomMargin(0.02);
    pad_4_2->SetTopMargin(0.02);
    pad_4_1->Draw();
    pad_4_2->Draw();

    pad_4_1->cd();
    frame_4->GetYaxis()->SetTitle("Events / (0.9 MeV/c^{2})");
    frame_4->GetYaxis()->SetTitleOffset(1.5); // Adjusted
    frame_4->GetYaxis()->SetTitleFont(62);
    frame_4->GetYaxis()->SetLabelFont(62);
    frame_4->GetYaxis()->SetTitleSize(0.06);
    frame_4->GetYaxis()->SetLabelSize(0.05);
    frame_4->Draw();
    TLatex* deltam_right2 = new TLatex(0.17, 0.95, "#Delta m #in [149, 152]MeV/c^{2}");
    deltam_right2->SetNDC();
    deltam_right2->SetTextFont(62);
    deltam_right2->SetTextSize(0.06);
    deltam_right2->SetTextColor(kBlack);
    deltam_right2->Draw();

    pad_4_2->cd();
    pull_frame_4->GetYaxis()->SetTitle("Pull");
    pull_frame_4->GetYaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_4->GetYaxis()->SetTitleOffset(0.45); // Adjusted
    pull_frame_4->GetYaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_4->GetXaxis()->SetLabelSize(0.14); // Adjusted
    pull_frame_4->GetXaxis()->SetTitleSize(0.14); // Adjusted
    pull_frame_4->GetXaxis()->SetTitleFont(62);
    pull_frame_4->GetYaxis()->SetTitleFont(62);
    pull_frame_4->GetYaxis()->SetLabelFont(62);
    pull_frame_4->GetXaxis()->SetLabelFont(62);
    pull_frame_4->GetYaxis()->SetRangeUser(-12, 12);
    pull_frame_4->GetYaxis()->SetNdivisions(4);
    pull_frame_4->GetXaxis()->SetNdivisions(5);
    pull_frame_4->Draw();
    TLine* zeroLine_4 = new TLine(1820, 0.0, 1910, 0.0);
    zeroLine_4->SetLineColor(kBlack);
    zeroLine_4->SetLineStyle(9);
    zeroLine_4->SetLineWidth(2);
    zeroLine_4->Draw();

    c3->Draw();
}
