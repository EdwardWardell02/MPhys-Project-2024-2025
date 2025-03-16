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
