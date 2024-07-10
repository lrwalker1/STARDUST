#include <iostream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TStyle.h>
#include <cmath>

void sortEvents(const char* inputFileName, const char* outputFileName) {
    // Toggle to switch between e and edet branches
    bool useEdet = true; // Set to false to use the e branch

    // Open the input ROOT file
    TFile* inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: could not open input file " << inputFileName << std::endl;
        return;
    }

    // Get the input tree
    TTree* inputTree = (TTree*)inputFile->Get("h");
    if (!inputTree) {
        std::cerr << "Error: could not find tree in input file" << std::endl;
        inputFile->Close();
        return;
    }

    // Set branch addresses
    int n;
    int evt;
    double x[1000], y[1000], z[1000], e[1000], t[1000], edet[1000];
    int pdg[1000], vlm[1000];

    inputTree->SetBranchAddress("n", &n);
    inputTree->SetBranchAddress("evt", &evt);
    inputTree->SetBranchAddress("x", x);
    inputTree->SetBranchAddress("y", y);
    inputTree->SetBranchAddress("z", z);
    inputTree->SetBranchAddress("e", e);
    inputTree->SetBranchAddress("t", t);
    inputTree->SetBranchAddress("pdg", pdg);
    inputTree->SetBranchAddress("vlm", vlm);

    if (useEdet) {
        inputTree->SetBranchAddress("edet", edet);
    }

    // Open the output ROOT file
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: could not open output file " << outputFileName << std::endl;
        return;
    }

    int nBins = 500;
    int hpgeBins = 16384;
    double maxEnergy = 400.0;
    double hpgeMaxE = 3000;
    double dt = 100;

    // Create a new 2D histogram for coincidences
    TH2D* coincidenceHist = new TH2D("gamma_beta", "Coincidences;#beta Energy (keV);#gamma Energy (keV)", nBins, 0, maxEnergy, nBins, 0, maxEnergy);

    // Create a new 1D histogram for energy distribution in detectors 1-21
    TH1D* energyHist = new TH1D("energy_gamma", "HPGe Edep;#gamma Energy (keV);Counts", (400/hpgeMaxE)*hpgeBins, 0, (400/hpgeMaxE)*hpgeMaxE);

    // Create a new 2D histogram for gamma-gamma coincidences
    TH2D* gammaGammaHist = new TH2D("gamma_gamma", "Gamma-Gamma Coincidences;#gamma_{1} Energy (keV);#gamma_{2} Energy (keV)", (400/hpgeMaxE)*hpgeBins, 0, (400/hpgeMaxE)*hpgeMaxE, (400/hpgeMaxE)*hpgeBins, 0, (400/hpgeMaxE)*hpgeMaxE);

    // Create a new 2D histogram for zoomed-in gamma-gamma coincidences
    TH2D* gammaGammaZoomHist = new TH2D("gamma_gamma_zoom", "Zoomed Gamma-Gamma Coincidences;#gamma_{1} Energy (keV);#gamma_{2} Energy (keV)", 100, 79.1, 82, 100, 79.1, 82);

    // Loop over the entries in the tree
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        inputTree->GetEntry(i);

        // Find coincidences within the event
        for (int j = 0; j < n; ++j) {
            double energy = useEdet ? edet[j] : e[j];

            // Fill the energy histogram for detectors 1-21
            if (vlm[j] >= 1 && vlm[j] <= 21) {
                energyHist->Fill(energy);
            }

            for (int k = j + 1; k < n; ++k) {
                double energy_k = useEdet ? edet[k] : e[k];
                if (std::abs(t[k] - t[j]) <= dt) {
                    // Check for detectors 801 and 802 vs detectors 1-21
                    if ((vlm[j] == 801 || vlm[j] == 802) && (vlm[k] >= 1 && vlm[k] <= 21)) {
                        coincidenceHist->Fill(energy, energy_k);
                    } else if ((vlm[k] == 801 || vlm[k] == 802) && (vlm[j] >= 1 && vlm[j] <= 21)) {
                        coincidenceHist->Fill(energy_k, energy);
                    }
                }
            }

            // Find gamma-gamma coincidences within the event
            if (vlm[j] >= 1 && vlm[j] <= 21) {
                for (int k = 0; k < n; ++k) {
                    if (k != j && vlm[k] >= 1 && vlm[k] <= 21 && std::abs(t[k] - t[j]) <= dt) {
                        double energy_k = useEdet ? edet[k] : e[k];
                        gammaGammaHist->Fill(energy, energy_k);
                        gammaGammaHist->Fill(energy_k, energy); // to make it symmetric

                        // Fill the zoomed-in gamma-gamma histogram
                        if (energy >= 75 && energy <= 85 && energy_k >= 75 && energy_k <= 85) {
                            gammaGammaZoomHist->Fill(energy, energy_k);
                            gammaGammaZoomHist->Fill(energy_k, energy); // to make it symmetric
                        }
                    }
                }
            }
        }
    }

    // Write the histograms to the output file
    outputFile->cd();
    coincidenceHist->Write();
    energyHist->Write();
    gammaGammaHist->Write();
    gammaGammaZoomHist->Write();

    // Set a predefined color palette
    gStyle->SetPalette(kCool);

    // Draw the histograms and save as images
    TCanvas* c1 = new TCanvas("c1", "Coincidences", 800, 600);
    coincidenceHist->SetStats(0); // Disable the statistics box
    coincidenceHist->SetTitle(""); // Remove the title
    coincidenceHist->GetYaxis()->SetRangeUser(0, 250); // Set Y-axis range to 250 keV
    coincidenceHist->Draw("COLZ");
    c1->SetLogz();
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15); // Set tight margins
    c1->SaveAs("gamma_beta_coincidences.png");

    TCanvas* c2 = new TCanvas("c2", "Energy Distribution", 800, 600);
    energyHist->SetStats(0);  // Disable the statistics box
    energyHist->SetTitle(""); // Remove the title
    energyHist->GetXaxis()->SetRangeUser(0, 85); // Set X-axis range to 85 keV
    energyHist->Draw();
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15); // Set tight margins
    c2->SaveAs("gamma_edep.png");

    TCanvas* c3 = new TCanvas("c3", "Gamma-Gamma Coincidences", 800, 600);
    gammaGammaHist->SetStats(0); // Disable the statistics box
    gammaGammaHist->SetTitle(""); // Remove the title
    gammaGammaHist->GetYaxis()->SetRangeUser(0, 250); // Set Y-axis range to 250 keV
    gammaGammaHist->GetXaxis()->SetRangeUser(0, 250); // Set X-axis range to 250 keV
    gammaGammaHist->Draw("COLZ"); // Draw as heatmap
    c3->SetLogz();
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15); // Set tight margins
    c3->SaveAs("gamma_gamma_coincidences.png");

    TCanvas* c4 = new TCanvas("c4", "Zoomed Gamma-Gamma Coincidences", 800, 600);
    gammaGammaZoomHist->SetStats(0); // Disable the statistics box
    gammaGammaZoomHist->SetTitle(""); // Remove the title
    gammaGammaZoomHist->GetXaxis()->SetTitle(""); // Remove X-axis title
    gammaGammaZoomHist->GetYaxis()->SetTitle(""); // Remove Y-axis title
    gammaGammaZoomHist->GetXaxis()->SetLabelSize(0.05); // Set X-axis label size
    gammaGammaZoomHist->GetYaxis()->SetLabelSize(0.05); // Set Y-axis label size
    gammaGammaZoomHist->Draw("COLZ"); // Draw as heatmap
    c4->SetLogz();
    gPad->SetMargin(0.15, 0.15, 0.15, 0.15); // Set tight margins
    c4->SaveAs("gamma_gamma_zoom.png");

    // Close the output file
    outputFile->Close();
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFile> <outputFile>" << std::endl;
        return 1;
    }

    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];

    sortEvents(inputFileName, outputFileName);

    return 0;
}

