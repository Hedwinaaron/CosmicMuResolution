#include <sys/stat.h>
#include <iostream>
#include <filesystem>
#include <TROOT.h>
#include <TH1F.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <string>
#include <cmath>
#include <vector>
#include <array>

namespace fs = std::filesystem;

// Global constants
const std::array<double, 9> PT_BINS = {20, 30, 40, 50, 70, 100, 150, 200, 1000};
const float binSize[8] = {5.0, 5.0, 5.0, 10.0, 15.0, 25.0, 25.0, 400.0};
const float binCenter[8] = {25.0, 35.0, 45.0, 60.0, 85.0, 125.0, 175.0, 600.0};
// Specify data type
std::string datatype = "DATA";

// Function declarations
void InitializeHistograms(std::vector<TH1F*>& ptHistograms, TH1F*& histPtTag, TH1F*& histEtaTag, TH1F*& histPhiTag, TH1F*& histChargeTag,
                          TH1F*& histPtProbe, TH1F*& histEtaProbe, TH1F*& histPhiProbe, TH1F*& histChargeProbe, TH1F*& histPtDGL);
void ProcessFiles(const std::vector<std::string>& files, std::vector<TH1F*>& ptHistograms, TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
                  TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL);
void FitAndDrawHistograms(const std::vector<TH1F*>& ptHistograms, TFile* outputFile);
void DrawControlPlots(TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
                      TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL, TFile* outputFile);

void DGM_sel() {


    // Create a ROOT file to save histograms
    TFile* outputFile = new TFile(("Cosmics_muons_" + datatype + ".root").c_str(), "RECREATE");

    // List of directories containing ROOT files
    std::vector<std::string> directories;
    if (datatype == "DATA") {
        directories = {
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023C_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171725/0000",
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023D_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171752/0000",
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023E_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171805/0000",
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/NoBPTX/CosmicsAnalysis_Run2023F_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test7/241209_171822/0000"
        };
    } else if (datatype == "MC") {
        directories = {
            "/eos/user/h/hencinas/Mu_efficiency_Analysis/Cosmics_Ntuples_Data-2023/UndergroundCosmiHPLooseMu_bottomPhiFilter/CosmicsAnalysis_Run2022_MC_MiniAOD-Ntuples_TnP_CMSSW_13_0_13_test10/241213_225101/0000"
        };
    } else {
        std::cerr << "Invalid datatype: " << datatype << std::endl;
        return;
    }

    // Vector to store all ROOT files from all directories
    std::vector<std::string> files;

    // Loop over directories to collect ROOT files
    for (const auto& dir : directories) {
        if (fs::is_directory(dir)) {
            for (const auto& entry : fs::directory_iterator(dir)) {
                if (entry.path().extension() == ".root") {
                    files.push_back(entry.path().string());
                    std::cout << "Found file: " << entry.path() << std::endl;
                }
            }
        } else {
            std::cerr << "Invalid directory: " << dir << std::endl;
        }
    }

    // Histograms for pT bins
    std::vector<TH1F*> ptHistograms;
    TH1F* histPtTag, *histEtaTag, *histPhiTag, *histChargeTag;
    TH1F* histPtProbe, *histEtaProbe, *histPhiProbe, *histChargeProbe;
    TH1F* histPtDGL;

    // Initialize histograms
    InitializeHistograms(ptHistograms, histPtTag, histEtaTag, histPhiTag, histChargeTag,
                         histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL);

    // Process files
    ProcessFiles(files, ptHistograms, histPtTag, histEtaTag, histPhiTag, histChargeTag,
                 histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL);

    // Fit and draw histograms
    FitAndDrawHistograms(ptHistograms, outputFile);

    // Draw control plots
    DrawControlPlots(histPtTag, histEtaTag, histPhiTag, histChargeTag,
                     histPtProbe, histEtaProbe, histPhiProbe, histChargeProbe, histPtDGL, outputFile);

    // Clean up
    for (auto hist : ptHistograms) delete hist;
    delete histPtTag;
    delete histPtProbe;

    // Close the output file
    outputFile->Close();
    delete outputFile;

    std::cout << "Finished processing and histogram saving." << std::endl;
}

void InitializeHistograms(std::vector<TH1F*>& ptHistograms, TH1F*& histPtTag, TH1F*& histEtaTag, TH1F*& histPhiTag, TH1F*& histChargeTag,
                          TH1F*& histPtProbe, TH1F*& histEtaProbe, TH1F*& histPhiProbe, TH1F*& histChargeProbe, TH1F*& histPtDGL) {
    // Initialize histograms for pT each bin
    for (size_t i = 0; i < PT_BINS.size() - 1; ++i) {
        ptHistograms.push_back(new TH1F(Form("pt_%g-%g", PT_BINS[i], PT_BINS[i + 1]),
                                        Form("pT %g-%g", PT_BINS[i], PT_BINS[i + 1]), 100, -0.3, 0.3));
    }

    // Initialize other histograms
    histPtTag = new TH1F("hist_pt_tag", "Tag muon pT", 100, -10, 1000);
    histEtaTag = new TH1F("hist_eta_tag", "Tag muon #eta", 60, -3, 3);
    histPhiTag = new TH1F("hist_phi_tag", "Tag muon #phi", 60, -3, -3);
    histChargeTag = new TH1F("hist_charge_tag", "Tag muon q", 10, -4, 4);

    histPtProbe = new TH1F("hist_pt_probe", "Probe muon pT", 100, -10, 1000);
    histEtaProbe = new TH1F("hist_eta_probe", "Probe muon #eta", 60, -3, 3);
    histPhiProbe = new TH1F("hist_phi_probe", "Probe muon #phi", 60, -3, -3);
    histChargeProbe = new TH1F("hist_charge_probe", "Probe muon q", 10, -4, 4);

    histPtDGL = new TH1F("hist_pt_dGL", "DGL muon pT", 40, 0, 33);
}

void ProcessFiles(const std::vector<std::string>& files, std::vector<TH1F*>& ptHistograms, TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
                  TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL) {
    // Loop through each file
    for (const auto& file : files) {
        std::cout << "Reading file: " << file << std::endl;
        TFile* f = TFile::Open(file.c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << file << std::endl;
            continue;
        }

        // Get the TTree
        TTree* tree = dynamic_cast<TTree*>(f->Get("Events"));
        if (!tree) {
            std::cerr << "Error: TTree 'Events' not found!" << std::endl;
            f->Close();
            delete f;
            continue;
        }

        // Variables to hold branch data
        const int maxMuons = 200; // Adjust this to the maximum expected size of ndmu
        int ndmu;
        int event;
        int dmuNumberOfChambersCSCorDT[maxMuons];
        float dmuDglPt[maxMuons];
        float dmuDglEta[maxMuons];
        float dmuDglPhi[maxMuons];
        float dmuDglDz[maxMuons];
        float dmuDglDxy[maxMuons];
        bool dmuDglPassTagID[maxMuons];
        bool dmuDglHasProbe[maxMuons];
        int dmuDglProbeID[maxMuons];
        int dmuIsDGL[maxMuons];
        float dmuDglCharge[maxMuons];
        bool hltL2Mu10NoVertexNoBPTX3BX = false;

        // Link branches
        tree->SetBranchAddress("ndmu", &ndmu);
        tree->SetBranchAddress("event", &event);
        tree->SetBranchAddress("dmu_numberOfChambersCSCorDT", dmuNumberOfChambersCSCorDT);
        tree->SetBranchAddress("dmu_dgl_pt", dmuDglPt);
        tree->SetBranchAddress("dmu_dgl_eta", dmuDglEta);
        tree->SetBranchAddress("dmu_dgl_phi", dmuDglPhi);
        tree->SetBranchAddress("dmu_dgl_dz", dmuDglDz);
        tree->SetBranchAddress("dmu_dgl_dxy", dmuDglDxy);
        tree->SetBranchAddress("dmu_dgl_charge", dmuDglCharge);
        tree->SetBranchAddress("dmu_dgl_passTagID", dmuDglPassTagID);
        tree->SetBranchAddress("dmu_dgl_probeID", dmuDglProbeID);
        tree->SetBranchAddress("dmu_dgl_hasProbe", dmuDglHasProbe);
        tree->SetBranchAddress("dmu_isDGL", dmuIsDGL);
        tree->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &hltL2Mu10NoVertexNoBPTX3BX);

        // Loop over all entries
        int nEntries = tree->GetEntries();
        for (int entry = 0; entry < nEntries; entry++) {
            tree->GetEntry(entry);

            // Only select events that pass the HLT path
            if(datatype == "MC") hltL2Mu10NoVertexNoBPTX3BX=true;
            if (!hltL2Mu10NoVertexNoBPTX3BX) continue;

            // Tag variables
            float tagPt = 0.0;
            float tagEta = 0.0;
            float tagPhi = 0.0;
            float tagDz = 0.0;
            float tagDxy = 0.0;
            float tagCharge = 0.0;

            // Probe variables
            float probePt = 0.0;
            float probeEta = 0.0;
            float probePhi = 0.0;
            float probeDz = 0.0;
            float probeDxy = 0.0;
            float probeCharge = 0.0;

            // Select events with 2 muons or more
            if (ndmu < 2) continue;

            // Loop over the number of muons (ndmu) and select the Tag muon
            for (int i = 0; i < ndmu; i++) {
                histPtDGL->Fill(dmuDglPt[i]);

                if (!dmuDglPassTagID[i]) continue;

                if (dmuDglPassTagID[i] && dmuIsDGL[i] == 1) {
                    // Fill tag histograms
                    tagPt = dmuDglPt[i];
                    tagEta = dmuDglEta[i];
                    tagPhi = dmuDglPhi[i];
                    tagDz = dmuDglDz[i];
                    tagDxy = dmuDglDxy[i];
                    tagCharge = dmuDglCharge[i];
                    histPtTag->Fill(dmuDglPt[i]);
                    histEtaTag->Fill(dmuDglEta[i]);
                    histPhiTag->Fill(dmuDglPhi[i]);
                    histChargeTag->Fill(dmuDglCharge[i]);

                    // Probe selection
                    if (dmuDglHasProbe[i]) {
                        int probeIndex = dmuDglProbeID[i];
                        probePt = dmuDglPt[probeIndex];
                        probeEta = dmuDglEta[probeIndex];
                        probePhi = dmuDglPhi[probeIndex];
                        probeDz = dmuDglDz[probeIndex];
                        probeDxy = dmuDglDxy[probeIndex];
                        probeCharge = dmuDglCharge[probeIndex];
                        histPtProbe->Fill(dmuDglPt[probeIndex]);
                        histEtaProbe->Fill(dmuDglEta[probeIndex]);
                        histPhiProbe->Fill(dmuDglPhi[probeIndex]);
                        histChargeProbe->Fill(dmuDglCharge[probeIndex]);

                        // Calculate the resolution
                        float invUppt = std::abs(dmuDglCharge[probeIndex] / dmuDglPt[probeIndex]);
                        float invDownpt = std::abs(dmuDglCharge[i] / dmuDglPt[i]);
                        double resolution = (invUppt - invDownpt) / (std::sqrt(2) * invDownpt);

                        // Fill pT histograms
                        for (size_t j = 0; j < PT_BINS.size() - 1; ++j) {
                            if (PT_BINS[j] < tagPt && tagPt < PT_BINS[j + 1]) {
                                ptHistograms[j]->Fill(resolution);
                            }
                        }
                    }
                }
            }
        }

        // Close the file
        f->Close();
        delete f;
    }
}

void FitAndDrawHistograms(const std::vector<TH1F*>& ptHistograms, TFile* outputFile) {
    // Mean and Sigma arrays
    float mean[9] = {};
    float meanErr[9] = {};
    float sigma[9] = {};
    float sigmaErr[9] = {};

    outputFile->cd();

    // Draw the histograms for pT
    for (size_t i = 0; i < ptHistograms.size(); ++i) {
        TCanvas* canvas = new TCanvas(Form("c_pt_%zu", i), Form("pT Histogram %g-%g", PT_BINS[i], PT_BINS[i + 1]), 800, 600);
        canvas->SetGrid();

        // Fit the histogram with a Gaussian
        TF1* gaussFit = new TF1("gaussFit", "gaus", -0.3, 0.3);
        ptHistograms[i]->Fit(gaussFit, "S");

        // Extract fit parameters
        double chi2 = gaussFit->GetChisquare();
        int ndf = gaussFit->GetNDF();
        double chi2PerNdf = chi2 / ndf;
        double meanVal = gaussFit->GetParameter(1);
        double stdDev = gaussFit->GetParameter(2);

        mean[i] = meanVal;
        meanErr[i] = gaussFit->GetParError(1);
        sigma[i] = stdDev;
        sigmaErr[i] = gaussFit->GetParError(2);

        // Print parameters
        std::cout << "Histogram: pT " << PT_BINS[i] << "-" << PT_BINS[i + 1]
                  << ", Mean: " << meanVal << ", Std Dev: " << stdDev << std::endl;

        // Draw histogram and fit
        ptHistograms[i]->Draw();
        ptHistograms[i]->GetYaxis()->SetTitle("Events");
        ptHistograms[i]->GetXaxis()->SetTitle("q/p_{T} Residuals");
        gaussFit->Draw("SAME");

        // Add legend with parameters
        TLegend* legend = new TLegend(0.6, 0.5, 0.9, 0.75);
        legend->SetTextSize(0.04);
        legend->SetFillStyle(0);
        legend->AddEntry((TObject*)nullptr, Form("Mean = %.4f", meanVal), "");
        legend->AddEntry((TObject*)nullptr, Form("#sigma = %.4f", stdDev), "");
        legend->AddEntry((TObject*)nullptr, Form("#chi^2 = %.4f", chi2), "");
        legend->AddEntry((TObject*)nullptr, Form("NDF = %d", ndf), "");
        legend->AddEntry((TObject*)nullptr, Form("#chi^2/ndf = %.4f", chi2PerNdf), "");
        legend->SetBorderSize(0);
        legend->Draw();

        ptHistograms[i]->Write();

        // Save canvas
        canvas->SaveAs(Form("pT_Plots/dmu_dgl_pt_%g-%g.png", PT_BINS[i], PT_BINS[i + 1]));
        delete gaussFit;
        delete canvas;
    }

    // Draw mean and sigma plots
    TCanvas* cc1 = new TCanvas("cc1", "cc1", 1000, 750);
    cc1->cd();
    cc1->SetGrid();
    cc1->SetLogx();
    cc1->SetLeftMargin(0.24);

    auto meanGraph = new TGraphErrors(8, binCenter, mean, binSize, meanErr);
    meanGraph->SetTitle("mean_total");
    meanGraph->SetLineWidth(2);
    meanGraph->SetMarkerStyle(8);
    meanGraph->GetXaxis()->SetTitle("p_{T} #mu_{ref} [GeV] ");
    meanGraph->GetYaxis()->SetTitle("Mean of q/p_{T} relative residual ");
    meanGraph->Draw("A*");
    meanGraph->Write("mean_total");
    cc1->SaveAs("mean_total_.png");

    TCanvas* cc2 = new TCanvas("cc2", "cc2", 1000, 750);
    cc2->cd();
    cc2->SetGrid();
    cc2->SetLogx();
    cc2->SetLeftMargin(0.24);

    auto sigmaGraph = new TGraphErrors(8, binCenter, sigma, binSize, sigmaErr);
    sigmaGraph->SetTitle("Sigma_total");
    sigmaGraph->SetLineWidth(2);
    sigmaGraph->SetMarkerStyle(8);
    sigmaGraph->GetXaxis()->SetTitle("p_{T} #mu_{ref} [GeV] ");
    sigmaGraph->GetYaxis()->SetTitle("#sigma of q/p_{T} relative residual ");
    sigmaGraph->Draw("A*");
    sigmaGraph->Write("Sigma_total");
    cc2->SaveAs("Sigma_total_.png");
}

void DrawControlPlots(TH1F* histPtTag, TH1F* histEtaTag, TH1F* histPhiTag, TH1F* histChargeTag,
                      TH1F* histPtProbe, TH1F* histEtaProbe, TH1F* histPhiProbe, TH1F* histChargeProbe, TH1F* histPtDGL, TFile* outputFile) {
    // Draw control plots for Tag muons
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
    histPtTag->Draw();
    histPtTag->GetYaxis()->SetTitle("Events");
    histPtTag->GetXaxis()->SetTitle("p_{T}");
    histPtTag->Write();
    c1->SaveAs("Control_plots/hist_pt_tag.png");

    TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
    histEtaTag->Draw();
    histEtaTag->GetYaxis()->SetTitle("Events");
    histEtaTag->GetXaxis()->SetTitle("#eta");
    histEtaTag->Write();
    c2->SaveAs("Control_plots/hist_eta_tag.png");

    TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
    histPhiTag->Draw();
    histPhiTag->GetYaxis()->SetTitle("Events");
    histPhiTag->GetXaxis()->SetTitle("#phi");
    histPhiTag->Write();
    c3->SaveAs("Control_plots/hist_phi_tag.png");

    TCanvas* c63 = new TCanvas("c63", "c63", 800, 600);
    histChargeTag->Draw();
    histChargeTag->GetYaxis()->SetTitle("Events");
    histChargeTag->GetXaxis()->SetTitle("charge");
    histChargeTag->Write();
    c63->SaveAs("Control_plots/hist_charge_tag.png");

    // Draw control plots for Probe muons
    TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
    histPtProbe->Draw();
    histPtProbe->GetYaxis()->SetTitle("Events");
    histPtProbe->GetXaxis()->SetTitle("p_{T}");
    histPtProbe->Write();
    c4->SaveAs("Control_plots/hist_pt_probe.png");

    TCanvas* c5 = new TCanvas("c5", "c5", 800, 600);
    histEtaProbe->Draw();
    histEtaProbe->GetYaxis()->SetTitle("Events");
    histEtaProbe->GetXaxis()->SetTitle("#eta");
    histEtaProbe->Write();
    c5->SaveAs("Control_plots/hist_eta_probe.png");

    TCanvas* c6 = new TCanvas("c6", "c6", 800, 600);
    histPhiProbe->Draw();
    histPhiProbe->GetYaxis()->SetTitle("Events");
    histPhiProbe->GetXaxis()->SetTitle("#phi");
    histPhiProbe->Write();
    c6->SaveAs("Control_plots/hist_phi_probe.png");

    TCanvas* c66 = new TCanvas("c66", "c66", 800, 600);
    histChargeProbe->Draw();
    histChargeProbe->GetYaxis()->SetTitle("Events");
    histChargeProbe->GetXaxis()->SetTitle("charge");
    histChargeProbe->Write();
    c66->SaveAs("Control_plots/hist_charge_probe.png");

    TCanvas* c7 = new TCanvas("c7", "c7", 800, 600);
    histPtDGL->Draw();
    c7->SaveAs("Control_plots/hist_pT_DGL.png");
}
