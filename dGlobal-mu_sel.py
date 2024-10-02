#!/usr/bin/python3
import os
import ROOT
import numpy as np

# Define constants
PATH = ""
is_MC=True
is_Data=False
data_type=""
if is_MC:
  data_type="MC"
  PATH = "/eos/user/h/hencinas/AOD/TKCosmic_38T_p20-3000/crab_TnP_ntuplizer_muon_Cosmics_Run2018_AOD_TKCosmic/220606_142613/0000/"
  print("is MC sample")
else:
  data_type="Data"
  PATH = "/eos/user/h/hencinas/AOD/NoBPTX"
  print("is Data sample")
  

PT_BINS = np.array([30, 40, 50, 70, 100, 150, 200, 350, 500, 1000])
DXY_BINS = np.array([0, 5, 10, 15, 20, 30, 40, 50, 80])
DZ_BINS = np.array([0, 5, 10, 20, 30, 45, 60, 80, 150])

# Function to create a histogram
def create_histograms():
    hist_params = [
        ("Down_dGl", "Down dGl", 1000, 0.0, 1000.0),
        ("Up_dGl", "Up dGl", 1000, 0.0, 1000.0),
        ("Resolution_hist", "Resolution hist", 200, -1, 1),
        ("Up_|dz|", "Up |dz|", 400, 0, 400),
        ("Up_|dxy|", "Up |dxy|", 400, 0, 120),
        ("Down_|dz|", "Down |dz|", 400, 0, 400),
        ("Down_|dxy|", "Down |dxy|", 400, 0, 120),
    ]
    
    # Create additional histograms for p_T, dz, and dxy ranges
    pt_histograms = [ROOT.TH1D(f"pt_{PT_BINS[i]}-{PT_BINS[i+1]}_hist", f"p_T {PT_BINS[i]}-{PT_BINS[i+1]}", 100, -0.3, 0.3) for i in range(len(PT_BINS) - 1)]
    dz_histograms = [ROOT.TH1D(f"|dz| {DZ_BINS[i]}-{DZ_BINS[i+1]}", f"|dz| {DZ_BINS[i]}-{DZ_BINS[i+1]}", 100, -0.2, 0.2) for i in range(len(DZ_BINS) - 1)]
    dxy_histograms = [ROOT.TH1D(f"|dxy| {DXY_BINS[i]}-{DXY_BINS[i+1]}", f"|dxy| {DXY_BINS[i]}-{DXY_BINS[i+1]}", 100, -0.2, 0.2) for i in range(len(DXY_BINS) - 1)]
    
    hist_list = [ROOT.TH1D(*params) for params in hist_params] + pt_histograms + dz_histograms + dxy_histograms
    
    return hist_list

# Function to save a histogram with a Gaussian fit
def save_histogram_with_fit(hist, canvas_name, filename):
    canvas = ROOT.TCanvas(canvas_name, "hist", 1000, 800)
    hist.GetXaxis().SetTitle('q/p_{T} relative residual')
    fit_result = hist.Fit("gaus", "S")
    
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.DrawText(0.7, 0.5, f"Fit mean = {round(fit_result.Parameter(1), 4)}")
    latex.DrawText(0.7, 0.45, f"Fit Std Dev = {round(fit_result.Parameter(2), 4)}")
    
    canvas.SaveAs(f"plots/histograms/{filename}.png")
    canvas.Close()

# Function to process events and fill histograms
def process_event(tree, hist_list, i):
    tree.GetEntry(i)
    tag_eta = abs(tree.tag_eta)
    
    if tree.HLT_L2Mu10_NoVertex_NoBPTX3BX_v == 1 and tree.tag_pt > 20 and tree.tag_isdGlobal == 1 and -2.6 < tree.tag_phi < -0.6 and tag_eta < 0.9 and tree.probe_isUpdGlobal:
        pair_dphi = abs(tree.tag_phi - tree.probe_updgl_phi)
        pair_dphi_reduced = abs(pair_dphi - np.pi * (1 + (pair_dphi - np.pi) / abs(pair_dphi - np.pi)))
        
        entry_down = i
        downdgl_pt = tree.tag_pt
        down_charge = tree.tag_charge
        down_absdz = abs(tree.tag_dz)
        down_absdxy = abs(tree.tag_dxy)
        
        if tree.probe_updgl_pt > 30:
            entry_up = i
            updgl_pt = tree.probe_updgl_pt
            up_charge = tree.tag_charge
            up_absdz = abs(tree.probe_updgl_dz)
            up_absdxy = abs(tree.probe_updgl_dxy)
            
            hist_list[0].Fill(downdgl_pt)   # h1
            hist_list[1].Fill(updgl_pt)     # h2
            inv_Uppt = up_charge / updgl_pt
            inv_Downpt = down_charge / downdgl_pt
            res = ((inv_Uppt) - (inv_Downpt)) / (np.sqrt(2) * (inv_Downpt))
            hist_list[2].Fill(res)  # h3
            hist_list[3].Fill(up_absdz) # h4
            hist_list[4].Fill(up_absdxy) # h5
            hist_list[5].Fill(down_absdz) # h6
            hist_list[6].Fill(down_absdxy) # h7

            # Fill histograms based on pt, dz, and dxy bins
            for j in range(len(PT_BINS) - 1):
                if PT_BINS[j] < updgl_pt < PT_BINS[j + 1]:
                    hist_list[j + 7].Fill(res)   # Corresponding p_T histogram
            
            for k in range(len(DZ_BINS) - 1):
                if DZ_BINS[k] < up_absdz < DZ_BINS[k + 1]:
                    hist_list[k + 16].Fill(res)  # Corresponding dz histogram
            
            for l in range(len(DXY_BINS) - 1):
                if DXY_BINS[l] < up_absdxy < DXY_BINS[l + 1]:
                    hist_list[l + 24].Fill(res)  # Corresponding dxy histogram
                    
                    
                    


def create_graph(title, x_title, y_title, n_data, x_data, y_data, x_errors, y_errors, canvas_name, file_name):
    graph = ROOT.TGraphErrors(n_data, x_data, y_data, x_errors, y_errors)
    graph.SetTitle(title)
    graph.GetXaxis().SetTitle(x_title)
    graph.GetYaxis().SetTitle(y_title)
    graph.Write(title)
    
    canvas = ROOT.TCanvas(canvas_name, "histograms", 1400, 1400)
    canvas.cd()
    canvas.SetGrid()
    canvas.SetLogx()
    graph.Draw('ap')
    canvas.SaveAs(f"{file_name}.png")



def fit_histogram(hist_name, index, average_resolution, sigma, res_error, sigma_error):
    hist = hist_name
    fit_result = hist.Fit("gaus", "S")
    mean = fit_result.Parameter(1)
    sigma_val = fit_result.Parameter(2)
    
    # Store the fit results
    average_resolution[index] = mean
    sigma[index] = sigma_val
    res_error[index] = fit_result.ParError(1)
    sigma_error[index] = fit_result.ParError(2)
    
    return mean, sigma_val



def main():

    if is_MC:
        files = os.listdir(PATH)
        hist_list = create_histograms()
        
        for file in files:
            inFileName = os.path.join(PATH, file)
            with ROOT.TFile.Open(inFileName, "READ") as inFile:
                tree = inFile.Get("muon/Events")
                print(f"Reading from {inFileName}")
                
                entries = tree.GetEntries()
                for i in range(entries):
                    process_event(tree, hist_list,i)
    else:
        directories = os.listdir(PATH)
        hist_list = create_histograms()
        for directory in directories:
            files = os.listdir( PATH+"/"+directory ) 
            data_PATH=os.path.join(PATH, directory )
            for file in files:
                inFileName = os.path.join(data_PATH, file)
                with ROOT.TFile.Open(inFileName, "READ") as inFile:
                    tree = inFile.Get("muon/Events")
                    print(f"Reading from {inFileName}")
                    
                    entries = tree.GetEntries()
                    for i in range(entries):
                        process_event(tree, hist_list,i)

    

    # Define bin size and center for pT
    pt_bin_size = np.array([5.0, 5.0, 10.0, 15.0, 25.0, 25.0, 75.0, 75.0, 250.0])
    pt_bin_centers = np.array([35.0, 45.0, 60.0, 85.0, 125.0, 175.0, 275.0, 425.0, 750.0])  
    
    # Initialize arrays to store the results for pT
    num_bins = len(pt_bin_size)
    average_resolution = np.zeros(num_bins)
    sigma = np.zeros(num_bins)
    res_error = np.zeros(num_bins)
    sigma_error = np.zeros(num_bins)
    
    # Define bin size and center for dz
    dz_bin_size=np.array([2.5,2.5,5.0,5.0,7.5,7.5,10.0,35.0])
    dz_bin_centers=np.array([2.5,7.5,15.0,25.0,37.5,52.5,70.0,115.0])

    # Initialize arrays to store the results for dz
    dz_num_bins = len(dz_bin_size)
    dz_average_resolution = np.zeros(dz_num_bins)
    dz_sigma = np.zeros(dz_num_bins)
    dz_res_error = np.zeros(dz_num_bins)
    dz_sigma_error = np.zeros(dz_num_bins)
    
    # Define bin size and center for dxy
    dxy_bin_size=np.array([2.5,2.5,2.5,2.5,5.0,5.0,5.0,15.0])
    dxy_bin_centers=np.array([2.5,7.5,12.5,17.5,25.0,35.0,45.0,65.0])
    
    # Initialize arrays to store the results for dxy
    dxy_num_bins = len(dxy_bin_size)
    dxy_average_resolution = np.zeros(dxy_num_bins)
    dxy_sigma = np.zeros(dxy_num_bins)
    dxy_res_error = np.zeros(dxy_num_bins)
    dxy_sigma_error = np.zeros(dxy_num_bins)
    
    
    # Fit histograms for pt bins
    for i in range(len(pt_bin_size)):
        fit_histogram(hist_list[i+7], i, average_resolution, sigma, res_error, sigma_error)
        
    # Fit histograms for dz bins
    for i in range(len(dz_bin_size)):
        fit_histogram(hist_list[i+16], i, dz_average_resolution, dz_sigma, dz_res_error, dz_sigma_error)
    
    # Fit histograms for dxy bins
    for i in range(len(dxy_bin_size)):
        fit_histogram(hist_list[i+24], i, dxy_average_resolution, dxy_sigma, dxy_res_error, dxy_sigma_error)
    

    

    


    # Write histograms to the output ROOT file
    with ROOT.TFile.Open(f"Res_hist_{data_type}.root", "RECREATE") as outHistFile:
        for hist in hist_list:
            hist.Write()
            # Create and save graphs
        create_graph(f'Average_Resolution_{data_type}', 'p_{T} \mu_{ref} [GeV]', 'Mean of q/p_{T} relative residual', num_bins, pt_bin_centers, average_resolution, pt_bin_size, res_error, "canvas1", f"Average_Resolution_{data_type}")
        create_graph(f'Sigma_{data_type}', 'p_{T} \mu_{ref} [GeV]', 'Width of q/p_{T} relative residual', num_bins, pt_bin_centers, sigma, pt_bin_size, sigma_error, "canvas2", f"Sigma_{data_type}")
        
        create_graph(f'Average_Resolution_dz_{data_type}', '|dz|', 'Mean of q/p_{T} relative residual', dz_num_bins, dz_bin_centers, dz_average_resolution, dz_bin_size, res_error, "canvas3", f"Average_Resolution_dz_{data_type}")
        create_graph(f'Sigma_dz_{data_type}', '|dz|', 'Width of q/p_{T} relative residual', dz_num_bins, dz_bin_centers, dz_sigma, dz_bin_size, dz_sigma_error, "canvas4", f"Sigma_dz_{data_type}")
        
        create_graph(f'Average_Resolution_dxy_{data_type}', '|dxy|', 'Mean of q/p_{T} relative residual', dxy_num_bins, dxy_bin_centers, dxy_average_resolution, dxy_bin_size, dxy_res_error, "canvas5", f"Average_Resolution_dxy_{data_type}")
        create_graph(f'Sigma_dxy_{data_type}', '|dxy|', 'Width of q/p_{T} relative residual', dxy_num_bins, dxy_bin_centers, dxy_sigma, dxy_bin_size, dxy_sigma_error, "canvas6", f"Sigma_dxy_{data_type}")

if __name__ == "__main__":
    main()
