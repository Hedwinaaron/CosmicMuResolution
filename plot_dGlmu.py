import ROOT

# Open the two ROOT files
file1_data = ROOT.TFile.Open("Cosmics_muons_DATA.root", "READ")
file2_mc = ROOT.TFile.Open("Cosmics_muons_MC.root", "READ")

def plot_comparison(plot1, plot2, tittle, x_axis_title, y_axis_title,y_low_limit,y_max_limit, file_name):
    # Retrieve graphs by name
    graph1_data = file1_data.Get(plot1)
    graph2_mc = file2_mc.Get(plot2)
    
    # Create a TMultiGraph
    multi_graph = ROOT.TMultiGraph()
    
    # Customize graph styles
    graph1_data.SetMarkerStyle(20)
    graph1_data.SetMarkerColor(ROOT.kBlack)
    graph1_data.SetLineColor(ROOT.kBlack)
    
    graph2_mc.SetMarkerStyle(21)
    graph2_mc.SetMarkerColor(ROOT.kRed)
    graph2_mc.SetLineColor(ROOT.kRed)
    
    # Add graphs to the TMultiGraph
    multi_graph.Add(graph1_data, "P")  # "P" means draw points
    multi_graph.Add(graph2_mc, "P")  # "L" means draw lines
    
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", "Comparison", 800, 600)
    canvas.SetGrid()
    canvas.SetLogx()
    canvas.SetLeftMargin(0.24)
    # Draw the TMultiGraph
    multi_graph.Draw("A")
    
    # Set axis titles
    multi_graph.SetTitle(f'{tittle};{x_axis_title};{y_axis_title}')
    # Adjust text sizes for axis titles and labels
    multi_graph.GetXaxis().SetTitleSize(0.04)  # Set X-axis title size
    multi_graph.GetYaxis().SetTitleSize(0.04)  # Set Y-axis title size
    multi_graph.GetXaxis().SetTitleOffset(1.2)
    #multi_graph.GetXaxis().SetLabelOffset(0.2)
    #multi_graph.GetYaxis().SetRangeUser(y_low_limit,y_max_limit)
    multi_graph.GetXaxis().SetRangeUser(0.0,1000)
    # Add a legend
    legend = ROOT.TLegend(0.9,0.6,1.0,0.8)
    legend.AddEntry(graph1_data, "Data ", "pl")
    legend.AddEntry(graph2_mc, "MC", "pl")
    legend.SetFillStyle(0)  # 0 means no fill
    legend.SetBorderSize(0)  # Remove the border
    legend.Draw()
    
    # Save the canvas
    canvas.SaveAs(f"{file_name}.png")

#pt
plot_comparison("mean_total", "mean_total", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Mean of q/p_{T} relative residual  ",-0.01,0.04, "Average_Resolution")
plot_comparison("Sigma_total", "Sigma_total", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Width of q/p_{T} relative residual  ",0 , 0.15, "Sigma")

#dz

#plot_comparison("Average_Resolution_dz_Data", "Average_Resolution_dz_MC", "", " |dz| ", "Mean of q/p_{T} relative residual  ", "Average_Resolution_dz")
#plot_comparison("Sigma_dz_Data", "Sigma_dz_MC", "", " |dz| ", "Width of q/p_{T} relative residual  ", "Sigma_dz")

#dxy
#plot_comparison("Average_Resolution_dxy_Data", "Average_Resolution_dxy_MC", "", " |dxy| ", "Mean of q/p_{T} relative residual  ", "Average_Resolution_dxy")
#plot_comparison("Sigma_dxy_Data", "Sigma_dxy_MC", "", " |dxy| ", "Width of q/p_{T} relative residual  ", "Sigma_dxy")


