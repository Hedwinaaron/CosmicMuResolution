import ROOT

# Open the two ROOT files
file1 = ROOT.TFile.Open("Res_hist_Data.root", "READ")
file2 = ROOT.TFile.Open("Res_hist_MC.root", "READ")

def plot_comparison(plot1, plot2, tittle, x_axis_title, y_axis_title, file_name):
    # Retrieve graphs by name
    graph1 = file1.Get(plot1)
    graph2 = file2.Get(plot2)
    
    # Create a TMultiGraph
    multi_graph = ROOT.TMultiGraph()
    
    # Customize graph styles
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerColor(ROOT.kRed)
    graph1.SetLineColor(ROOT.kRed)
    
    graph2.SetMarkerStyle(21)
    graph2.SetMarkerColor(ROOT.kBlue)
    graph2.SetLineColor(ROOT.kBlue)
    
    # Add graphs to the TMultiGraph
    multi_graph.Add(graph1, "P")  # "P" means draw points
    multi_graph.Add(graph2, "P")  # "L" means draw lines
    
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

    # Add a legend
    legend = ROOT.TLegend(0.8,0.6,0.9,0.8)
    legend.AddEntry(graph1, "Data ", "pl")
    legend.AddEntry(graph2, "MC", "pl")
    legend.Draw()
    
    # Save the canvas
    canvas.SaveAs(f"{file_name}.png")

#pt
plot_comparison("Average_Resolution_Data", "Average_Resolution_MC", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Mean of q/p_{T} relative residual  ", "Average_Resolution")
plot_comparison("Sigma_Data", "Sigma_MC", "", "p_{T} $\mu_{ref}$ [GeV]  ", "Width of q/p_{T} relative residual  ", "Sigma")

#dz

plot_comparison("Average_Resolution_dz_Data", "Average_Resolution_dz_MC", "", " |dz| ", "Mean of q/p_{T} relative residual  ", "Average_Resolution_dz")
plot_comparison("Sigma_dz_Data", "Sigma_dz_MC", "", " |dz| ", "Width of q/p_{T} relative residual  ", "Sigma_dz")

#dxy
plot_comparison("Average_Resolution_dxy_Data", "Average_Resolution_dxy_MC", "", " |dxy| ", "Mean of q/p_{T} relative residual  ", "Average_Resolution_dxy")
plot_comparison("Sigma_dxy_Data", "Sigma_dxy_MC", "", " |dxy| ", "Width of q/p_{T} relative residual  ", "Sigma_dxy")


