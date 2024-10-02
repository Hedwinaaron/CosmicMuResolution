# CosmicMuResolution
=====================

## Selection
To apply the selection over the Data or MC samples, simply run:

    python3 dGlobal-mu_sel.py

This will generate a ROOT file named either `Res_hist_Data.root` or `Res_hist_MC.root`, depending on the option specified in the code.

## Plot Muon Resolution

To plot a resolution comparison between the MC and Data samples, both `Res_hist_[type].root` files must exist in the same directory. If both files are present, run the following command:

    python3 plot_dGlmu.py

This will generate several comparison plots.
