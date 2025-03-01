# CosmicMuResolution

## Selection
To apply the selection over the Data or MC samples, simply run:

    root -q -b DGM_sel.C

This will generate a ROOT file named either `Cosmics_muons_[data type]_[Muon type].root` or `Cosmics_muons_MC_DGL.root`, depending on the option specified in the code.
The `[Data type]` can be `DATA` or `MC` and is specified here:
https://github.com/Hedwinaaron/CosmicMuResolution/blob/af11f4074f2fbc550810ac6a8c441f435de4c9eb/DGM_sel.C#L29
This determines which files will be read based on the selected option. The list of files for `DATA` and `MC` can be modified here:
https://github.com/Hedwinaaron/CosmicMuResolution/blob/af11f4074f2fbc550810ac6a8c441f435de4c9eb/DGM_sel.C#L70-L82
The `[Muon type]` can be `DGL` or `DSA` and it is specified here:
https://github.com/Hedwinaaron/CosmicMuResolution/blob/af11f4074f2fbc550810ac6a8c441f435de4c9eb/DGM_sel.C#L31
This determines if the muons used are Displaced global muons `[DGL]`, or standalone muons `[DSA]`.
## Plot Muon Resolution

To plot a resolution comparison between the MC and Data samples, both `Cosmics_muons_[data type]_[Muon type].root` files must exist in the same directory. If both files are present, run the following command:

    python3 plot_dGlmu.py

This will generate several comparison plots.
