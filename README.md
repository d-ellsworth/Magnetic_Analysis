# Magnetic_Analysis
Various Matlab analysis tools for magnetics research.

## FMR data fitting
#### FMRdatafitting.m
'FMRdatafitting.m' is used for fitting the Lorentzian derivative to ferromagnetic resonance measurement data and extracting the FMR field and FMR linewidth. It will calculate Kittel equation and damping fits for frequency data.

#### PolarAngleFitting.m
`PolarAngleFitting.m` is used for fitting angle dependent FMR data. This file will use the output of 'FMRdatafitting.m' to run, see the code comments. It will fit Hfmr vs angle and fit Gilbert damping as a function of angle.

## PSV
`PSVanalysis.m` is used for cleaning up and analyzing time domain data from the photo-spin-voltaic effect.
