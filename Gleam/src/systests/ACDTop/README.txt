Feb. 14, 2003

This system test concerns only the ACD.  Currently we are solely checking the
Monte Carlo deposition within a single ACD tile - in this case top tile 0022.

To Run this test, execute systest_ACDTop.  
This will start up 5 Gleam jobs and produce a set of 5 Monte Carlo and Digi 
ROOT files.  After the runs are completed, ROOT will start up and execute the 
following macro:  RUN_linux (or RUN_Windows if on a windows machine).  This in 
turn will execute the RUN_mcDep_linux (or RUN_mcDep_Windows) macro.  The 5 ROOT
files will be chained together, run through RootTreeAnalysis producing a 
Histograms.root file.  The Histograms file is then accessed and one TProfile
plot is drawn pertaining to the energy deposition versus angle of incidence.
