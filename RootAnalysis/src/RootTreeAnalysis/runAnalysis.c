{

     gROOT->Reset();    // roll back CINT context to last Save
     gROOT->LoadMacro("RootTreeAnalysis.cxx");  // load your code
     // Update for your Root File

    RootTreeAnalysis* m= new RootTreeAnalysis(
	       "nsbf_r000053_20010804_072159_ivte-raw.root",
	       "",
	       "");
// Latest Balloon Flight raw file before the flight. 
// It can be found at glast.05.
}







