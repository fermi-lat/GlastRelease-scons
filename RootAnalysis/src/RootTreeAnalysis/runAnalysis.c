{

     gROOT->Reset();    // roll back CINT context to last Save
     gROOT->LoadMacro("RootTreeAnalysis.cxx");  // load your code
     // Update for your Root File

    RootTreeAnalysis* m= new RootTreeAnalysis(
               "digi.root",
	       "",
	       "");
}







