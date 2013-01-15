void runTreqCALGamma(ULong_t run = 377604976)
 {
gROOT->Reset();    // roll back CINT context to last Save
gROOT->ProcessLine(".x $WORKDIR/RooLogon.C");

const char digifile[200], reconfile[200],calfile[200], outHTML[256], outROOT[256];
int tkrDelay=0, acdDelay=0;

if(run==238465830){
  // L&EO 238465830
  // root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/merit/r0238465830_v000_merit.root 
  sprintf(digifile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/digi/r0%d_v000_digi.root", run);
  sprintf(reconfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/recon/r0%d_v000_recon.root", run);
  sprintf(calfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/cal/r0%d_v000_cal.root", run);
  tkrDelay=5;
  acdDelay=16;
} else if(run==377604976) {
  // New 377604976
  //root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/4.1/merit/r0377604976_v001_merit.root 
  sprintf(digifile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/4.1/digi/r0%d_v001_digi.root", run);
  sprintf(reconfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/4.1/recon/r0%d_v000_recon.root", run);
  sprintf(calfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/4.1/cal/r0%d_v001_cal.root", run);
  tkrDelay=4;
  acdDelay=15;
} else {
break;
}
// Creating an empty TimeStruct Class
treqCAL *r = new treqCAL();
// Initialize the input files (Note: this doesn't initialize histograms)
r->Init(digifile,reconfile);
r->initCalTuple(calfile);
r->useGammaCut(true);
// No need to initialize histograms for this job
r->inithistos();
r->setParameters(tkrDelay, 0);

// Process
r->Go(2000000);

// Save results
char outputdir[256]=gSystem->ExpandPathName("$OUTPUTDIR");
sprintf(outHTML,"%s/r0%d_TReqCALGamma.html", outputdir, run);
sprintf(outROOT,"%s/r0%d_TReqCALGamma.root", outputdir, run);

r->writeoutresults(outHTML,outROOT);
}
