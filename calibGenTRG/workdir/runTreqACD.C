void runTreqACD(ULong_t run = 377616981)
{
gROOT->Reset();    // roll back CINT context to last Save
gROOT->ProcessLine(".x $WORKDIR/RooLogon.C");

const char digifile[200], reconfile[200], outHTML[256], outROOT[256];
int tkrDelay=0, acdDelay=0;

if(run==238477768){
  // L&EO 238477768
  //root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/merit/r0238477768_v002_merit.root 
  sprintf(digifile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/digi/r0%d_v002_digi.root", run);
  sprintf(reconfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.61/recon/r0%d_v002_recon.root", run);
  tkrDelay=5;
  acdDelay=31;
}
else if(run==377616981){
  // New 377616981
  sprintf(digifile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/4.1/digi/r0%d_v000_digi.root", run);
  sprintf(reconfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/4.1/recon/r0%d_v000_recon.root", run);
  tkrDelay=4; 
  acdDelay=31;
}
else{
printf("Unknow run");
break;
}
// Creating an empty TimeStruct Class
treqACD *r = new treqACD();
// Initialize the input files (Note: this doesn't initialize histograms)
r->Init(digifile,reconfile);
// No need to initialize histograms for this job
r->inithistos();
r->setParameters(tkrDelay, acdDelay);

// Process
r->Go(2000000);

// Save results
char outputdir[256]=gSystem->ExpandPathName("$OUTPUTDIR"); 
sprintf(outHTML,"%s/r0%d_TReqACD.html", outputdir, run);
sprintf(outROOT,"%s/r0%d_TReqACD.root", outputdir, run);
r->writeoutresults(outHTML,outROOT);
}
