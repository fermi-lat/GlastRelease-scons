void runTakedata_1(ULong_t run = 236084237)
{
gROOT->Reset();    // roll back CINT context to last Save

// Take DATA runs: 236084237
//ULong_t run = 236084237; 

const char digifile[256], reconfile[256],calfile[256], outTXT[256], outROOT[256];

sprintf(digifile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.56/digi/r0%d_v000_digi.root", run);
sprintf(reconfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.56/recon/r0%d_v000_recon.root", run);
sprintf(calfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.56/cal/r0%d_v000_cal.root", run);

 
// Creating an empty TimeStruct Class
takedata_mt *r = new takedata_mt();

// Initialize the input files (Note: this doesn't initialize histograms)
r->Init(digifile,reconfile);
r->initCalTuple(calfile);
r->setParameters(1,45,7,11);
r->useMipCut(true);
r->useKalmanCut(true);
r->useToTCut(true);
r->useOneTrack(true);
r->setTotCut(0.9, 1.5); //instead of 1, 1.6

r->inithistos();
r->Go(10000);

sprintf(outROOT,"r0%d_acdtest-1.root", run);
sprintf(outTXT,"r0%d_acdtest_1.txt", run);
r->savehistos(outROOT);
r->writeoutresults(outTXT);

}
