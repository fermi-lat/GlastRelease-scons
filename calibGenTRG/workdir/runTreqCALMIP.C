void runTreqCALMIP(ULong_t run = 236198321)
 {
gROOT->Reset();    // roll back CINT context to last Save

// TReq_TkrCal_Mip runs: 236409225, 236225703, 236198321

const char digifile[200], reconfile[200],calfile[200], outHTML[256], outROOT[256];

sprintf(digifile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.57/digi/r0%d_v000_digi.root", run);
sprintf(reconfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.57/recon/r0%d_v000_recon.root", run);
sprintf(calfile,"root://glast-rdr.slac.stanford.edu//glast/Data/Flight/Level1/LPA/prod/1.57/cal/r0%d_v000_cal.root", run);

// Creating an empty TimeStruct Class
treqCAL *r = new treqCAL();
// Initialize the input files (Note: this doesn't initialize histograms)
r->Init(digifile,reconfile);
r->initCalTuple(calfile);
r->useMipCut(true);
r->useKalmanCut(true);
r->useToTCut(true);
r->useOneTrack(true);

r->setTotCut(0.9, 1.4); //instead of 1, 1.6
r->setMipCut(4);
// No need to initialize histograms for this job
r->inithistos();
r->setParameters(5,0);

// Set the tower mask 
 r->Go(1000000);
 
sprintf(outHTML,"r0%d_TReqCALMIP.html", run);
sprintf(outROOT,"r0%d_TReqCALMIP.root", run);

r->writeoutresults(outHTML,outROOT);
// r->writeoutresults("test_TReqCALMIP.html","test_TReqCALMIP.root");

}
