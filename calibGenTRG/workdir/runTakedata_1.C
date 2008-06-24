//
// Example of using RTAutil to dump event info 
//
{
gROOT->Reset();    // roll back CINT context to last Save

Int_t irun = 77016135;

const char digifile[200], reconfile[200], calfile[200];
sprintf(digifile, "/nfs/farm/g/glast/u52/Integration/rootData/077016874/v13r6p1/digi/digitization-v4r0p2_077016874_digi_DIGI.root");
sprintf(reconfile, "/nfs/farm/g/glast/u52/Integration/rootData/077016874/v13r6p1/recon/recon-v4r0p2_077016874_recon_RECON.root");
sprintf(calfile,"/nfs/farm/g/glast/u52/Integration/rootData/077016874/v13r6p1/recon/recon-v4r0p2_077016874_cal_ntuple.root");

 
// Creating an empty TimeStruct Class
takedata_mt *r = new takedata_mt();

// Initialize the input files (Note: this doesn't initialize histograms)
r->Init(digifile,reconfile);
r->initCalTuple(calfile);
r->setParameters(1,45,7,11);
r->useKalmanCut(true);
r->useToTCut(true);
 r->inithistos();
r->Go(10000000);
r->savehistos("acdtest-1.root");
r->writeoutresults("acdtest_1.txt");
}
