
//
// Example of using RTAutil to dump event info 
//
 {
gROOT->Reset();    // roll back CINT context to last Save


const char digifile[200], reconfile[200];

sprintf(digifile,"/nfs/farm/g/glast/u52/Integration/rootData/077017751/v13r6p1/digi/digitization-v4r0p2_077017751_digi_DIGI.root");
sprintf(reconfile,"/nfs/farm/g/glast/u52/Integration/rootData/077017751/v13r6p1/recon/recon-v4r0p2_077017751_recon_RECON.root");

// Creating an empty TimeStruct Class
treqCAL *r = new treqCAL();
// Initialize the input files (Note: this doesn't initialize histograms)
r->Init(digifile,reconfile);
// No need to initialize histograms for this job
r->inithistos();
r->setParameters(5,0);

// Set the tower mask 

 r->Go(1000000);
 r->writeoutresults("testhist.html","testhist.root");
}
