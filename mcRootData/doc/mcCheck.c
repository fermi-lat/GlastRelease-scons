// THis little macro written by Dan Flath will do a very basic check on the
// contents of a RMcEvent written to disk.
//
// Note: you will have to change the locations of
//       the dll and of the root file.

{

  gSystem->Load("mcRootData.dll");// change me
  
  TFile *f = new TFile("pdr-mc.root");// change me
  TTree *t = (TTree*)f->Get("T;1");

  RMcEvent *e = new RMcEvent();
  t->SetBranchAddress("RMcEvent", &e);

  for (int i=0; i<t->GetEntries(); i++) {
     e->Clean();
     e->Create();
     t->GetEvent(i);
     printf("%u\n\n", i);
     e->Dump();
     printf("----------------\n");

     printf("i: %u,   event: %u,    particles: %u\n", i, e->getEvent(), e->getVertCount());
     printf("Incident Particle (%5.3f, %5.3f, %5.3f)\n", e->getVert(0)->getPosInitial()->X(),
		e->getVert(0)->getPosInitial()->Y(), e->getVert(0)->getPosInitial()->Z());
	
  }
}