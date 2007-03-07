// Test routine to make sure there are no memory leaks when creating
// and deleting mcRootData objects.
// To Run:
// 1) Make sure mcRootData.dll is in your ROOT library path
//    a) This can be done by either modifying your local .rootrc file 
//       Unix.*.Root.DynamicPath:    $(ROOTANALYSIS)/lib
//       WinNT.*.Root.DynamicPath:   $(ROOTANALYSIS)/lib
//   OR
//    b) Copy mcRootData.so (or .dll) into the directory from where you start ROOT.
// 2) You want to make sure that memory statistics are kept in ROOT by modifying your
//    local .rootrc file, and setting:
//    Root.MemStat:            1
//    Root.ObjectStat:         1
// 3) Start ROOT
// 4) At the ROOT prompt, type:  ".x leaktest.cxx" 

{
    UInt_t numEvents = 500;
    UInt_t numParticles = 100;
    UInt_t runNum = 1;
	Int_t sourceId = 7;
	UInt_t sequence = 4;

    gObjectTable->Print();
    
    if (!strcmp(gSystem->GetName(), "WinNT")) {
      gSystem->Load("mcRootData.dll");
    } else {
      gSystem->Load("libmcRootData.so");
    }
    TFile *f =  new TFile("mc.root", "RECREATE");
    TTree *t = new TTree("Mc", "Mc");
    McEvent *ev = new McEvent();
    t->Branch("McEvent", "McEvent", &ev, 64000, 1);
    
    gObjectTable->Print();
    TRandom randGen ;
    RAND_NUM = randGen.Rndm() ;

    Int_t ievent;
    for (ievent = 0; ievent < numEvents; ievent++) {
        ev->Fake(ievent,runNum,RAND_NUM) ;

        t->Fill();
        ev->Clear();
    }
    
    delete ev;
    
    printf("Here is the object table after creating events, storing them to file and deleting the objects\n");
    gObjectTable->Print();
    
    f->Write();
    f->Close();
    delete f;
    
}


