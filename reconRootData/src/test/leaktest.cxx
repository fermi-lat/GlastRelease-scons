// Test routine to make sure there are no memory leaks when creating
// and deleting reconRootData objects.
// To Run:
// 1) Make sure reconRootData.dll is in your ROOT library path
//    a) This can be done by either modifying your local .rootrc file 
//       Unix.*.Root.DynamicPath:    $(ROOTANALYSIS)/lib
//       WinNT.*.Root.DynamicPath:   $(ROOTANALYSIS)/lib
//   OR
//    b) Copy reconRootData.so (or .dll) into the directory from where you start ROOT.
// 2) You want to make sure that memory statistics are kept in ROOT by modifying your
//    local .rootrc file, and setting:
//    Root.MemStat:            1
//    Root.ObjectStat:         1
// 3) Start ROOT
// 4) At the ROOT prompt, type:  ".x leaktest.cxx" 

{
    UInt_t numEvents = 500;
    UInt_t numXtals = 100;
    UInt_t runNum = 1;

    gObjectTable->Print();
    
    gSystem->Load("reconRootData.dll");
    TFile *f =  new TFile("recon.root", "RECREATE");
    TTree *t = new TTree("Recon", "Recon");
    ReconEvent *ev = new ReconEvent();
    t->Branch("ReconEvent", "ReconEvent", &ev, 64000, 1);
    
    gObjectTable->Print();
    

    Int_t ievent, ixtal;
    for (ievent = 0; ievent < numEvents; ievent++) {

       }
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


