// Test routine to make sure there are no memory leaks when creating
// and deleting digiRootData objects.
// To Run:
// 1) Make sure digiRootData.dll is in your ROOT library path
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
    UInt_t numXtals = 100;
    UInt_t runNum = 1;

    gObjectTable->Print();
    
    gSystem->Load("digiRootData.dll");
    TFile *f =  new TFile("digi.root", "RECREATE");
    TTree *t = new TTree("Digi", "Digi");
    DigiEvent *ev = new DigiEvent();
    t->Branch("DigiEvent", "DigiEvent", &ev, 64000, 1);
    
    gObjectTable->Print();
    

    Int_t ievent, ixtal;
    for (ievent = 0; ievent < numEvents; ievent++) {

        ev->initialize(ievent, runNum);

        for (ixtal = 0; ixtal < numXtals; ixtal++) {
            CalDigi *cal = new CalDigi();
            CalXtalId::CalTrigMode mode = CalXtalId::BESTRANGE;
            Short_t tower = 5;
            Short_t layer = 4;
            Short_t col = 3;
            CalXtalId xtalId(tower, layer, col);
            cal->initialize(mode, xtalId);
            ev->addCalDigi(cal);
            Char_t rangeM = CalXtalId::LEX8;
            Char_t rangeP = CalXtalId::HEX8;
            UShort_t adcM = 4095;
            UShort_t adcP = 4095;
            cal->addReadout(rangeP, adcP, rangeM, adcM);
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


