// Test routine to make sure there are no memory leaks when creating
// and deleting digiRootData objects.
// To Run:
// 1) Make sure digiRootData.dll is in your ROOT library path
//    a) This can be done by either modifying your local .rootrc file 
//       Unix.*.Root.DynamicPath:    $(ROOTANALYSIS)/lib
//       WinNT.*.Root.DynamicPath:   $(ROOTANALYSIS)/lib
//   OR
//    b) Copy digiRootData.so (or .dll) into the directory from where you start ROOT.
// 2) You want to make sure that memory statistics are kept in ROOT by modifying your
//    local .rootrc file, and setting:
//    Root.MemStat:            1
//    Root.ObjectStat:         1
// 3) Start ROOT
// 4) At the ROOT prompt, type:  ".x leaktest.cxx" 

{
    UInt_t numEvents = 500;
    UInt_t numXtals = 100;
    UInt_t numDigi = 150;
    UInt_t numAcd = 50;
    UInt_t runNum = 1;

    gObjectTable->Print();
    
    gSystem->Load("digiRootData.dll");
    TFile *f =  new TFile("digi.root", "RECREATE");
    TTree *t = new TTree("Digi", "Digi");
    DigiEvent *ev = new DigiEvent();
    t->Branch("DigiEvent", "DigiEvent", &ev, 64000, 1);
    
    gObjectTable->Print();
    

    Int_t ievent, ixtal;
    CalXtalId xtalId;

    for (ievent = 0; ievent < numEvents; ievent++) {

        ev->initialize(ievent, runNum, true);

        for (ixtal = 0; ixtal < numXtals; ixtal++) {
            CalDigi *cal = new CalDigi();
            CalXtalId::CalTrigMode mode = CalXtalId::BESTRANGE;
            Short_t tower = 5;
            Short_t layer = 4;
            Short_t col = 3;
            xtalId.init(tower, layer, col);
            cal->initialize(mode, xtalId);
            ev->addCalDigi(cal);
            Char_t rangeM = CalXtalId::LEX8;
            Char_t rangeP = CalXtalId::HEX8;
            UShort_t adcM = 4095;
            UShort_t adcP = 4095;
            cal->addReadout(rangeP, adcP, rangeM, adcM);
       }

        UInt_t idigi;
        for (idigi = 0; idigi < numDigi; idigi++) {
            TkrDigi *tkr = new TkrDigi();
            Int_t tot[2] = {idigi, idigi+1};
            TowerId id(3, 2);
            tkr->initialize(idigi, GlastAxis::Y, id, tot);
            UInt_t istrip;
            for (istrip = 0; istrip < idigi; istrip++) {
                tkr->addC0Hit(istrip);
            }
            tkr->addC1Hit(idigi*2);
            ev->addTkrDigi(tkr);
        }

        for (idigi = 0; idigi < numAcd; idigi++) {
            AcdId acd_id(0, 2, 5, 4);
            Float_t energy = ievent;
            Bool_t veto[2] = {kFALSE, kTRUE};
            Bool_t low[2] = {kTRUE, kTRUE};
            Bool_t high[2] = {kFALSE, kTRUE};
            UShort_t pha[2] = {4095, 1};
            ev->addAcdDigi(acd_id, energy, pha, veto, low, high);
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


