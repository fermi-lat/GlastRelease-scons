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

    gObjectTable->Print();
    
    gSystem->Load("mcRootData.dll");
    TFile *f =  new TFile("mc.root", "RECREATE");
    TTree *t = new TTree("Mc", "Mc");
    McEvent *ev = new McEvent();
    t->Branch("McEvent", "McEvent", &ev, 64000, 1);
    
    gObjectTable->Print();
    
    VolumeIdentifier id;

    Int_t ievent, ipart;
    for (ievent = 0; ievent < numEvents; ievent++) {

        ev->initialize(ievent, runNum);

        for (ipart = 0; ipart < numParticles; ipart ++) {
            McParticle *mcPart = new McParticle();
            Float_t rand = gRandom->Rndm();
            TLorentzVector initMom(1., 2., 3.);
            TLorentzVector finalMom(ipart*rand+ievent, rand*ievent+ipart, 1., 2.);
            TVector3 finalPos(ipart*rand, ipart*rand, ipart*rand);
            mcPart->initialize(mcPart, 0, 0, initMom, finalMom, finalPos);
            ev->addMcParticle(mcPart);
            if (ipart % 2) {
                McPositionHit *posHit = new McPositionHit();
                TVector3 entry(1., 1., 1.);
                TVector3 exit(0., 0.5, 0.1);
                id.initialize(0, 1, 1);
                posHit->initialize(rand, id, entry, exit, mcPart, 0, rand*0.1, rand*0.4, 0);
                ev->addMcPositionHit(posHit);
            } else {
                McIntegratingHit *intHit = new McIntegratingHit();
                id.initialize(0, 0, 0);
                intHit->initialize(id);
                intHit->addEnergyItem(1.5, mcPart, TVector3(0.5, 0.5, 0.5));
                ev->addMcIntegratingHit(intHit);
            }
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


