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
    
    VolumeIdentifier id;

    Int_t ievent, ipart;
    Double_t timestamp=0.0;
    for (ievent = 0; ievent < numEvents; ievent++) {

        ev->initialize(ievent, runNum, sourceId, sequence, timestamp);

        for (ipart = 0; ipart < numParticles; ipart ++) {
            McParticle *mcPart = new McParticle();
            Float_t rand = gRandom->Rndm();
            TLorentzVector initMom(1., 2., 3.);
            TLorentzVector finalMom(ipart*rand+ievent, rand*ievent+ipart, 1., 2.);
            TVector3 initPos(ipart*rand, ipart*2.0*rand, ipart*4.0*rand);
            TVector3 finalPos(ipart*rand, ipart*rand, ipart*rand);
            mcPart->initialize(mcPart, 0, 0, initMom, finalMom, initPos, finalPos, "pair");
            ev->addMcParticle(mcPart);
            if (ipart % 2) {
                McPositionHit *posHit = new McPositionHit();
                TVector3 entry(1., 1., 1.);
                TVector3 exit(0., 0.5, 0.1);
                TVector3 gEntry(3., 3., 3.);
                TVector3 gExit(0., 1.5, 0.2);
                id.initialize(0, 1, 1);
                posHit->initialize(ipart, 0, rand, id, 
                    entry, exit, gEntry, gExit,
                    0, 0, rand*0.1, rand*0.4, 0);
                ev->addMcPositionHit(posHit);
            } else {
                McIntegratingHit *intHit = new McIntegratingHit();
                id.initialize(0, 0, 0);
                intHit->initialize(id);
                TVector3 pos = mcPart->getFinalPosition();
                intHit->addEnergyItem(1.5, mcPart, pos);
                intHit->addEnergyItem(1.5, McIntegratingHit::PRIMARY, pos);
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


