/*! Monte Carlo Analysis Histogram Defintions */

void RootTreeAnalysis::McHistDefine() {
    
    histFile->cd();

    float maxEnergy = 2000;
    
    TH1F *PARTCOUNTMC = new TH1F("PARTCOUNTMC", "McParticle Count",
        20, 0, 20);
    PARTCOUNTMC->SetXTitle("Number");

    TH1F *MCX = new TH1F("MCX", "MC x (mm) of original particle",
        20, -200, 200);
    MCX->SetXTitle("position (mm)");

    TH1F *MCY = new TH1F("MCY", "MC y (mm) of original particle",
        20, -200, 200);
    MCY->SetXTitle("position (mm)");

    TH1F *INTCOUNTMC = new TH1F("INTCOUNTMC", "McIntegratingHit Count",
        25, 0, 250);
    INTCOUNTMC->SetXTitle("Number");

    TH1F *INTENERGYDEPCAL = new TH1F("INTENERGYDEPCAL", "McIntegratingHit Ehit",
        100, 0, 10);
    INTENERGYDEPCAL->SetXTitle("Energy (MeV)");

    TH1F *INTENERGYTOT = new TH1F("INTENERGYTOT", "McIntegratingHit Total Energy",
        100, 0, maxEnergy);
    INTENERGYTOT->SetXTitle("Energy (MeV)");

    TH1F *INTSEGMENT = new TH1F("INTSEGMENT", "MC Int Xtal Segment Number",
        12, 0, 12);
    INTSEGMENT->SetXTitle("Segment Number");

    TH1F *INTMOMPOS = new TH1F("INTMOMPOS", "McIntegratingHit Position from moments",
        200, -200, 200);
    INTMOMPOS->SetXTitle("position (mm)");

    TH1F *ENERGYMC = new TH1F("ENERGYMC", "McParticle Initial Energy (MeV)",
        500, 0, 2000);   
    ENERGYMC->SetXTitle("Energy (MeV)");
    
    TH1F *MCPARTBITS = new TH1F("MCPARTBITS", "McParticle Status Bits",
        20, 0, 20);   
    MCPARTBITS->SetXTitle("Bit Number");
    
    TH1F *MCTERMZ = new TH1F("MCTERMZ", "MC Particle Termination z",
        300, -400, 800);   
    MCTERMZ->SetXTitle("position (mm)");
    
    TH1F *POSCOUNTMC = new TH1F("POSCOUNTMC", "McPositionHit Count",
        20, 0, 100);
    POSCOUNTMC->SetXTitle("Number");

    TH1F *POSENERGYDEP = new TH1F("POSENERGYDEP", "McPositionHit Deposited Energy",
        25, 0, 0.5);
    POSENERGYDEP->SetXTitle("Energy (MeV)");
    
    TH1F *POSMCTYPE = new TH1F("POSMCTYPE", "MC Pos MC pTypes",
        40, -20, 20);
    POSMCTYPE->SetXTitle("Particle type");

    TH1F *POSMCZ = new TH1F("POSMCZ", "McPositionHit local z position",
        200, -0.5, 0.5);
    POSMCZ->SetXTitle("position (mm)");

    TH1F *POSMCTKRTRAY = new TH1F("POSMCTKRTRAY", 
				  "McPositionHit Tkr Tray Number from Volume ID",
        20, 0, 20);
    POSMCTKRTRAY->SetXTitle("Tray Number");



    
}

/* Process the Monte Carlo Data
Called by Go()
*/


void RootTreeAnalysis::McData() {
    
// get multiplicities

    ((TH1F*)GetObjectPtr("PARTCOUNTMC"))->Fill((Float_t)mc->getMcParticleCol()->GetEntries());
    ((TH1F*)GetObjectPtr("INTCOUNTMC"))->Fill((Float_t)mc->getMcIntegratingHitCol()->GetEntries());
    ((TH1F*)GetObjectPtr("POSCOUNTMC"))->Fill((Float_t)mc->getMcPositionHitCol()->GetEntries());

// look at McParticles

     TObjArray* iP = mc->getMcParticleCol();
     int nP = iP->GetEntries();

     for (int ip=0; ip < nP; ip++) {
       McParticle* p = (McParticle*)iP->At(ip);

       ((TH1F*)GetObjectPtr("ENERGYMC"))->Fill(p->getInitialFourMomentum().E());
       ((TH1F*)GetObjectPtr("MCTERMZ"))->Fill( p->getFinalPosition().Pz());

       // check the status bits in StatusFlags. 

       int statusBits = p->getStatusFlags();
       for (int ib=0; ib<20; ib++) {
	 if (statusBits & 1<<ib) ((TH1F*)GetObjectPtr("MCPARTBITS"))->Fill(ib);
       }

       // get the x position from the initial particle 

       if (ip==0) {
	 float mcX = p->getInitialPosition().Px();
	 ((TH1F*)GetObjectPtr("MCX"))->Fill(mcX);
	 float mcY = p->getInitialPosition().Py();
	 ((TH1F*)GetObjectPtr("MCY"))->Fill(mcY);
       }
     }

// look at McIntegratingHits

     TObjArray* iHL = mc->getMcIntegratingHitCol();

     int nH = iHL->GetEntries();
     float eTotCal = 0.;
     float eTotAcd = 0.;

     for (int ih=0; ih < nH; ih++) {
       McIntegratingHit* hit = (McIntegratingHit*)iHL->At(ih);

       float eHit = hit->getTotalEnergy();

// use volume ID to pick ACD or CAL

       VolumeIdentifier volId = hit->getVolumeId();


       if ((volId[0] == 0) &&    // CAL identifier
	   (volId[3] == 0)){ 

	 eTotCal += eHit;

	 // get segment number from volumeID

	 int segn = volId[8];
	 float moment = hit->getMoment1().Px();
	 ((TH1F*)GetObjectPtr("INTSEGMENT"))->Fill((float)segn);

	 ((TH1F*)GetObjectPtr("INTENERGYDEPCAL"))->Fill(eHit);

	 // recalculate position from the segment and moment. Requires knowledge of geometry!

	 float momPos = 326.*(segn+0.5)/12 + moment - 163.;
	 ((TH1F*)GetObjectPtr("INTMOMPOS"))->Fill(momPos); 
       }
    }
     ((TH1F*)GetObjectPtr("INTENERGYTOT"))->Fill(eTotCal);

// look at McPositionHit

     TObjArray* iPL = mc->getMcPositionHitCol();
     int nP = iPL->GetEntries();
     for (int ip=0; ip < nP; ip++) {

       McPositionHit* pHit = (McPositionHit*)iPL->At(ip);
       ((TH1F*)GetObjectPtr("POSENERGYDEP"))->Fill((Float_t)pHit->getDepositedEnergy());

       int mcType = pHit->getMcParticleId();
       ((TH1F*)GetObjectPtr("POSMCTYPE"))->Fill(mcType);
       float z = pHit->getExitPosition()->z();
       ((TH1F*)GetObjectPtr("POSMCZ"))->Fill(z);
       VolumeIdentifier vID = pHit->getVolumeId();
       int tray = vID[4];
       ((TH1F*)GetObjectPtr("POSMCTKRTRAY"))->Fill(tray);
     }


}

