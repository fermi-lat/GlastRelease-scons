/*! Monte Carlo Analysis Histogram Defintions */

void RootTreeAnalysis::McHistDefine() {
    
    histFile->cd();
    
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
        100, 0, 1000);
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
    


    TH2F *INTMOMPOSMCX = new TH2F("INTCOUNTMCX", "McIntegratingHits vs MC x",
        200,-200,200,5, 0, 5);
    INTMOMPOSMCX->SetXTitle("position (mm)");
    INTMOMPOSMCX->SetYTitle("Number McIntegratingHits");
    
    TH2F *INTMOMMCXZ = new TH2F("INTMOMMCXZ", "x moment vs MC x",
        200,-200,200,200, -20, 20);
    INTMOMMCXZ->SetXTitle("MC position (mm)");
    INTMOMMCXZ->SetYTitle("moment (mm)");
    
    TH2F *INTSEGMENTMCX = new TH2F("INTSEGMENTMCX", "MC Int Xtal Segment vs MC x",
        20, 0, 12, 50, -200, 200);
    INTSEGMENTMCX->SetXTitle("MC position (mm)");
    INTSEGMENTMCX->SetYTitle("Segment number");

    TH2F *INTMCEDEPMCX = new TH2F("INTMCEDEPMCX", "MC Energy vs MC x",
         50, -200, 200, 50, 0, 50);
    INTMCEDEPMCX->SetXTitle("MC position (mm)");
    INTMCEDEPMCX->SetYTitle("Energy (MeV)");
        
}

/* Process the Monte Carlo Data
Called by Go()
*/
 float mcX=0.;
 int eXtal = 0;

void RootTreeAnalysis::McData() {
    
// get multiplicities

    ((TH1F*)GetObjectPtr("PARTCOUNTMC"))->Fill((Float_t)mc->getMcParticleCol()->GetEntries());
    ((TH1F*)GetObjectPtr("INTCOUNTMC"))->Fill((Float_t)mc->getMcIntegratingHitCol()->GetEntries());

// look at McParticles

     TObjArray* iP = mc->getMcParticleCol();
     int nP = iP->GetEntries();

     for (int ip=0; ip < nP; ip++) {
       McParticle* p = (McParticle*)iP->At(ip);
       ((TH1F*)GetObjectPtr("ENERGYMC"))->Fill(p->getInitialFourMomentum().E());
       if (ip==0) {
	 mcX = p->getInitialPosition().Px();
	 ((TH1F*)GetObjectPtr("MCX"))->Fill(mcX);
       }
     }

// look at McIntegratingHits

     TObjArray* iHL = mc->getMcIntegratingHitCol();

     int nH = iHL->GetEntries();
     float eTotCal = 0.;
     float eTotAcd = 0.;
     ((TH2F*)GetObjectPtr("INTCOUNTMCX"))->Fill(mcX,nH);

     for (int ih=0; ih < nH; ih++) {
       McIntegratingHit* hit = (McIntegratingHit*)iHL->At(ih);
       float eHit = hit->getTotalEnergy();

// use volume ID to pick ACD or CAL

       VolumeIdentifier volId = hit->getVolumeId();

       int segn = volId[8];
       eXtal = volId[7];
       float moment = hit->getMoment1().Px();
       ((TH1F*)GetObjectPtr("INTSEGMENT"))->Fill((float)segn);
       ((TH2F*)GetObjectPtr("INTSEGMENTMCX"))->Fill((float)segn,mcX);
       ((TH2F*)GetObjectPtr("INTMOMMCXZ"))->Fill(mcX,moment);

       if ((volId[0] == 0) &&    // CAL identifier
	   (volId[3] == 0)){ 
	 eTotCal += eHit;
	 ((TH1F*)GetObjectPtr("INTENERGYDEPCAL"))->Fill(eHit);
	 float momPos = 326.*(segn+0.5)/12 + moment - 163.;
	 ((TH1F*)GetObjectPtr("INTMOMPOS"))->Fill(momPos); 
       }
    }
     ((TH1F*)GetObjectPtr("INTENERGYTOT"))->Fill(eTotCal);
     ((TH2F*)GetObjectPtr("INTMCEDEPMCX"))->Fill(mcX,eTotCal);

}

