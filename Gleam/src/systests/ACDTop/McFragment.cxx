
void RootTreeAnalysis::McHistDefine() {
  // Purpose and Method:  Monte Carlo histogram definitions
  // Must cd into the histFile to be sure the histograms are created within
  // this histogram file.
    histFile->cd();
    
    TH1F *EVENTMC = new TH1F("EVENTMC", "MC Event Id",
        2000, 0, 2000);
    EVENTMC->SetBit(TH1::kCanRebin);
    
    TH1F *RUNMC = new TH1F("RUNMC", "MC Run Id",
        2000, 0, 2000);
    
    TH1F *MCPARTDAUGHTER = new TH1F("MCPARTDAUGHTER", "Number of Daughters",
        101, -0.5, 100.5);
    MCPARTDAUGHTER->SetBit(TH1::kCanRebin);
  
    TH1F *NUMELECPEREVENT = new TH1F("NUMELECPEREVENT", "Number of Electrons produced in an event", 51, -0.5, 50.5);
    NUMELECPEREVENT->SetBit(TH1::kCanRebin);

    TH1F *MCPARTTYPES = new TH1F("MCPARTTYPES", "All Monte Carlo Particle Types",     61, -30.5, 30.5);
    MCPARTTYPES->SetBit(TH1::kCanRebin);

    TH1F *MCELECENG = new TH1F("MCELECENG", "Energy of Secondary Electrons",
      200, 0, 10);
    MCELECENG->SetBit(TH1::kCanRebin);

    TH1F *ACDMCHIT22 = new TH1F("ACDMCHIT22", "Total Energy Deposition in Tile 022", 100, 0, 5);
    ACDMCHIT22->SetBit(TH1::kCanRebin);

    TH1F *ACDELECENG22 = new TH1F("ACDELECENG22", "Total Energy Deposition of Secondary Electrons in Tile 0022", 100, 0, 5);
    ACDELECENG22->SetBit(TH1::kCanRebin);

    TH1F *POSCOUNTMC = new TH1F("POSCOUNTMC", "MC PosHit Count",
        101, -0.5, 100.5);
    
    TH1F *POSENERGYDEP = new TH1F("POSENERGYDEP", "MC Pos Edep",
        25, 0, 0.5);
    
    TH1F *ENERGYMC = new TH1F("ENERGYMC", "MC Part Energy (MeV)",
        10001, -0.5, 10000.5);   
    
}


void RootTreeAnalysis::McData() {
    // Purpose and Method:  Process on Monte Carlo event
    
    // Update histograms which keep simple counts
    ((TH1F*)GetObjectPtr("EVENTMC"))->Fill((Float_t)mcEventId);
    ((TH1F*)GetObjectPtr("RUNMC"))->Fill((Float_t)mcRunNum);
    ((TH1F*)GetObjectPtr("POSCOUNTMC"))->Fill((Float_t)mc->getMcPositionHitCol()->GetEntries());
    
    McParticle *mom = (McParticle*)mc->getMcParticleCol()->At(0);
    ((TH1F*)GetObjectPtr("ENERGYMC"))->Fill(mom->getInitialFourMomentum().E());

    /// The 2nd McParticle in the list that represents the 
    // initial particle that interacts in the system
    McParticle *mom1 = (McParticle*)mc->getMcParticleCol()->At(1);
    ((TH1F*)GetObjectPtr("MCPARTDAUGHTER"))->Fill(mom1->getDaughterList().GetEntries());

    UInt_t numElec = 0;

    // Count the particle types
    const TObjArray* mcPartCol = mc->getMcParticleCol();
    TIter mcHitIt(mcPartCol);
    McParticle *mcPart = 0;
    Int_t flag = 0;
    while (mcPart = (McParticle*)mcHitIt.Next()) {
      if (flag) {
        ((TH1F*)GetObjectPtr("MCPARTTYPES"))->Fill(mcPart->getParticleId());
        if (mcPart->getParticleId() == 11) {  // electron
          ++numElec;
          ((TH1F*)GetObjectPtr("MCELECENG"))->Fill(mcPart->getInitialFourMomentum().E());

        }
           
      } else {
        flag = 1;
      }
    } // end McParticleCol while



    ((TH1F*)GetObjectPtr("NUMELECPEREVENT"))->Fill(numElec);

    const TObjArray* posHitCol = mc->getMcPositionHitCol();
    TIter posHitIt(posHitCol);

    McPositionHit* posHit = 0;
    Float_t totE = 0.0, totElecEng = 0.0;

    // Iterate through all McPositionHits to determine the total amount
    // of energy deposited in a specific ACD tile
    while( posHit = (McPositionHit*)posHitIt.Next()) {
      VolumeIdentifier volId = posHit->getVolumeId();
      if( volId[0] != 1 ) continue;  // Skip non ACD tiles

      // Check that this is the center ACD top tile
      if ((volId[1] == 0) && (volId[2] == 2) && (volId[3] == 2) ){
          totE += posHit->getDepositedEnergy();
          if (posHit->getMcParticleId() == 11) totElecEng += posHit->getDepositedEnergy();

      } // end tile 0022

    }

  
  ((TH1F*)GetObjectPtr("ACDMCHIT22"))->Fill(totE);
  ((TH1F*)GetObjectPtr("ACDELECENG22"))->Fill(totElecEng);


}

