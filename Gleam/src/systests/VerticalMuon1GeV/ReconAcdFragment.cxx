void RootTreeAnalysis::ReconAcdHistDefine() {
    
    histFile->cd();
    
    float maxEnergy = 20.;

    TH1F *ACDTILECOUNT = new TH1F("ACDTILECOUNT", "Acd Tile Count",
        5, 0, 5);
    ACDTILECOUNT->SetXTitle("Number");

    TH1F *ACDENERGY = new TH1F("ACDENERGY", "Acd Total Energy Deposit",
        50, 0,maxEnergy );
    ACDENERGY->SetXTitle("Energy (MeV)");

    TH1F *ACDDOCA = new TH1F("ACDDOCA", "Acd DOCA",
        50, 0, 500);
    ACDDOCA->SetXTitle("Distance (mm)");

    TH1F *ACDACTDIST = new TH1F("ACDACTDIST", "Acd Active Distance",
        50, -2500, 250);
    ACDACTDIST->SetXTitle("Distance (mm)");

    TH1F *ACDGAMMADOCA = new TH1F("ACDGAMMADOCA", "Acd Gamma DOCA",
        50, 0, 250);
    ACDGAMMADOCA->SetXTitle("Distance (mm)");



}

void RootTreeAnalysis::ReconAcd() {
    
  AcdRecon* acd = rec->getAcdRecon();

// get multiplicities

  ((TH1F*)GetObjectPtr("ACDTILECOUNT"))->Fill(acd->getTileCount());
  ((TH1F*)GetObjectPtr("ACDENERGY"))->Fill(acd->getEnergy());
  ((TH1F*)GetObjectPtr("ACDDOCA"))->Fill(acd->getDoca());
  ((TH1F*)GetObjectPtr("ACDACTDIST"))->Fill(acd->getActiveDist());
  ((TH1F*)GetObjectPtr("ACDGAMMADOCA"))->Fill(acd->getGammaDoca());

};
