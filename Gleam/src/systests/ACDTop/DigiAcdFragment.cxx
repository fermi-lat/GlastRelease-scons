void RootTreeAnalysis::DigiAcdHistDefine() {
    
    histFile->cd();

    gStyle->SetLabelSize(0.03);
    gStyle->SetTitleSize(0.03);
    gStyle->SetTitleFont(42);

    TH1F *ACDTILE0022PHA_A = new TH1F("ACDTILE0022PHA_A", "Tile 0022 PHA for PMT A", 817, -0.5, 4095.5);
    ACDTILE0022PHA_A->SetBit(TH1::kCanRebin);
    TH1F *ACDTILE0022PHA_B = new TH1F("ACDTILE0022PHA_B", "Tile 0022 PHA for PMT A", 817, -0.5, 4095.5);
    ACDTILE0022PHA_B->SetBit(TH1::kCanRebin);
    TH1F *ACDTILE0022VETO_A = new TH1F("ACDTILE0022VETO_A", "Tile 0022 Veto Disrim PMT A", 2, -0.5, 1.5);
    ACDTILE0022VETO_A->SetBit(TH1::kCanRebin);
    TH1F *ACDTILE0022VETO_B = new TH1F("ACDTILE0022VETO_B", "Tile 0022 Veto Disrim PMT B", 2, -0.5, 1.5);
    ACDTILE0022VETO_B->SetBit(TH1::kCanRebin);
    TH1F *ACDTILE0022VETO = new TH1F("ACDTILE0022VETO", "Tile 0022 Veto Disrim PMT A or B", 2, -0.5, 1.5);
    ACDTILE0022VETO->SetBit(TH1::kCanRebin);

};

void RootTreeAnalysis::DigiAcd() {
    
    const TObjArray* acdDigiCol = evt->getAcdDigiCol();
    if (!acdDigiCol) return;

    UShort_t pha0022_A = 0, pha0022_B = 0;
    // Create an id for tile 0022
    AcdId id22(0, 0, 2, 2);

    TIter acdDigiIter(acdDigiCol);
    AcdDigi *acdDigiItem = 0;

    while (acdDigiItem = (AcdDigi*)acdDigiIter.Next()) {
        AcdId id = acdDigiItem->getId();
        if (id.getId() == id22.getId()) {
            ((TH1F*)GetObjectPtr("ACDTILE0022PHA_A"))->Fill(acdDigiItem->getPulseHeight(AcdDigi::A)); 
            ((TH1F*)GetObjectPtr("ACDTILE0022PHA_B"))->Fill(acdDigiItem->getPulseHeight(AcdDigi::B));
            ((TH1F*)GetObjectPtr("ACDTILE0022VETO_A"))->Fill(acdDigiItem->getVeto(AcdDigi::A));
            ((TH1F*)GetObjectPtr("ACDTILE0022VETO_B"))->Fill(acdDigiItem->getVeto(AcdDigi::B));
            ((TH1F*)GetObjectPtr("ACDTILE0022VETO"))->Fill(acdDigiItem->getVeto(AcdDigi::A) || acdDigiItem->getVeto(AcdDigi::B));
        
        }
    }


};

