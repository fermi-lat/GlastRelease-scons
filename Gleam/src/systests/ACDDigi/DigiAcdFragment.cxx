void RootTreeAnalysis::DigiAcdHistDefine() {
    
    histFile->cd();
    TH1F *ACDNUMDIGI = new TH1F("ACDNUMDIGI", "Number of Acd Digis",
        101, -0.5, 100.5);


};

void RootTreeAnalysis::DigiAcd() {
    
    ((TH1F*)GetObjectPtr("ACDNUMDIGI"))->Fill(evt->getAcdDigiCol()->GetEntries());
};

