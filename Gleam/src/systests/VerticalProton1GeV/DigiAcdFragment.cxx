void RootTreeAnalysis::DigiAcdHistDefine() {
    
    histFile->cd();

    float maxAdc = 400.;   // histo max range
    
    TH1F *ACDDIGICOUNT = new TH1F("ACDDIGICOUNT", "AcdDigi multiplicity",
        50, 0, 50);
    ACDDIGICOUNT->SetXTitle("Number");

    TH1F *ACDADC = new TH1F("ACDADC", "Acd Digi ADC",
        200, 0, maxAdc);
    ACDADC->SetXTitle("ADC");

    TH1F *ACDADCTOT = new TH1F("ACDADCTOT", "Acd Digi ADC summed over all tiles",
        200, 0, maxAdc);
    ACDADCTOT->SetXTitle("ADC");

    TH1F *ACDLAYER = new TH1F("ACDLAYER", "Cal Digi Layer",
        15, 0, 15);
    ACDLAYER->SetXTitle("Layer");

    TH1F *ACDCOLUMN = new TH1F("ACDCOLUMN", "Cal Digi Column",
        20, 0, 20);
    ACDCOLUMN->SetXTitle("Column");

    TH1F *ACDFACE = new TH1F("ACDFACE", "Cal Digi Face",
        20, 0, 20);
    ACDFACE->SetXTitle("Face");

    TH1F *ACDROW = new TH1F("ACDROW", "Cal Digi Row",
        20, 0, 20);
    ACDROW->SetXTitle("Row");


};

void RootTreeAnalysis::DigiAcd() {
    

// get multiplicities

    TObjArray* cL = evt->getAcdDigiCol();
    int nAcdDigi = cL->GetEntries();
    ((TH1F*)GetObjectPtr("ACDDIGICOUNT"))->Fill((Float_t)nAcdDigi);

    float eTotal = 0.;

    for (int ic=0; ic < nAcdDigi; ic++) {
      AcdDigi* c=(AcdDigi*)cL->At(ic);

      AcdId id = c->getId();

      int layer = id->getLayer();
      int column = id->getColumn();
      int face = id->getFace();
      int row = id->getRow();

      ((TH1F*)GetObjectPtr("ACDLAYER"))->Fill(layer);
      ((TH1F*)GetObjectPtr("ACDFACE"))->Fill(face);
      ((TH1F*)GetObjectPtr("ACDCOLUMN"))->Fill(column);
      ((TH1F*)GetObjectPtr("ACDROW"))->Fill(row);

      // get readout

      float adcA = c->getPulseHeight(AcdDigi::A);
      float adcB = c->getPulseHeight(AcdDigi::B);

      ((TH1F*)GetObjectPtr("ACDADC"))->Fill(adcA);
      ((TH1F*)GetObjectPtr("ACDADC"))->Fill(adcB);

      eTotal += adcA + adcB;
    }
      
    if (eTotal == 0) ((TH1F*)GetObjectPtr("ACDADCTOT"))->Fill(eTotal);

};

