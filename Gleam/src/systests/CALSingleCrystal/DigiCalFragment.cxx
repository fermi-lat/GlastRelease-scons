void RootTreeAnalysis::DigiCalHistDefine() {
    
  histFile->cd();

  //float maxAdc = 4000.;   // electrons histo max range
            float maxAdc = 1000.;   // muons histo max range
    
    TH1F *CALDIGICOUNT = new TH1F("CALDIGICOUNT", "CalDigi multiplicity",
        50, 0, 50);
    CALDIGICOUNT->SetXTitle("Number");

    TH1F *CALADC = new TH1F("CALADC", "Cal Digi ADC - both faces",
        200, 0, maxAdc);
    CALADC->SetXTitle("ADC");

    TH1F *CALADCN = new TH1F("CALADCN", "Cal Digi ADC - NEG",
        200, 0, maxAdc);
    CALADCN->SetXTitle("ADC");

    TH1F *CALADCP = new TH1F("CALADCP", "Cal Digi ADC - POS",
        200, 0, maxAdc);
    CALADCP->SetXTitle("ADC");

    TH2F *CALADCPM = new TH2F("CALADCPM", "Cal Digi ADC - NEG vs POS",
        200,0,maxAdc,200, 0, maxAdc);
    CALADCPM->SetXTitle("position (mm)");
    CALADCPM->SetYTitle("position (mm)");

    TH1F *CALRANGE = new TH1F("CALRANGE", "Cal Digi Hit Range - both faces",
        10, 0, 10);
    CALRANGE->SetXTitle("Range");

    TH1F *CALEAVE = new TH1F("CALEAVE", "CalDigi Energy - sum faces/2",
        200, 0, maxAdc);
    CALEAVE->SetXTitle("ADC");

    TH1F *CALEAVETOTAL = new TH1F("CALEAVETOTAL", "Cal Digi Total Energy - faces/2 summed",
        200, 0, 10*maxAdc);
    CALEAVETOTAL->SetXTitle("ADC");

    TH1F *CALLAYER = new TH1F("CALLAYER", "Cal Digi Layer",
        15, 0, 15);
    CALLAYER->SetXTitle("Layer");

    TH1F *CALTOWER = new TH1F("CALTOWER", "Cal Digi Tower",
        20, 0, 20);
    CALTOWER->SetXTitle("Tower");

    TH1F *CALCOLUMN = new TH1F("CALCOLUMN", "Cal Digi Column",
        20, 0, 20);
    CALCOLUMN->SetXTitle("Column");


    TH2F *CALADCNMCX = new TH2F("CALADCNMCX", "Cal Digi ADC - NEG vs MC x",
        50,-200,200,200, 0, maxAdc);
    CALADCNMCX->SetXTitle("position (mm)");
    CALADCNMCX->SetYTitle("ADC");

    TH2F *CALADCPMCX = new TH2F("CALADCPMCX", "Cal Digi ADC - POS vs MC x",
        50,-200,200,200, 0, maxAdc);
    CALADCPMCX->SetXTitle("position (mm)");
    CALADCPMCX->SetYTitle("ADC");

    TH1F *RESIDUAL = new TH1F("RESIDUAL", "Cal Asy - MC POS",
        200, -200, 200);
    RESIDUAL->SetXTitle("Residual (mm)");

    TH2F *RESIDUALMCX = new TH2F("RESIDUALMCX", "Cal Residual vs MC x",
        200,-200,200,200, -200, 200);
    RESIDUALMCX->SetXTitle("position (mm)");
    RESIDUALMCX->SetYTitle("Asymmetry");

    TH2F *CALASYMCX = new TH2F("CALASYMCX", "Cal Light asy vs MC x",
        200,-200,200,200, -0.4, 0.4);
    CALASYMCX->SetXTitle("position (mm)");
    CALASYMCX->SetYTitle("Asymmetry");

};

void RootTreeAnalysis::DigiCal() {
    

// get multiplicities

    TObjArray* cL = evt->getCalDigiCol();
    TObjArray* iP = mc->getMcParticleCol();
    int nCalDigi = cL->GetEntries();
    ((TH1F*)GetObjectPtr("CALDIGICOUNT"))->Fill((Float_t)nCalDigi);

    int nLayer[8]={0,0,0,0,0,0,0,0};
    float eLayer[8]={0.,0.,0.,0.,0.,0.,0.,0.};
    float eTotal = 0.;

    for (int ic=0; ic < nCalDigi; ic++) {
      CalDigi* c=(CalDigi*)cL->At(ic);
      McParticle* p = (McParticle*)iP->At(ic);

      CalXtalId id = c->getPackedId();
      int layer = id->getLayer();

      //      if (layer != 1) continue;

      CalXtalReadout* cRo=c->getXtalReadout(0);
      float adcP = cRo->getAdc(CalXtalId::POS); 
      float adcN = cRo->getAdc(CalXtalId::NEG); 

      float asy = (adcP-adcN)/(adcP+adcN-200.);   // subtract pedestal!
      float mcX = p->getInitialPosition().Px();
      float resid = 1452.*asy - mcX;
      ((TH1F*)GetObjectPtr("RESIDUAL"))->Fill(resid);
      ((TH2F*)GetObjectPtr("RESIDUALMCX"))->Fill(mcX,resid);

      float eAve = (adcN+adcP)/2.;
      ((TH1F*)GetObjectPtr("CALEAVE"))->Fill(eAve);
      ((TH1F*)GetObjectPtr("CALADC"))->Fill(adcN);
      ((TH1F*)GetObjectPtr("CALADC"))->Fill(adcP);
      ((TH1F*)GetObjectPtr("CALADCN"))->Fill(adcN);
      ((TH2F*)GetObjectPtr("CALADCNMCX"))->Fill(mcX,adcN);
      ((TH2F*)GetObjectPtr("CALADCPM"))->Fill(adcP,adcN);
      ((TH1F*)GetObjectPtr("CALADCP"))->Fill(adcP);
      ((TH2F*)GetObjectPtr("CALADCPMCX"))->Fill(mcX,adcP);
      ((TH2F*)GetObjectPtr("CALASYMCX"))->Fill(mcX,asy);
      ((TH1F*)GetObjectPtr("CALRANGE"))->Fill(cRo->getRange(0));
      ((TH1F*)GetObjectPtr("CALRANGE"))->Fill(cRo->getRange(1));

      int tower = id->getTower();
      int column = id->getColumn();
      ((TH1F*)GetObjectPtr("CALLAYER"))->Fill(layer);
      ((TH1F*)GetObjectPtr("CALTOWER"))->Fill(tower);
      ((TH1F*)GetObjectPtr("CALCOLUMN"))->Fill(column);

      nLayer[layer] += 1;
      eLayer[layer] += eAve;
      eTotal += eAve;
    }
      
    ((TH1F*)GetObjectPtr("CALEAVETOTAL"))->Fill(eTotal);



};


