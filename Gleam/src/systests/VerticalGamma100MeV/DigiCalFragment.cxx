void RootTreeAnalysis::DigiCalHistDefine() {
    
  histFile->cd();

        float maxAdc = 5000.;   // electrons histo max range
	//    float maxAdc = 1000.;   // muons histo max range
    
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

    TH1F *CALNLAYER0 = new TH1F("CALNLAYER0", "Cal Digi Hits per layer 0",
        20, 0, 20);
    CALNLAYER0->SetXTitle("Number");
    TH1F *CALNLAYER1 = new TH1F("CALNLAYER1", "Cal Digi Hits per layer 1",
        20, 0, 20);
    CALNLAYER1->SetXTitle("Number");
    TH1F *CALNLAYER2 = new TH1F("CALNLAYER2", "Cal Digi Hits per layer 2",
        20, 0, 20);
    CALNLAYER2->SetXTitle("Number");
    TH1F *CALNLAYER3 = new TH1F("CALNLAYER3", "Cal Digi Hits per layer 3",
        20, 0, 20);
    CALNLAYER3->SetXTitle("Number");
    TH1F *CALNLAYER4 = new TH1F("CALNLAYER4", "Cal Digi Hits per layer 4",
        20, 0, 20);
    CALNLAYER4->SetXTitle("Number");
    TH1F *CALNLAYER5 = new TH1F("CALNLAYER5", "Cal Digi Hits per layer 5",
        20, 0, 20);
    CALNLAYER5->SetXTitle("Number");
    TH1F *CALNLAYER6 = new TH1F("CALNLAYER6", "Cal Digi Hits per layer 6",
        20, 0, 20);
    CALNLAYER6->SetXTitle("Number");
    TH1F *CALNLAYER7 = new TH1F("CALNLAYER7", "Cal Digi Hits per layer 7",
        20, 0, 20);
    CALNLAYER7->SetXTitle("Number");

    TH1F *CALELAYER0 = new TH1F("CALELAYER0", "Cal Digi ADC summed per layer 0",
        200, 0, maxAdc);
    CALELAYER0->SetXTitle("ADC");
    TH1F *CALELAYER1 = new TH1F("CALELAYER1", "Cal Digi ADC summed per layer 1",
        200, 0, maxAdc);
    CALELAYER1->SetXTitle("ADC");
    TH1F *CALELAYER2 = new TH1F("CALELAYER2", "Cal Digi ADC summed per layer 2",
        200, 0, maxAdc);
    CALELAYER2->SetXTitle("ADC");
    TH1F *CALELAYER3 = new TH1F("CALELAYER3", "Cal Digi ADC summed per layer 3",
        200, 0, maxAdc);
    CALELAYER3->SetXTitle("ADC");
    TH1F *CALELAYER4 = new TH1F("CALELAYER4", "Cal Digi ADC summed per layer 4",
        200, 0, maxAdc);
    CALELAYER4->SetXTitle("ADC");
    TH1F *CALELAYER5 = new TH1F("CALELAYER5", "Cal Digi ADC summed per layer 5",
        200, 0, maxAdc);
    CALELAYER5->SetXTitle("ADC");
    TH1F *CALELAYER6 = new TH1F("CALELAYER6", "Cal Digi ADC summed per layer 6",
        200, 0, maxAdc);
    CALELAYER6->SetXTitle("ADC");
    TH1F *CALELAYER7 = new TH1F("CALELAYER7", "Cal Digi ADC summed per layer 7",
        200, 0, maxAdc);
    CALELAYER7->SetXTitle("ADC");


};

void RootTreeAnalysis::DigiCal() {
    

// get multiplicities

    TObjArray* cL = evt->getCalDigiCol();
    int nCalDigi = cL->GetEntries();
    ((TH1F*)GetObjectPtr("CALDIGICOUNT"))->Fill((Float_t)nCalDigi);

    int nLayer[8]={0,0,0,0,0,0,0,0};
    float eLayer[8]={0.,0.,0.,0.,0.,0.,0.,0.};
    float eTotal = 0.;

    for (int ic=0; ic < nCalDigi; ic++) {
      CalDigi* c=(CalDigi*)cL->At(ic);

      CalXtalId id = c->getPackedId();
      int layer = id->getLayer();

      // get readout from zeroth range. Set to BEST

      CalXtalReadout* cRo=c->getXtalReadout(0);
      float adcP = cRo->getAdc(CalXtalId::POS); 
      float adcN = cRo->getAdc(CalXtalId::NEG);

      // light asymmetry

      float asy = (adcP-adcN)/(adcP+adcN-200.);   // subtract pedestal!

      float eAve = (adcN+adcP)/2.;

      ((TH1F*)GetObjectPtr("CALEAVE"))->Fill(eAve);
      ((TH1F*)GetObjectPtr("CALADC"))->Fill(adcN);
      ((TH1F*)GetObjectPtr("CALADC"))->Fill(adcP);
      ((TH1F*)GetObjectPtr("CALADCN"))->Fill(adcN);
      ((TH1F*)GetObjectPtr("CALADCP"))->Fill(adcP);
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

    ((TH1F*)GetObjectPtr("CALNLAYER0"))->Fill(nLayer[0]);
    ((TH1F*)GetObjectPtr("CALNLAYER1"))->Fill(nLayer[1]);
    ((TH1F*)GetObjectPtr("CALNLAYER2"))->Fill(nLayer[2]);
    ((TH1F*)GetObjectPtr("CALNLAYER3"))->Fill(nLayer[3]);
    ((TH1F*)GetObjectPtr("CALNLAYER4"))->Fill(nLayer[4]);
    ((TH1F*)GetObjectPtr("CALNLAYER5"))->Fill(nLayer[5]);
    ((TH1F*)GetObjectPtr("CALNLAYER6"))->Fill(nLayer[6]);
    ((TH1F*)GetObjectPtr("CALNLAYER7"))->Fill(nLayer[7]);

    ((TH1F*)GetObjectPtr("CALELAYER0"))->Fill(eLayer[0]);
    ((TH1F*)GetObjectPtr("CALELAYER1"))->Fill(eLayer[1]);
    ((TH1F*)GetObjectPtr("CALELAYER2"))->Fill(eLayer[2]);
    ((TH1F*)GetObjectPtr("CALELAYER3"))->Fill(eLayer[3]);
    ((TH1F*)GetObjectPtr("CALELAYER4"))->Fill(eLayer[4]);
    ((TH1F*)GetObjectPtr("CALELAYER5"))->Fill(eLayer[5]);
    ((TH1F*)GetObjectPtr("CALELAYER6"))->Fill(eLayer[6]);
    ((TH1F*)GetObjectPtr("CALELAYER7"))->Fill(eLayer[7]);



};

