void RootTreeAnalysis::ReconCalHistDefine() {
    
    histFile->cd();
    
    float maxEnergy = 50.; //MeV

    TH1F *CALCLUSCOUNT = new TH1F("CALCLUSCOUNT", "Cal Cluster Count",
        10, 0, 10);
    CALCLUSCOUNT->SetXTitle("Number");

    TH1F *CALRMSTRANS = new TH1F("CALRMSTRANS", "Cal Cluster Transverse RMS",
        50, 0, 250000);
    CALRMSTRANS->SetXTitle("units??");

    TH1F *CALRMSLONG = new TH1F("CALRMSLONG", "Cal Cluster Longitudinal RMS",
        50, 0, 40000);
    CALRMSLONG->SetXTitle("units??");

    TH1F *CALTRANSVOFFSET = new TH1F("CALTRANSVOFFSET", "Cal Cluster Transverse Offset",
        50, 0, 250);
    CALTRANSVOFFSET->SetXTitle("units??");

    TH1F *CALXPOS = new TH1F("CALXPOS", "Cal Cluster x position",
        50, 50, 250);
    CALXPOS->SetXTitle("position (mm)");
    TH1F *CALYPOS = new TH1F("CALYPOS", "Cal Cluster x position",
        50, 50, 250);
    CALYPOS->SetXTitle("position (mm)");
    TH1F *CALZPOS = new TH1F("CALZPOS", "Cal Cluster x position",
        50, -200, 0);
    CALZPOS->SetXTitle("position (mm)");

    TH1F *CALDRC = new TH1F("CALDRC", "Cal Cluster direction cosine",
        50, -1, 1);
    CALDRC->SetXTitle("direction cosine");



    TH1F *CALRECESUM = new TH1F("CALRECESUM", "Cal Cluster Energy Sum",
        50, 0, maxEnergy);
    CALRECESUM->SetXTitle("Energy (MeV)");

    TH1F *CALRECELEAK = new TH1F("CALRECELEAK", "Cal Cluster Leakage Energy",
        50, 0, maxEnergy);
    CALRECELEAK->SetXTitle("Energy (MeV)");

    TH1F *CALRECECORR = new TH1F("CALRECECORR", "Cal Cluster Corrected Energy",
        50, 0, maxEnergy);
    CALRECECORR->SetXTitle("Energy (MeV)");

    TH1F *CALRECEFIT = new TH1F("CALRECEFIT", "Cal Cluster Fitted Energy",
        50, 0, maxEnergy);
    CALRECEFIT->SetXTitle("Energy (MeV)");

}

void RootTreeAnalysis::ReconCal() {
    
  CalRecon* cRH = rec->getCalRecon();

// get multiplicities

    TObjArray* clusterCol = cRH->getCalClusterCol();
    int nCalClus = clusterCol->GetEntries();
    ((TH1F*)GetObjectPtr("CALCLUSCOUNT"))->Fill((Float_t)nCalClus);

    for (int rc=0; rc < nCalClus; rc++) {
      CalCluster* c=(CalCluster*)clusterCol->At(rc);
      ((TH1F*)GetObjectPtr("CALRECESUM"))->Fill(c->getEnergySum());
      ((TH1F*)GetObjectPtr("CALRECELEAK"))->Fill(c->getEnergyLeak());
      ((TH1F*)GetObjectPtr("CALRECECORR"))->Fill(c->getEnergyCorrected());
      ((TH1F*)GetObjectPtr("CALRECEFIT"))->Fill(c->getFitEnergy());
      ((TH1F*)GetObjectPtr("CALRMSTRANS"))->Fill(c->getRmsTrans());
      ((TH1F*)GetObjectPtr("CALRMSLONG"))->Fill(c->getRmsLong());
      ((TH1F*)GetObjectPtr("CALTRANSVOFFSET"))->Fill(c->getTransvOffset());
      ((TH1F*)GetObjectPtr("CALXPOS"))->Fill(c->getPosition().Px());
      ((TH1F*)GetObjectPtr("CALYPOS"))->Fill(c->getPosition().Py());
      ((TH1F*)GetObjectPtr("CALZPOS"))->Fill(c->getPosition().Pz());
      ((TH1F*)GetObjectPtr("CALDRC"))->Fill(c->getDirection().CosTheta());
    }
};
