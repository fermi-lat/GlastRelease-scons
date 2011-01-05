#ifndef RootTreeAnalysis_cxx
#define RootTreeAnalysis_cxx 1

#include "RootTreeAnalysis.h"

UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum, reconRunNum, mcRunNum;

namespace {
    bool verbose = false;
}

void RootTreeAnalysis::McHistDefine() {
    // Purpose and Method:  Monte Carlo histogram definitions
    
    // Must cd into the histFile to be sure the histograms are created within
    // this histogram file.
    histFile->cd();

    // uncomment the statement if you want to use the pointer
    
    /* TH1F *EVENTMC = */ new TH1F("EVENTMC", "MC Event Id",
        2000, 0, 2000);
    
    /* TH1F *RUNMC = */ new TH1F("RUNMC", "MC Run Id",
        2000, 0, 2000);
    
    /* TH1F *PARTCOUNTMC = */ new TH1F("PARTCOUNTMC", "MC Part Count",
        20, 0, 600);
    /* TH1F *POSCOUNTMC = */ new TH1F("POSCOUNTMC", "MC Pos Count",
        20, 0, 100);
    
    /* TH1F *POSENERGYDEP = */ new TH1F("POSENERGYDEP", "MC Pos Edep",
        25, 0, 0.5);
    
    /* TH1F *POSMCTYPE = */ new TH1F("POSMCTYPE", "MC Pos MC pTypes",
        40, -20, 20);
    
    /* TH1F *INTCOUNTMC = */ new TH1F("INTCOUNTMC", "MC Int Count",
        20, 0, 40);
    
    /* TH1F *INTENERGYDEPCAL = */ new TH1F("INTENERGYDEPCAL", "MC Int Xtal Etot",
        20, 0, 20);
    /* TH1F *INTCALMCTYPE = */ new TH1F("INTCALMCTYPE", "MC Int MC pTypes",
        100, -100, 100);
    
    /* TH1F *INTENERGYDEPACD = */ new TH1F("INTENERGYDEPACD", "MC Int Scint Etot",
        20, 0, 20);
    /* TH1F *INTENERGYTOT = */ new TH1F("INTENERGYTOT", "MC Int Etot",
        20, 0, 500);
    
    /* TH1F *ENERGYMC = */ new TH1F("ENERGYMC", "MC Part Energy (MeV)",
        500, 0, 50000);   
    
}

void RootTreeAnalysis::DigiHistDefine() {
    // Purpose and Method:  Digitization histogram definitions
    // Must cd into the histFile to be sure these new histograms are created in
    // this histogram ROOT file
    histFile->cd();
    
    /* TNtuple *tkrDigiTup= */ new TNtuple("tkrDigiTup", "example Ntuple", "TkrDigiCount:TotalHits");
    /* TNtuple *calDigiTup= */ new TNtuple("calDigiTup", "example Ntuple", "CalDigiCount:CalEAveTotal");

    /* TH1F *NUMTKRDIGI = */ new TH1F("NUMTKRDIGI", "Number of Tkr Digis",
        200, 0, 200);

    /* TH1F *TKRSTRIPSLYR5 = */ new TH1F("TKRSTRIPSLYR5", "Hit Tkr Strips BiLayer 5", 800, 0, 1600);

    /* TH1F *CALDIGICOUNT = */ new TH1F("CALDIGICOUNT", "Cal Digi multiplicity",
        50, 0, 50);
    
    /* TH1F *CALADC = */ new TH1F("CALADC", "Cal Digi ADC - both faces",
        200, 0, 1000);
    /* TH1F *CALRANGE = */ new TH1F("CALRANGE", "Cal Digi Range - both faces",
        10, 0, 10);
    /* TH1F *CALEAVE = */ new TH1F("CALEAVE", "Cal Digi Energy - faces/2",
        200, 0, 1000);
    
    /* TH1F *CALEAVETOTAL = */ new TH1F("CALEAVETOTAL", "Cal Digi Energy - faces/2 summed",
        200, 0, 10000);
    
    /* TH1F *CALLAYER = */ new TH1F("CALLAYER", "Cal Digi Layer",
        15, 0, 15);
    /* TH1F *CALTOWER = */ new TH1F("CALTOWER", "Cal Digi Tower",
        20, 0, 20);
    /* TH1F *CALCOLUMN = */ new TH1F("CALCOLUMN", "Cal Digi Column",
        20, 0, 20);
    
    /* TH1F *CALELYR0 = */ new TH1F("CALELYR0", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR1 = */ new TH1F("CALELYR1", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR2 = */ new TH1F("CALELYR2", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR3 = */ new TH1F("CALELYR3", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR4 = */ new TH1F("CALELYR4", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR5 = */ new TH1F("CALELYR5", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR6 = */ new TH1F("CALELYR6", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    /* TH1F *CALELYR7 = */ new TH1F("CALELYR7", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    
    
    /* TH1F *CALNLYR0 = */ new TH1F("CALNLYR0", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR1 = */ new TH1F("CALNLYR1", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR2 = */ new TH1F("CALNLYR2", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR3 = */ new TH1F("CALNLYR3", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR4 = */ new TH1F("CALNLYR4", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR5 = */ new TH1F("CALNLYR5", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR6 = */ new TH1F("CALNLYR6", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    /* TH1F *CALNLYR7 = */ new TH1F("CALNLYR7", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);

    /* TH1F *ACDPHATILE41 = */ new TH1F ("ACDPHATILE41", "PHA for Tile 041", 
        50, 0, 500);

    /* TH1F *ACDTOTE = */ new TH1F("ACDTOTE", "ACD Total Energy",
        20, 0, 10);
    
    /* TH1F *GEMCONDSUM = */ new TH1F("GEMCONDSUM", "Gem Condition Summary",
        129, -0.5, 128.5);
}

void RootTreeAnalysis::ReconHistDefine() {
    // Purpose and Method:  Reconstruction histogram defintions
    // Must cd into the histFile to be sure these new histograms are created
    // within this histogram ROOT file.
    histFile->cd();

    /* TNtuple *TKRRECONTUP = */ new TNtuple("TKRRECONTUP", "example Ntuple", "TkrNumFitTracks:TkrNumClusters");

    /* TNtuple *CALRECONTUP = */ new TNtuple("CALRECONTUP", "example Ntuple", "NumClusters:NumXtalRec");

    /* TH1F *TKRNUMFITTRACKS = */ new TH1F("TKRNUMFITTRACKS", "Number of Fit Tracks",
        30, 0, 30);

    /* TH1F *TKRNUMHITSPERTRACK = */ new TH1F("TKRNUMHITSPERTRACK", "Number of Hits/Tracks", 40, 0, 40);


    /* TH1F *CALXTALCOUNT = */ new TH1F("CALXTALCOUNT", "Cal Xtal multiplicity",
        50, 0, 200);
    /* TH1F*CALXTALTOTE = */ new TH1F("CALXTALTOTE","Total Xtal Energy in CAL",50,0,20000);
    
    /* TH1F*CALXTALE = */ new TH1F("CALXTALE","Xtal Energy in CAL",50,0,7500);


    /* TH1F*CALNUMCLUS = */ new TH1F("CALNUMCLUS","Number of clusters",20,0,20);
    /* TH1F*CALECLUS = */ new TH1F("CALECLUS","Energy per cluster",50,0,20000);
    /* TH1F*CALTOTE = */ new TH1F("CALTOTE","Total Energy in CAL",50,0,20000);
    /* TH1F*CALEFRAC = */ new TH1F("CALEFRAC","Energy fraction per cluster",50,0,1.);
    TH1F *CALRECESUM = new TH1F("CALRECESUM", "Cal Cluster Energy Sum",
        50, 0, 25000);
    CALRECESUM->SetXTitle("Energy (MeV)");

    TH1F *CALRECELEAK = new TH1F("CALRECELEAK", "Cal Cluster Leakage Energy",
        50, 0, 25000 );
    CALRECELEAK->SetXTitle("Energy (MeV)");

    TH1F *CALRECECORR = new TH1F("CALRECECORR", "Cal Cluster Corrected Energy",
        50, 0, 25000);
    CALRECECORR->SetXTitle("Energy (MeV)");


    /* TH1F *ACDDOCA = */ new TH1F("ACDDOCA", "ACD DOCA",
        40, 0, 200);

    /* TH1F *ACDACTDIST = */ new TH1F("ACDACTDIST", "ACD Active Distance",
        40, 0, 200);

    return;
}


void RootTreeAnalysis::McData() {
    // Purpose and Method:  Process on Monte Carlo event
    
    // Update histograms which keep simple counts
    ((TH1F*)GetObjectPtr("EVENTMC"))->Fill((Float_t)mcEventId);
    ((TH1F*)GetObjectPtr("RUNMC"))->Fill((Float_t)mcRunNum);
    ((TH1F*)GetObjectPtr("PARTCOUNTMC"))->Fill((Float_t)mc->getMcParticleCol()->GetEntries());
    ((TH1F*)GetObjectPtr("POSCOUNTMC"))->Fill((Float_t)mc->getMcPositionHitCol()->GetEntries());
    ((TH1F*)GetObjectPtr("INTCOUNTMC"))->Fill((Float_t)mc->getMcIntegratingHitCol()->GetEntries());
    
    // look at McIntegratingHits
    const TObjArray* iHL = mc->getMcIntegratingHitCol();
    
    int nH = iHL->GetEntries();
    float eTotCal = 0.;
    float eTotAcd = 0.;
    
    for (int ih=0; ih < nH; ih++) {
        McIntegratingHit* hit = (McIntegratingHit*)iHL->At(ih);
        float eHit = hit->getTotalEnergy();
        
        // use volume ID to pick ACD or CAL
        VolumeIdentifier volId = hit->getVolumeId();
        if (volId[0] == 0 &&    // CAL identifier
            volId[3] == 0){ 
            eTotCal += eHit;
            ((TH1F*)GetObjectPtr("INTENERGYDEPCAL"))->Fill(eHit);
            
            
            // examine energy map - note: does not seem to work in CINT right now! Disabled
            
            /*	 const McIntegratingHit::energyDepositMap mcPartMap = hit->getItemizedEnergy();
            McIntegratingHit::energyDepositMap::const_iterator it;
            for (it = mcPartMap.begin(); it != mcPartMap.end(); it++) {
            int mcPartType = (McParticle*(it->first))->getParticleId();
            ((TH1F*)GetObjectPtr("INTCALMCTYPE"))->Fill(mcPartType);
            }
            */
        } else {    // ACD
            eTotAcd += eHit;
            ((TH1F*)GetObjectPtr("INTENERGYDEPACD"))->Fill(eHit);
        }
    }
    ((TH1F*)GetObjectPtr("INTENERGYTOT"))->Fill(eTotCal);
    
    // look at McPositionHit
    const TObjArray* iPL = mc->getMcPositionHitCol();
    int nP = iPL->GetEntries();
    
    for (int ip=0; ip < nP; ip++) {
        McPositionHit* pHit = (McPositionHit*)iPL->At(ip);
        ((TH1F*)GetObjectPtr("POSENERGYDEP"))->Fill((Float_t)pHit->getDepositedEnergy());
        
        // get particle type that caused this hit. 
        Int_t mcType = pHit->getMcParticleId();
        ((TH1F*)GetObjectPtr("POSMCTYPE"))->Fill(mcType);
    }
}

void RootTreeAnalysis::DigiTkr() {
    // Purpose and Method: Process on TKR digi event
    
    // The full collection of TkrDigis for this event
    const TObjArray* tkrDigiCol = evt->getTkrDigiCol();
    if (!tkrDigiCol) return;

    // Total number of TkrDigis for this event
    Int_t numTkrDigi = tkrDigiCol->GetEntries();
    Int_t totalHits = 0;

    // Loop over all TkrDigis
    TIter tkrIter(tkrDigiCol);
    TkrDigi *tkr = 0;
    while ((tkr = (TkrDigi*)tkrIter.Next())) {
      // Identify the tower and layer
      Int_t tower = tkr->getTower().id();
      Int_t layer = tkr->getBilayer();

      // get to ToT
      Int_t tot0 = tkr->getToT(0);
      Int_t tot1 = tkr->getToT(1);
     
      Int_t lastController0 = tkr->getLastController0Strip();

      // Returns the orientation of the strips
      GlastAxis::axis view = tkr->getView();

      UInt_t numHits = tkr->getNumHits();
      totalHits += numHits;
      // Loop through collection of hit strips for this TkrDigi
      UInt_t ihit;
      for (ihit = 0; ihit < numHits; ihit++) {
        // Retrieve the strip number
          Int_t stripNum = tkr->getStrip(ihit);
          // this is mainly to use the local vars
          if (verbose) {
              std::cout << "DigiTkr " << tower << " " << tot0 <<  " " << tot1
                  << " " << lastController0 << " " << view << std::endl;
          }
          if (layer == 5) ((TH1F*)GetObjectPtr("TKRSTRIPSLYR5"))->Fill(stripNum);
      }
    }

    ((TH1F*)GetObjectPtr("NUMTKRDIGI"))->Fill((Float_t)numTkrDigi);

    // example of filling a simple tree
    Float_t digiArr[2] = { (Float_t)numTkrDigi, (Float_t)totalHits };
    ((TNtuple*)GetObjectPtr("tkrDigiTup"))->Fill(digiArr);

    return;
}

void RootTreeAnalysis::DigiCal() {
    // Purpose and Method:  Process on CAL digi event
        
    const TObjArray* calDigiCol = evt->getCalDigiCol();
    if (!calDigiCol) return;

    Int_t nCalDigi = calDigiCol->GetEntries();
    ((TH1F*)GetObjectPtr("CALDIGICOUNT"))->Fill((Float_t)nCalDigi);
    
    Int_t nLayer[8]={0,0,0,0,0,0,0,0};
    Float_t eLayer[8]={0.,0.,0.,0.,0.,0.,0.,0.};
    Float_t eTotal = 0.;

    TIter calDigiIter(calDigiCol);
    CalDigi *c = 0;
    
    while ((c = (CalDigi*)calDigiIter.Next())) {

        // Assuming Best Range and checking only the first readout
        const CalXtalReadout* cRo=c->getXtalReadout(0);
        Float_t adcP = cRo->getAdc(CalXtalId::POS);
        Float_t adcN = cRo->getAdc(CalXtalId::NEG);

        // ligh asymmetry
        Float_t asy = 1.;
        if (adcP + adcN == 200.)
          asy = 1.;
        else 
          asy = (adcP-adcN)/(adcP+adcN-200.);

        Float_t eAve = (adcP + adcN)/2.;

        ((TH1F*)GetObjectPtr("CALEAVE"))->Fill(eAve);
        ((TH1F*)GetObjectPtr("CALADC"))->Fill((float)cRo->getAdc(CalXtalId::POS));
        ((TH1F*)GetObjectPtr("CALADC"))->Fill((float)cRo->getAdc(CalXtalId::NEG));
        ((TH1F*)GetObjectPtr("CALRANGE"))->Fill(cRo->getRange(CalXtalId::POS));
        ((TH1F*)GetObjectPtr("CALRANGE"))->Fill(cRo->getRange(CalXtalId::NEG));

        // Retrieve the identifer for this crystal
        CalXtalId id = c->getPackedId();
        Int_t layer = id.getLayer();
        Int_t tower = id.getTower();
        Int_t column = id.getColumn();
        ((TH1F*)GetObjectPtr("CALLAYER"))->Fill(layer);
        ((TH1F*)GetObjectPtr("CALTOWER"))->Fill(tower);
        ((TH1F*)GetObjectPtr("CALCOLUMN"))->Fill(column);
        
        nLayer[layer] += 1;
        eLayer[layer] += eAve;
        eTotal += eAve;
    }
    
    ((TH1F*)GetObjectPtr("CALEAVETOTAL"))->Fill(eTotal);
    
    ((TH1F*)GetObjectPtr("CALELYR0"))->Fill(eLayer[0]);
    ((TH1F*)GetObjectPtr("CALELYR1"))->Fill(eLayer[1]);
    ((TH1F*)GetObjectPtr("CALELYR2"))->Fill(eLayer[2]);
    ((TH1F*)GetObjectPtr("CALELYR3"))->Fill(eLayer[3]);
    ((TH1F*)GetObjectPtr("CALELYR4"))->Fill(eLayer[4]);
    ((TH1F*)GetObjectPtr("CALELYR5"))->Fill(eLayer[5]);
    ((TH1F*)GetObjectPtr("CALELYR6"))->Fill(eLayer[6]);
    ((TH1F*)GetObjectPtr("CALELYR7"))->Fill(eLayer[7]);
    
    ((TH1F*)GetObjectPtr("CALNLYR0"))->Fill(nLayer[0]);
    ((TH1F*)GetObjectPtr("CALNLYR1"))->Fill(nLayer[1]);
    ((TH1F*)GetObjectPtr("CALNLYR2"))->Fill(nLayer[2]);
    ((TH1F*)GetObjectPtr("CALNLYR3"))->Fill(nLayer[3]);
    ((TH1F*)GetObjectPtr("CALNLYR4"))->Fill(nLayer[4]);
    ((TH1F*)GetObjectPtr("CALNLYR5"))->Fill(nLayer[5]);
    ((TH1F*)GetObjectPtr("CALNLYR6"))->Fill(nLayer[6]);
    ((TH1F*)GetObjectPtr("CALNLYR7"))->Fill(nLayer[7]);
    
    // example of filling a simple tree
    Float_t digiArr[2] = { (Float_t)nCalDigi, eTotal };
    ((TNtuple*)GetObjectPtr("calDigiTup"))->Fill(digiArr);
}

void RootTreeAnalysis::DigiAcd() {
    // Purpose and Method:  Process on ACD digi event
    //   Determines the total energy deposited and
    //   store the PHA for tile 0000 in a histogram

    const TObjArray* acdDigiCol = evt->getAcdDigiCol();
    if (!acdDigiCol) return;

    Double_t totE = 0.0;
    UShort_t pha0 = 0;
    // Create an id for tile 0041
    AcdId id41(0, 0, 4, 1);

    TIter acdDigiIter(acdDigiCol);
    AcdDigi *acdDigiItem = 0;
    
    while ((acdDigiItem = (AcdDigi*)acdDigiIter.Next())) {
        AcdId id = acdDigiItem->getId();
        if (id.getId() == id41.getId()) {
            pha0 = acdDigiItem->getPulseHeight(AcdDigi::A) + 
                acdDigiItem->getPulseHeight(AcdDigi::B);
        }
    }

    ((TH1F*)GetObjectPtr("ACDPHATILE41"))->Fill(pha0);
    ((TH1F*)GetObjectPtr("ACDTOTE"))->Fill(totE);

    return;
}

void RootTreeAnalysis::DigiGem() {
    const Gem& gem = evt->getGem();

    ((TH1F*)GetObjectPtr("GEMCONDSUM"))->Fill((Float_t)gem.getConditionSummary());

    UShort_t deadZone =  gem.getMissed();
  
    UInt_t liveTime = gem.getLiveTime();

    const GemCondArrivalTime& condArrTime = gem.getCondArrTime();
    if(verbose) {
        std::cout << "DigiGem " << deadZone << " " << liveTime << " " 
            << condArrTime.external() << std::endl;
    }
}

void RootTreeAnalysis::DigiDiagnostic() {

    const TClonesArray *calDiagCol = evt->getCalDiagnosticCol();
    const TClonesArray *tkrDiagCol = evt->getTkrDiagnosticCol();

    if (!calDiagCol && !tkrDiagCol) return;

    if (calDiagCol) {

      TIter calDiagIt(calDiagCol);
      CalDiagnosticData *cal = 0;
      while( (cal = (CalDiagnosticData*)calDiagIt.Next()) ) {

      }

    }

    if (tkrDiagCol) {
      TIter tkrDiagIt(tkrDiagCol);
      TkrDiagnosticData *tkr = 0;
      while( (tkr = (TkrDiagnosticData*)tkrDiagIt.Next()) ) {

      }
   
    }

}

void RootTreeAnalysis::ReconTkr() {
    // Purpose and Method:  Process one TkrRecon event

    TkrRecon *tkrRecon = rec->getTkrRecon();
    
    // If no TRKRECON data is available then return
    if (!tkrRecon)  return;

    // Retrieve collection of clusters
    TObjArray* clusterCol = tkrRecon->getClusterCol();

    // Retreive collection of vertices
    TObjArray* vertexCol = tkrRecon->getVertexCol();

    // Retrieve collection of tracks    
    TObjArray* trackCol = tkrRecon->getTrackCol();

     if (trackCol) {
        ((TH1F*)GetObjectPtr("TKRNUMFITTRACKS"))->Fill(trackCol->GetEntries());

        // loop over all tracks
        TIter trackIter(trackCol);
        TkrTrack *track = 0;
        while ((track = (TkrTrack*)trackIter.Next())) {
           ((TH1F*)GetObjectPtr("TKRNUMHITSPERTRACK"))->Fill(track->Size());

            TVector3 initPos = track->getInitialPosition();
            TVector3 dir = track->getInitialDirection();
            Double_t costh = -dir.Z();

            Double_t quality = track->getQuality();
            Double_t kalEnergy = track->getKalEnergy();

            if(verbose) {
                std::cout << costh << " " << quality << " " << kalEnergy << std::endl;
            }

            // loop over hits in this track
            TkrTrackHit *hit = 0;
            UInt_t i;
            for (i=0; i < track->Size(); i++ ) {
                hit = (TkrTrackHit*)track->At(i);
                Double_t planeZ = hit->getZPlane();

                const TkrCluster *cluster = ((*hit).getClusterPtr());
                if (!cluster) continue;

                Int_t layer = cluster->getLayer();

                commonRootData::TkrId id = cluster->getTkrId();
                Int_t plane  = cluster->getPlane();
                Int_t view   = id.getView();

                // get the tower from the TkrId
                Int_t towerX = id.getTowerX();
                Int_t towerY = id.getTowerY();
                if (verbose) {
                    std::cout << planeZ << " " << layer << " " << plane << " "
                        << view << " " << towerX << " " << towerY << " " << std::endl;
                }

           }

        }
    }

    if (trackCol && clusterCol) {
        Float_t recArr[2] = {(Float_t)trackCol->GetEntries(), (Float_t)clusterCol->GetEntries()};
        ((TNtuple*)GetObjectPtr("TKRRECONTUP"))->Fill(recArr);
    }
}

void RootTreeAnalysis::ReconCal() {
    // Purpose and Method:  Process on CalRecon event
    
    CalRecon *calRec = rec->getCalRecon();
    if (!calRec) return;

    float totXE = 0.;
    
    TObjArray *xtalRecCol = calRec->getCalXtalRecCol();
    if (xtalRecCol) {
        TIter xtalIter(xtalRecCol);
        CalXtalRecData *xtal = 0;
        while ((xtal = (CalXtalRecData*)xtalIter.Next())) {
            Double_t xtalEnergy = xtal->getEnergy();
            if (xtalEnergy > 2000) {
                const CalXtalId id = xtal->getPackedId();
	        int lyr = id.getLayer();
	        int twr = id.getTower();
	        int col = id.getColumn();
	        const CalRangeRecData* rData = xtal->getRangeRecData(0);
	        CalXtalId::AdcRange range = (CalXtalId::AdcRange) (rData->getRange(CalXtalId::POS));
	        double ph0 = xtal->getEnergySelectedRange(range,CalXtalId::POS);
            if (verbose) {
                std::cout << "ReconCal " << lyr << " " << twr << " " << col << " " 
                    << ph0 << std::endl;
            }
  	        continue;
             }
          totXE += xtalEnergy;

        }


        ((TH1F*)GetObjectPtr("CALXTALTOTE"))->Fill(totXE);

    }
    

    TObjArray*  clusCol = calRec->getCalClusterCol();
    if (clusCol) {
        Int_t numClus = clusCol->GetEntries();
        ((TH1F*)GetObjectPtr("CALNUMCLUS"))->Fill(numClus);

        float totE = 0.;
        for (int jc=0;jc<numClus; jc++) {
            CalCluster* c1 = (CalCluster*)clusCol->At(jc);
            float clusterEnergy = c1->getMomParams().getEnergy();
            totE += clusterEnergy;
            ((TH1F*)GetObjectPtr("CALRECESUM"))->Fill(clusterEnergy);
        }
        ((TH1F*)GetObjectPtr("CALTOTE"))->Fill(totE);

        if (totE > 0.) {
            for (int ic=0;ic<numClus; ic++) {
	        CalCluster* c2 = (CalCluster*)clusCol->At(ic);
	        float clusterEnergy = c2->getMomParams().getEnergy();
	        ((TH1F*)GetObjectPtr("CALECLUS"))->Fill(clusterEnergy);
                float eFraction = clusterEnergy/totE;
                ((TH1F*)GetObjectPtr("CALEFRAC"))->Fill(eFraction);
            }
        }
    }



    if (clusCol && xtalRecCol) {
        Float_t recArr[2] = {clusCol->GetEntries(), xtalRecCol->GetEntries()};
        ((TNtuple*)GetObjectPtr("CALRECONTUP"))->Fill(recArr);
    }

    return;
}

void RootTreeAnalysis::ReconAcd() {
    // Purpose and Method:  Processes ACD reconstruction data

    AcdRecon *acdRec = rec->getAcdRecon();
    if (!acdRec) return;

    ((TH1F*)GetObjectPtr("ACDDOCA"))->Fill(acdRec->getDoca());
    ((TH1F*)GetObjectPtr("ACDACTDIST"))->Fill(acdRec->getActiveDist());

    return;
}

void RootTreeAnalysis::HistDefine() {
    // Purpose and Method:  Setup Histograms
    
    gStyle->SetOptStat(111111);
    
    histFile = new TFile(m_histFileName,"RECREATE");
    
    McHistDefine();
    DigiHistDefine();
    ReconHistDefine();
    
}

void RootTreeAnalysis::Go(Long64_t numEvents)
{    
    // Purpose and Method:  Event Loop
    //   All analysis goes here
    
    //  To read only selected branches - saves processing time
    //  Comment out any branches you are not interested in.
    
    // mc branches:
    if (mcTree) {
        mcTree->SetBranchStatus("*", 1);    // enable all branches
        /*  Use the following if you want to save I/O and only read a portion
            of the data
        mcTree->SetBranchStatus("*", 0);    // disable all branches
        // Activate desired branches...
        mcTree->SetBranchStatus("m_eventId", 1);
        mcTree->SetBranchStatus("m_particleCol", 1);
        mcTree->SetBranchStatus("m_runId", 1);        
        mcTree->SetBranchStatus("m_integratingHitCol", 1);        
        mcTree->SetBranchStatus("m_positionHitCol", 1);        
       */
    }
    
    if (digiTree) {
        digiTree->SetBranchStatus("*",1);  // enable all branches
        /* Use the following if you want to save I/O and read only a portion of the data
        digiTree->SetBranchStatus("*",0);  // disable all branches
        // activate desired brances
        digiTree->SetBranchStatus("m_cal*",1);  
        digiTree->SetBranchStatus("m_tkr*",1);  
        digiTree->SetBranchStatus("m_acd*",1);
        digiTree->SetBranchStatus("m_eventId", 1); 
        digiTree->SetBranchStatus("m_runId", 1);
        digiTree->SetBranchStatus("m_gem", 1);
        */
    }
    
    if (reconTree) {
        reconTree->SetBranchStatus("*",1);  // enable all branches
        /* Use the following if you want to save I/O and read only a portion of the data
        reconTree->SetBranchStatus("*",0);  // disable all branches
        // activate desired branches
        reconTree->SetBranchStatus("m_cal*", 1);  
        reconTree->SetBranchStatus("m_tkr*", 1);
        reconTree->SetBranchStatus("m_acd*", 1);
        reconTree->SetBranchStatus("m_eventId", 1); 
        reconTree->SetBranchStatus("m_runId", 1);
        */
    }
        
    // determine how many events to process
    Long64_t nentries = GetEntries();
    std::cout << "\nNum Events in File is: " << nentries << std::endl;
    Int_t curI;
    Int_t nMax = TMath::Min(numEvents+m_StartEvent,nentries);
    if (m_StartEvent == nentries) {
        std::cout << " all events in file read" << std::endl;
        return;
    }
    if (nentries <= 0) return;
    
    // Keep track of how many bytes we have read in from the data files
    Int_t nbytes = 0, nb = 0;
    
    // BEGINNING OF EVENT LOOP
    for (Long64_t ievent=m_StartEvent; ievent<nMax; ievent++, curI=ievent) {
        
        if (mc) mc->Clear();
        if (evt) evt->Clear();
        if (rec) rec->Clear();
        
        digiEventId = 0; reconEventId = 0; mcEventId = 0;
        digiRunNum = 0; reconRunNum = 0; mcRunNum = 0;
        
        nb = GetEvent(ievent);
        nbytes += nb;
        
        // Monte Carlo ONLY analysis
        if (mc) {  // if we have mc data process it
            mcEventId = mc->getEventId();
            mcRunNum = mc->getRunId();
            McData();
        } 
        
        // Digi ONLY analysis
        if (evt) {
            digiEventId = evt->getEventId(); 
            digiRunNum = evt->getRunId();
            
            DigiTkr();
            DigiCal();
            DigiAcd();

            //DigiGem();
            //DigiDiagnostic();
        }
        
        
        // RECON ONLY analysis
        if (rec) {  // if we have recon data proccess it            
            ReconTkr();
            ReconCal();
            ReconAcd();
        } 
        
        
    }  // end analysis code in event loop
    
    m_StartEvent = curI;
}

#endif
