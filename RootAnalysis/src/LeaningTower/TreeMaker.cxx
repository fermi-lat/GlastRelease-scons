#ifndef TreeMaker_cxx
#define TreeMaker_cxx 1

#include "TreeMaker.h"
#include <vector>

TObjArray tkrLayerHistArr;
UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum, reconRunNum, mcRunNum;

const long int TOWERNUM = 1;
const long int LAYERNUM = 5;
const long int VIEWNUM  = 2;
const TString TS = "";
#define DEBUG 0

void makeIndices(long int t, long int l, long int v, TString idx[3]) 
{
  std::cout << "a" << std::endl;
  idx[2] = TS + "[" + t + "][" + l + "][" + v + "]";
  std::cout << "b" << std::endl;
  idx[1] = idx[2] + "[1]";
  std::cout << "c" << std::endl;
  idx[0] = idx[2] + "[0]";
  std::cout << "d" << std::endl;
}

void TreeMaker::McHistDefine() {
    // Purpose and Method:  Monte Carlo histogram definitions
    
	// Must cd into the histFile to be sure the histograms are created within
	// this histogram histFile.
    histFile->cd();
    
    TH1F *EVENTMC = new TH1F("EVENTMC", "MC Event Id",        2000, 0, 2000);
    
    TH1F *RUNMC = new TH1F("RUNMC", "MC Run Id",        2000, 0, 2000);
    
    TH1F *PARTCOUNTMC = new TH1F("PARTCOUNTMC", "MC Part Count",        20, 0, 600);
    TH1F *POSCOUNTMC = new TH1F("POSCOUNTMC", "MC Pos Count",        20, 0, 100);
    
    TH1F *POSENERGYDEP = new TH1F("POSENERGYDEP", "MC Pos Edep",        25, 0, 0.5);
    
    TH1F *POSMCTYPE = new TH1F("POSMCTYPE", "MC Pos MC pTypes",        40, -20, 20);
    
    TH1F *INTCOUNTMC = new TH1F("INTCOUNTMC", "MC Int Count",        20, 0, 40);
    
    TH1F *INTENERGYDEPCAL = new TH1F("INTENERGYDEPCAL", "MC Int Xtal Etot",        20, 0, 20);
    TH1F *INTCALMCTYPE = new TH1F("INTCALMCTYPE", "MC Int MC pTypes",        100, -100, 100);
    
    TH1F *INTENERGYDEPACD = new TH1F("INTENERGYDEPACD", "MC Int Scint Etot",        20, 0, 20);
    TH1F *INTENERGYTOT = new TH1F("INTENERGYTOT", "MC Int Etot",        20, 0, 500);
    
    TH1F *ENERGYMC = new TH1F("ENERGYMC", "MC Part Energy (MeV)",        500, 0, 50000);   
    
}

void TreeMaker::TkrDigiHistDefine() {
    // Purpose and Method:  Tkr Digitization histogram definitions
    histFile->cd();

    TH1F* TkrDigiCount = new TH1F("TkrDigiCount", "TkrDigi multiplicity;# of TkrDigi per event", 128, -0.5, 127.5);
    TH1F* TkrDigiToT[TOWERNUM][LAYERNUM][VIEWNUM][2];
    TH1F* TkrDigiHitCount[TOWERNUM][LAYERNUM][VIEWNUM][3];
    TH1F* TkrDigiHitMap[TOWERNUM][LAYERNUM][VIEWNUM];

    for ( long int t=0; t<TOWERNUM; ++t )
        for ( long int l=0; l<LAYERNUM; ++l )
            for ( long int v=0; v<VIEWNUM; ++v ) {
                TString idx[3];
                //                makeIndices(t, l, v, idx);
                idx[2] = TS + "[" + t + "][" + l + "][" +  v + "]";
                idx[1] = idx[2] + "[1]";
                idx[0] = idx[2] + "[0]";
                TkrDigiHitCount[t][l][v][2] = new TH1F("TkrDigiHitCount"+idx[2], "TkrDigi hit counter"+idx[2]+";# of Tkr hits per plane", 256, -0.5, 255.5);
                TkrDigiHitMap[t][l][v] = new TH1F("TkrDigiHitMap"+idx[2], "TkrDigi hit map"+idx[2]+";strip id", 1536, -0.5, 1535.5);
                for ( long int roc=0; roc<2; ++roc ) {
		  TkrDigiToT[t][l][v][roc] = new TH1F("TkrDigiToT"+idx[roc], "TkrDigi ToT"+idx[roc]+";ADC channels", 256, -0.5, 255.5);
		  TkrDigiHitCount[t][l][v][roc] = new TH1F("TkrDigiHitCount"+idx[roc], "TkrDigi hit counter"+idx[roc]+";# of Tkr hits per ROC", 256, -0.5, 255.5);
                }
            }
}

void TreeMaker::CalDigiHistDefine() {
    // Purpose and Method:  Cal Digitization histogram definitions
    histFile->cd();
    TH1F *CALDIGICOUNT = new TH1F("CALDIGICOUNT", "Cal Digi multiplicity",        50, 0, 50);
    
    TH1F *CALADC = new TH1F("CALADC", "Cal Digi ADC - both faces",        200, 0, 1000);
    TH1F *CALRANGE = new TH1F("CALRANGE", "Cal Digi Range - both faces",        10, 0, 10);
    TH1F *CALEAVE = new TH1F("CALEAVE", "Cal Digi Energy - faces/2",        200, 0, 1000);
    
    TH1F *CALEAVETOTAL = new TH1F("CALEAVETOTAL", "Cal Digi Energy - faces/2 summed",        200, 0, 10000);
    
    TH1F *CALLAYER = new TH1F("CALLAYER", "Cal Digi Layer",        15, 0, 15);
    TH1F *CALTOWER = new TH1F("CALTOWER", "Cal Digi Tower",        20, 0, 20);
    TH1F *CALCOLUMN = new TH1F("CALCOLUMN", "Cal Digi Column",        20, 0, 20);
    
    TH1F *CALELYR0 = new TH1F("CALELYR0", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR1 = new TH1F("CALELYR1", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR2 = new TH1F("CALELYR2", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR3 = new TH1F("CALELYR3", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR4 = new TH1F("CALELYR4", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR5 = new TH1F("CALELYR5", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR6 = new TH1F("CALELYR6", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    TH1F *CALELYR7 = new TH1F("CALELYR7", "Cal Digi Energy - faces/2 - Layer 0",        200, 0, 1000);
    
    
    TH1F *CALNLYR0 = new TH1F("CALNLYR0", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR1 = new TH1F("CALNLYR1", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR2 = new TH1F("CALNLYR2", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR3 = new TH1F("CALNLYR3", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR4 = new TH1F("CALNLYR4", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR5 = new TH1F("CALNLYR5", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR6 = new TH1F("CALNLYR6", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);
    TH1F *CALNLYR7 = new TH1F("CALNLYR7", "Cal Digi Energy - faces/2 - Layer 0",        10, 0, 10);

    TH1F *ACDPHATILE41 = new TH1F ("ACDPHATILE41", "PHA for Tile 041",         50, 0, 500);

    TH1F *ACDTOTE = new TH1F("ACDTOTE", "ACD Total Energy",        20, 0, 10);
}

void TreeMaker::AcdDigiHistDefine() {
    // Purpose and Method:  AcdDigitization histogram definitions
    histFile->cd();
}

void TreeMaker::DigiHistDefine() {
    // Purpose and Method:  Digitization histogram definitions
    histFile->cd();
    TH1F* DigiEventID = new TH1F("DigiEventID", "Digi Event Id;event id", 2000, -0.5, 1999.5);
    TH1F* DigiRunID   = new TH1F("DigiRunID",   "Digi Run Id;run number",   2000, 2120000.5, 2122000.5);

    TkrDigiHistDefine();
    //    CalDigiHistDefine();
    //    AcdDigiHistDefine();
}

void TreeMaker::ReconHistDefine() {
    // Purpose and Method:  Reconstruction histogram defintions
    
	// Must cd into the histFile to be sure these new histograms are created
	// within this histogram ROOT file.
    histFile->cd();

    TH1F *TKRNUMFITTRACKS = new TH1F("TKRNUMFITTRACKS", "Number of Fit Tracks",        30, 0, 30);

    TH1F *TKRNUMHITSPERTRACK = new TH1F("TKRNUMHITSPERTRACK", "Number of Hits/Tracks",        40, 0, 40);

    TH1F *ACDDOCA = new TH1F("ACDDOCA", "ACD DOCA",        40, 0, 200);

    TH1F *ACDACTDIST = new TH1F("ACDACTDIST", "ACD Active Distance",        40, 0, 200);

    return;
}


void TreeMaker::McData() {
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

void TreeMaker::DigiTkr() {
    // Purpose and Method: Process on TKR digi event

  const TObjArray* tkrDigiCol = 0;
    tkrDigiCol = evt->getTkrDigiCol();
    if (!tkrDigiCol)
        return;

    const int nTkrDigi = tkrDigiCol->GetEntries();
    if ( DEBUG ) std::cout << "num tkrDigi " << nTkrDigi << std::endl;
    ((TH1F*)GetObjectPtr("TkrDigiCount"))->Fill(nTkrDigi);

    TIter tkrDigiIter(tkrDigiCol);
    TkrDigi* pTkrDigi = 0;

    while ( pTkrDigi = (TkrDigi*)tkrDigiIter.Next() ) {

        const long int t = pTkrDigi->getTower().id();
        const long int l = pTkrDigi->getBilayer();
        const long int v = pTkrDigi->getView();
        if ( DEBUG ) std::cout << "tower, bilayer, view: " << t << " " << l << " " << v << std::endl;
        if ( t>=TOWERNUM ) {
            std::cout << "no histograms defined for tower " << t << std::endl;
            exit(0);
        }
        if ( l>=LAYERNUM ) {
            std::cout << "no histograms defined for layer " << l << std::endl;
            exit(0);
        }
        if ( v>=VIEWNUM ) {
            std::cout << "no histograms defined for view " << v << std::endl;
            exit(0);
        }

        TString idx[3];
        idx[2] = TS + "[" + t + "][" + l + "][" + v + "]";
        idx[1] = idx[2] + "[1]";
        idx[0] = idx[2] + "[0]";
        
        const int ToT[2] = { pTkrDigi->getToT(0), pTkrDigi->getToT(1) };
        if ( DEBUG ) std::cout << "ToT: " << ToT[0] << " " << ToT[1] << std::endl;
        ((TH1F*)GetObjectPtr("TkrDigiToT"+idx[0]))->Fill(ToT[0]);
        ((TH1F*)GetObjectPtr("TkrDigiToT"+idx[1]))->Fill(ToT[1]);
	
        int last0strip = pTkrDigi->getLastController0Strip();
        if ( DEBUG ) std::cout << "last strip of controller 0: " << last0strip << std::endl;
        
        // I would have liked to use the hitCol, but I cannot use the size() method from within root:
        // Error: Can't call vector<unsigned int,allocator<unsigned int> >::size() in current scope ...
        //        std::vector<UInt_t> hitCol = pTkrDigi->getHitCol(); 
        //        const int nHitCol = hitCol.size();
        const int nHitCol = pTkrDigi->getNumHits();
        if ( DEBUG ) std::cout << "num hits: " << nHitCol << std::endl;
        int nHits[3] = { 0, 0, 0 };
        if ( nHitCol > 0 ) 
	  {
            if ( DEBUG ) std::cout << "list of hits:";
            //            for ( std::vector<UInt_t>::iterator it = hitCol.begin(); it<hitCol.end(); ++it )
            //                if ( DEBUG ) std::cout << " " << it;
	    
	    for ( int i=0; i<nHitCol; ++i ) {
	      int hit = pTkrDigi->getHit(i);
	      if ( DEBUG ) std::cout << " " << pTkrDigi->getHit(i);
	      //                ++nHits[2];  trivial, would be a check only if we could use iterators
	      ((TH1F*)GetObjectPtr("TkrDigiHitMap"+idx[2]))->Fill(hit);
	      if ( hit > last0strip )
		++nHits[1];
	      else
		++nHits[0];
	    }
	    
            if ( DEBUG ) std::cout << std::endl;
	  }
	
        ((TH1F*)GetObjectPtr("TkrDigiHitCount"+idx[2]))->Fill(nHitCol);
        ((TH1F*)GetObjectPtr("TkrDigiHitCount"+idx[1]))->Fill(nHits[1]);
        ((TH1F*)GetObjectPtr("TkrDigiHitCount"+idx[0]))->Fill(nHits[0]);
    }
    
    return;
}

void TreeMaker::DigiCal() {
    // Purpose and Method:  Process on CAL digi event
        
    const TObjArray* calDigiCol = evt->getCalDigiCol();
    if (!calDigiCol) return;

    int nCalDigi = calDigiCol->GetEntries();
    ((TH1F*)GetObjectPtr("CALDIGICOUNT"))->Fill((Float_t)nCalDigi);
    
    int nLayer[8]={0,0,0,0,0,0,0,0};
    float eLayer[8]={0.,0.,0.,0.,0.,0.,0.,0.};
    float eTotal = 0.;

    TIter calDigiIter(calDigiCol);
    CalDigi *c = 0;
    
    while (c = (CalDigi*)calDigiIter.Next()) {
        const CalXtalReadout* cRo=c->getXtalReadout(0);
        float eAve = (cRo->getAdc(CalXtalId::POS)+cRo->getAdc(CalXtalId::NEG))/2.;
        ((TH1F*)GetObjectPtr("CALEAVE"))->Fill(eAve);
        ((TH1F*)GetObjectPtr("CALADC"))->Fill((float)cRo->getAdc(CalXtalId::POS));
        ((TH1F*)GetObjectPtr("CALADC"))->Fill((float)cRo->getAdc(CalXtalId::NEG));
        ((TH1F*)GetObjectPtr("CALRANGE"))->Fill(cRo->getRange(CalXtalId::POS));
        ((TH1F*)GetObjectPtr("CALRANGE"))->Fill(cRo->getRange(CalXtalId::NEG));
        CalXtalId id = c->getPackedId();
        int layer = id.getLayer();
        int tower = id.getTower();
        int column = id.getColumn();
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
    
}

void TreeMaker::DigiAcd() {
    // Purpose and Method:  Process on ACD digi event
    //   Determines the total energy deposited and
    //   store the PHA for tile 0000 in a histogram

    const TObjArray* acdDigiCol = evt->getAcdDigiCol();
    if (!acdDigiCol) return;

    Double_t totE = 0.0;
    UShort_t pha0 = 0;
    // Create an id for tile 0000
    AcdId id41(0, 0, 4, 1);

    TIter acdDigiIter(acdDigiCol);
    AcdDigi *acdDigiItem = 0;
    
    while (acdDigiItem = (AcdDigi*)acdDigiIter.Next()) {
        totE = acdDigiItem->getEnergy();
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


void TreeMaker::ReconTkr() {
  // Purpose and Method:  Process one TkrRecon event
  TkrRecon *tkrRec = rec->getTkrRecon();
  
  // If no TRKRECON data is available then return
  std::cout<<"TKR REC 2 = "<<tkrRec<<std::endl;
  
  if (!tkrRec)  return;
  
  std::cout<<"**************************************************"<<std::endl;
  
  ((TH1F*)GetObjectPtr("TKRNUMFITTRACKS"))->Fill(tkrRec->getTrackCol()->GetEntries());
  
  const TObjArray *trackCol = tkrRec->getTrackCol();
  TIter trackIter(trackCol);
  TkrKalFitTrack *track = 0;
  
    while (track = (TkrKalFitTrack*)trackIter.Next()) {
        ((TH1F*)GetObjectPtr("TKRNUMHITSPERTRACK"))->Fill(track->getNumHits());
    }
}

void TreeMaker::ReconCal() {
    // Purpose and Method:  Process on CalRecon event
    
    CalRecon *calRec = rec->getCalRecon();
    
    if (!calRec) return;
    
    
    return;
}

void TreeMaker::ReconAcd() {
    // Purpose and Method:  Processes ACD reconstruction data

    AcdRecon *acdRec = rec->getAcdRecon();
    if (!acdRec) return;

    ((TH1F*)GetObjectPtr("ACDDOCA"))->Fill(acdRec->getDoca());
    ((TH1F*)GetObjectPtr("ACDACTDIST"))->Fill(acdRec->getActiveDist());

    return;
}

void TreeMaker::HistDefine() {
    // Purpose and Method:  Setup Histograms
    
    gStyle->SetOptStat(111111);
    
    histFile = new TFile(m_histFileName,"RECREATE");
    
    //    McHistDefine();
    DigiHistDefine();
    //    ReconHistDefine();
    
}

void TreeMaker::Go(Int_t numEvents)
{    
    // Purpose and Method:  Event Loop
    //   All analysis goes here
    
    //  To read only selected branches - saves processing time
    //  Comment out any branches you are not interested in.

    // mc branches:
  if (mcTree) {
    mcTree->SetBranchStatus("*", 0);    // disable all branches
    // Activate desired branches...
        //        mcTree->SetBranchStatus("m_eventId", 1);
        //        mcTree->SetBranchStatus("m_particleCol", 1);
        //        mcTree->SetBranchStatus("m_runId", 1);        
        //        mcTree->SetBranchStatus("m_integratingHitCol", 1);        
        //        mcTree->SetBranchStatus("m_positionHitCol", 1);        
    }
    
    if (digiTree) {
        digiTree->SetBranchStatus("*",0);  // disable all branches
        // activate desired brances
        digiTree->SetBranchStatus("m_cal*",1);  
        digiTree->SetBranchStatus("m_tkr*",1);  
        digiTree->SetBranchStatus("m_acd*",1);
        digiTree->SetBranchStatus("m_eventId", 1); 
        digiTree->SetBranchStatus("m_runId", 1);
    }
    
    if (reconTree) {
        reconTree->SetBranchStatus("*",0);  // disable all branches
        // activate desired branches
	reconTree->SetBranchStatus("m_cal*", 1);  
	reconTree->SetBranchStatus("m_tkr*", 1);
	reconTree->SetBranchStatus("m_acd*", 1);
	reconTree->SetBranchStatus("m_eventId", 1); 
	reconTree->SetBranchStatus("m_runId", 1);
    }
        
    // determine how many events to process
    Int_t nentries = GetEntries();
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
    
    int target_fraction = 0;

    // BEGINNING OF EVENT LOOP
    for (Int_t ievent=m_StartEvent; ievent<nMax; ievent++, curI=ievent) {
      
      // more or less stolen from FluxSvc
      int fraction = (int)( (nMax-m_StartEvent)>0 ? 100.0*(ievent-m_StartEvent)/(nMax-m_StartEvent) : 0 );
      //        std::cout << "hurz " << ievent << " " << m_StartEvent << " " << nMax << " " << fraction << " " << target_fraction << std::endl;
      if ( fraction > target_fraction ) {
	target_fraction = fraction;
	if ( fraction >= 10 )
	  target_fraction = target_fraction/10*10 + 9;
	  std::cout << fraction << "% complete: event "<< ievent << std::endl;
        }
        
        if (mc)
	  mc->Clear();
        if (evt)
	  evt->Clear();
        if (rec)
	  rec->Clear();
        
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
	  if ( DEBUG ) std::cout << "run number: " << digiRunNum << "  event Id: " << digiEventId << std::endl;
	  
	  //////////////////////////////////////////////////
	  //	  DigiTkr();
	  const TObjArray* tkrDigiCol = 0;
	  tkrDigiCol = evt->getTkrDigiCol();
	  
	  if (tkrDigiCol)
	    {
	      
	      const int nTkrDigi = tkrDigiCol->GetEntries();
	      if ( DEBUG ) std::cout << "num tkrDigi " << nTkrDigi << std::endl;
	      
	      TIter tkrDigiIter(tkrDigiCol);
	      TkrDigi* pTkrDigi = 0;
	      
	      while ( pTkrDigi = (TkrDigi*)tkrDigiIter.Next() ) {
		
		const long int t = pTkrDigi->getTower().id();
		const long int l = pTkrDigi->getBilayer();
		const long int v = pTkrDigi->getView();
		if ( DEBUG ) std::cout << "tower, bilayer, view: " << t << " " << l << " " << v << std::endl;
		
		TString idx[3];
		idx[2] = TS + "[" + t + "][" + l + "][" + v + "]";
		idx[1] = idx[2] + "[1]";
		idx[0] = idx[2] + "[0]";
		
		const int ToT[2] = { pTkrDigi->getToT(0), pTkrDigi->getToT(1) };
		if ( DEBUG ) std::cout << "ToT: " << ToT[0] << " " << ToT[1] << std::endl;
				
		int last0strip = pTkrDigi->getLastController0Strip();
		if ( DEBUG ) std::cout << "last strip of controller 0: " << last0strip << std::endl;
		
		const int nHitCol = pTkrDigi->getNumHits();
		if ( DEBUG ) std::cout << "num hits: " << nHitCol << std::endl;
		
		if ( nHitCol > 0 ) 
		  {
		    if ( DEBUG ) std::cout << "list of hits:";
		    //            for ( std::vector<UInt_t>::iterator it = hitCol.begin(); it<hitCol.end(); ++it )
		    //                if ( DEBUG ) std::cout << " " << it;
		    
		    for ( int i=0; i<nHitCol; ++i ) 
		      {
			int hit = pTkrDigi->getHit(i);
		      }
		    
		    if ( DEBUG ) std::cout << std::endl;
		  }
		
	      }
	    }
	  //////////////////////////////////////////////////
	  // DigiTkr();
	  //            DigiCal();
	  //            DigiAcd();
        }
	
        // RECON ONLY analysis
        if (rec) 
	  { 
	    // if we have recon data proccess it
	    std::cout << "rec?"<< std::endl;
	    ReconTkr();
	    ReconCal();
	    ReconAcd();
	  } 
        
    }  // end analysis code in event loop
    
    m_StartEvent = curI;
}


void TreeMaker::CreateTree(Int_t numEvents)
{    

  
  long EventId,RunId;  
  ////////// Digi //////////    
  int TkrTotalNumHits;
  int TkrNumHits[18][2];
  int TkrHits[18][2][128];
  int ToT[2][18][2];
  int last0strip[18][2];

  ////////// clusters
  int TkrNumClus;
  int TkrClusPlane[128]; // !!
  int TkrClusView[128];
  double TkrClusX[128];
  double TkrClusY[128];
  double TkrClusZ[128];
  
  //////////  tracks
  
  int TkrNumTracks;
  int TkrTrk1Clusters[128];
  double TkrTheta[128];
  
  ////////// vertex //////////    
  int TkrNumVtx;
  double TkrVtxX, TkrVtxY, TkrVtxZ;
  //////////////////// Recon //////////////////////////////
  double Theta,Phi,ThetaXZ,ThetaYZ;
  TTree *myTree = new TTree("myTree","myTree");
  
  // mc branches:
  if (mcTree) {
    mcTree->SetBranchStatus("*", 0);    // disable all branches
    // Activate desired branches...
    //        mcTree->SetBranchStatus("m_eventId", 1);
    //        mcTree->SetBranchStatus("m_particleCol", 1);
    //        mcTree->SetBranchStatus("m_runId", 1);        
    //        mcTree->SetBranchStatus("m_integratingHitCol", 1);        
    //        mcTree->SetBranchStatus("m_positionHitCol", 1);        
  }
  
  if (digiTree) {
      digiTree->SetBranchStatus("*",0);  // disable all branches
      // activate desired brances
      digiTree->SetBranchStatus("m_cal*",1);  
      digiTree->SetBranchStatus("m_tkr*",1);  
      digiTree->SetBranchStatus("m_acd*",1);
      digiTree->SetBranchStatus("m_eventId", 1); 
      digiTree->SetBranchStatus("m_runId", 1);
      
      myTree->Branch("EventId",&EventId,"EventId/I");
      myTree->Branch("RunId",&RunId,"RunId/I");
      
      myTree->Branch("ToT",ToT,"ToT[2][18][2]/I");
      
      myTree->Branch("TkrTotalNumHits",&TkrTotalNumHits,"TkrTotalNumHits/I");
      myTree->Branch("TkrNumHits",TkrNumHits,"TkrNumHits[18][2]/I");
      myTree->Branch("TkrHits",TkrHits,"TkrHits[18][2][128]/I");
      myTree->Branch("last0strip",last0strip,"last0striplast0strip[18][2]/I");
  }
  
  if (reconTree) {
      reconTree->SetBranchStatus("*",0);  // disable all branches
      // activate desired branches
      reconTree->SetBranchStatus("m_cal*", 1);  
      reconTree->SetBranchStatus("m_tkr*", 1);
      reconTree->SetBranchStatus("m_acd*", 1);
      reconTree->SetBranchStatus("m_eventId", 1); 
      reconTree->SetBranchStatus("m_runId", 1);
      
      ////////// Clusters: ////////////////////////////////////////
      myTree->Branch("TkrNumClus",&TkrNumClus,"TkrNumClus/I");
      myTree->Branch("TkrClusX",TkrClusX,"TkrClusX[TkrNumClus]/D");
      myTree->Branch("TkrClusY",TkrClusY,"TkrClusY[TkrNumClus]/D");
      myTree->Branch("TkrClusZ",TkrClusZ,"TkrClusZ[TkrNumClus]/D");
      myTree->Branch("TkrClusPlane",TkrClusPlane,"TkrClusPlane[TkrNumClus]/I");
      myTree->Branch("TkrClusView",TkrClusView,"TkrClusView[TkrNumClus]/I");
      /////////// Tracks: ////////////////////////////////////////
      myTree->Branch("TkrNumTracks",&TkrNumTracks,"TkrNumTracks/I");
      myTree->Branch("TkrTrk1Clusters",TkrTrk1Clusters,"TkrTrk1Clusters[TkrNumClus]/I");
      myTree->Branch("TkrTheta",TkrTheta,"TkrTheta[TkrNumTracks]/D");
      
      /////////// Vertex: ////////////////////////////////////////
      myTree->Branch("TkrNumVtx",&TkrNumVtx,"TkrNumVtx/I");
      myTree->Branch("TkrVtxX",&TkrVtxX,"TkrVtxX/D");
      myTree->Branch("TkrVtxY",&TkrVtxY,"TkrVtxY/D");
      myTree->Branch("TkrVtxZ",&TkrVtxZ,"TkrVtxZ/D");
      
      //////////Recon //////////
      myTree->Branch("Theta",&Theta,"Theta/D");
      myTree->Branch("Phi",&Phi,"Phi/D");
      myTree->Branch("ThetaXZ",&ThetaXZ,"ThetaXZ/D");
      myTree->Branch("ThetaYZ",&ThetaYZ,"ThetaYZ/D");
  }
  
  // determine how many events to process
  Int_t nentries = GetEntries();
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
  
  int target_fraction = 0;
  //////////////////////////////////////////////////
  // BEGINNING OF EVENT LOOP
  for (Int_t ievent=m_StartEvent; ievent<nMax; ievent++, curI=ievent) {
      for(int l=0;l<18; l++) {
	  for(int v=0;v<2; v++) {
	      TkrTotalNumHits=0;    
	      TkrNumHits[l][v]=0;
	      
	      for(int h=0;h<128; h++)
                  TkrHits[l][v][h]=-1;
	      ToT[0][l][v]=-1;
	      ToT[1][l][v]=-1;
	      last0strip[l][v]=-1;
          }
      }
      for ( int h=0; h<128; ++h )
          TkrTrk1Clusters[h] = -1;

      // more or less stolen from FluxSvc
      int fraction = (int)( (nMax-m_StartEvent)>0 ? 100.0*(ievent-m_StartEvent)/(nMax-m_StartEvent) : 0 );
      if ( fraction > target_fraction ) {
	target_fraction = fraction;
	if ( fraction >= 10 )
	  target_fraction = target_fraction/10*10 + 9;
	std::cout << fraction << "% complete: event "<< ievent << std::endl;
      }
      
      if (mc)
	mc->Clear();
      if (evt)
	evt->Clear();
      if (rec)
	rec->Clear();
      
      digiEventId = 0; reconEventId = 0; mcEventId = 0;
      digiRunNum = 0; reconRunNum = 0; mcRunNum = 0;
      
      nb = GetEvent(ievent);
      
      nbytes += nb;
      //////////////////////////////////////////////////
      // Monte Carlo ONLY analysis
      if (mc) {  // if we have mc data process it
	mcEventId = mc->getEventId();
	mcRunNum  = mc->getRunId();
	McData();
      } 
      
      // Digi ONLY analysis
      if (evt) 
	{
	  digiEventId = evt->getEventId(); 
	  digiRunNum = evt->getRunId();
	  if ( DEBUG ) std::cout << "run number: " << digiRunNum << "  event Id: " << digiEventId <<" ievent ="<< ievent <<std::endl;
	  EventId = ievent+1;
	  RunId   = digiRunNum;
	  //////////////////////////////////////////////////
	  // digitkr:
	  
	  const TObjArray* tkrDigiCol = 0;
	  tkrDigiCol = evt->getTkrDigiCol();
	  //////////////////////////////////////////////////
	  
	  if (tkrDigiCol)
	    {
	      
	      const int nTkrDigi = tkrDigiCol->GetEntries();
	      if ( DEBUG ) std::cout << "num tkrDigi " << nTkrDigi << std::endl;
	      
	      TIter tkrDigiIter(tkrDigiCol);
	      TkrDigi* pTkrDigi = 0;
	      

	      while ( pTkrDigi = (TkrDigi*)tkrDigiIter.Next() ) {
		
		int Tower = pTkrDigi->getTower().id();
		int Layer = pTkrDigi->getBilayer();
		int View  = pTkrDigi->getView();

		if ( DEBUG ) std::cout << "tower, bilayer, view: " << Tower << " " << Layer << " " << View << std::endl;
		const int TOT[2] = { pTkrDigi->getToT(0), pTkrDigi->getToT(1) };
		
		ToT[0][Layer][View]  = TOT[0];
		ToT[1][Layer][View]  = TOT[1];
		
		if ( DEBUG ) std::cout << "ToT: " << ToT[0][Layer][View]  << " " << ToT[1][Layer][View]  << std::endl;
		
		//last0strip = pTkrDigi->getLastController0Strip();
		//		if ( DEBUG ) std::cout << "last strip of controller 0: " << last0strip << std::endl;
		TkrNumHits[Layer][View]  = pTkrDigi->getNumHits();
		TkrTotalNumHits += TkrNumHits[Layer][View];

		if ( DEBUG ) std::cout << "num hits: " << TkrNumHits[Layer][View]  << std::endl;
				
		if ( TkrNumHits[Layer][View] > 0 ) 
		  {
		    for ( int i=0; i<TkrNumHits[Layer][View] ; ++i ) 
		      TkrHits[Layer][View][i] = pTkrDigi->getHit(i);
		    
		    
		    
		    if ( DEBUG ) 
		      { 
			std::cout << "list of hits:";
			for ( int i=0; i<TkrNumHits[Layer][View] ; ++i ) std::cout << " " << TkrHits[Layer][View][i];
			std::cout << std::endl;
		      }
		    
		    
		  }
		
	      }
	    }
	  
	  
	  //////////////////////////////////////////////////
	  //////////////////////////////////////////////////
	  //            DigiCal();
	  //            DigiAcd();
	  
	}
      // RECON ONLY analysis
      if (rec) 
	{  // if we have recon data proccess it
	  TkrRecon *tkrRec = rec->getTkrRecon();
	  if (tkrRec)
	    {
	      ////////////////// TKR CLUSTER: ////////////////////////////////
	      const TObjArray *clusCol = tkrRec->getClusterCol();
	      TkrNumClus = clusCol->GetEntries();

	      TIter tkrClusIter(clusCol);
	      TkrCluster* pTkrClus = 0;
	      int clusIdx=0;
	      while ( pTkrClus = (TkrCluster*) tkrClusIter.Next() ) 
		{
		  TkrClusX[clusIdx]     = pTkrClus->getPosition().X();
		  TkrClusY[clusIdx]     = pTkrClus->getPosition().Y();
		  TkrClusZ[clusIdx]     = pTkrClus->getPosition().Z();
		  TkrClusPlane[clusIdx] = pTkrClus->getPlane();
		  TkrClusView[clusIdx]  = pTkrClus->getView();
		  clusIdx++;
		}
	      ////////////////// TKR TRACKS: ////////////////////////////////
	      const TObjArray *tkrCol = tkrRec->getTrackCol();
	      TkrNumTracks = tkrCol->GetEntries();
 
              if ( TkrNumTracks > 0 ) {
                  const TkrKalFitTrack* track1 = (TkrKalFitTrack*)tkrCol->First();
                  const int nHits = track1->getNumHits();
                  for ( int i=0; i<nHits; ++i )
                      TkrTrk1Clusters[i] = track1->getHitPlane(i)->getIdHit();
              }
 	      
	      ////////////////// TLKR VERTEX: ////////////////////////////////

	      const TObjArray *vtxCol = tkrRec->getVertexCol();

	      TkrNumVtx = vtxCol->GetEntries();

	      TkrVertex *vtx = (TkrVertex*) vtxCol->First();
	      
	      if(vtx)
		{
		  TVector3 dir = vtx->getDirection();
		  
		  TkrVtxX = dir.X();
		  TkrVtxY = dir.Y();
		  TkrVtxZ = dir.Z();
		  Theta = dir.Theta();
		  Phi =dir.Phi();
		  if(dir.Z()==0)
		    {

		  ThetaXZ = 0.0;
		  ThetaYZ = 0.0;
		    }
		  else
		    {
		      ThetaXZ = dir.X()/dir.Z();
		      ThetaYZ = dir.Y()/dir.Z();
		    }
		  
		    if ( DEBUG ) std::cout<<TkrVtxX<<" "<<TkrVtxY<<" "<<TkrVtxZ<<std::endl;
		}
	    }
	  
	  //	  ReconTkr();
	  //	  ReconCal();
	  //      ReconAcd();
	  
	} 
      myTree->Fill();
    }  // end analysis code in event loop
  
  m_StartEvent = curI;
  
  TFile myFile("MyRootFile.root","RECREATE");
  myTree->Write();
  myFile.Close();
  
}

#endif
