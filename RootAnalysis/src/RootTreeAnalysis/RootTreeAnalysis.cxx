#ifndef RootTreeAnalysis_cxx
#define RootTreeAnalysis_cxx 1

#include "RootTreeAnalysis.h"

TObjArray tkrLayerHistArr;
UInt_t digiEventId, reconEventId, mcEventId;
UInt_t digiRunNum, reconRunNum, mcRunNum;

void RootTreeAnalysis::McHistDefine() {
    // Purpose and Method:  Monte Carlo histogram definitions
    
    // Must cd into the histFile to be sure the histograms are created within
    // this histogram file.
    histFile->cd();
    
    TH1F *EVENTMC = new TH1F("EVENTMC", "MC Event Id",
        2000, 0, 2000);
    
    TH1F *RUNMC = new TH1F("RUNMC", "MC Run Id",
        2000, 0, 2000);
    
    TH1F *PARTCOUNTMC = new TH1F("PARTCOUNTMC", "MC Part Count",
        20, 0, 600);
    TH1F *POSCOUNTMC = new TH1F("POSCOUNTMC", "MC Pos Count",
        20, 0, 100);
    
    TH1F *POSENERGYDEP = new TH1F("POSENERGYDEP", "MC Pos Edep",
        25, 0, 0.5);
    
    TH1F *POSMCTYPE = new TH1F("POSMCTYPE", "MC Pos MC pTypes",
        40, -20, 20);
    
    TH1F *INTCOUNTMC = new TH1F("INTCOUNTMC", "MC Int Count",
        20, 0, 40);
    
    TH1F *INTENERGYDEPCAL = new TH1F("INTENERGYDEPCAL", "MC Int Xtal Etot",
        20, 0, 20);
    TH1F *INTCALMCTYPE = new TH1F("INTCALMCTYPE", "MC Int MC pTypes",
        100, -100, 100);
    
    TH1F *INTENERGYDEPACD = new TH1F("INTENERGYDEPACD", "MC Int Scint Etot",
        20, 0, 20);
    TH1F *INTENERGYTOT = new TH1F("INTENERGYTOT", "MC Int Etot",
        20, 0, 500);
    
    TH1F *ENERGYMC = new TH1F("ENERGYMC", "MC Part Energy (MeV)",
        500, 0, 50000);   
    
}

void RootTreeAnalysis::DigiHistDefine() {
    // Purpose and Method:  Digitization histogram definitions
    // Must cd into the histFile to be sure these new histograms are created in
    // this histogram ROOT file
    histFile->cd();
    
    TH1F *CALDIGICOUNT = new TH1F("CALDIGICOUNT", "Cal Digi multiplicity",
        50, 0, 50);
    
    TH1F *CALADC = new TH1F("CALADC", "Cal Digi ADC - both faces",
        200, 0, 1000);
    TH1F *CALRANGE = new TH1F("CALRANGE", "Cal Digi Range - both faces",
        10, 0, 10);
    TH1F *CALEAVE = new TH1F("CALEAVE", "Cal Digi Energy - faces/2",
        200, 0, 1000);
    
    TH1F *CALEAVETOTAL = new TH1F("CALEAVETOTAL", "Cal Digi Energy - faces/2 summed",
        200, 0, 10000);
    
    TH1F *CALLAYER = new TH1F("CALLAYER", "Cal Digi Layer",
        15, 0, 15);
    TH1F *CALTOWER = new TH1F("CALTOWER", "Cal Digi Tower",
        20, 0, 20);
    TH1F *CALCOLUMN = new TH1F("CALCOLUMN", "Cal Digi Column",
        20, 0, 20);
    
    TH1F *CALELYR0 = new TH1F("CALELYR0", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR1 = new TH1F("CALELYR1", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR2 = new TH1F("CALELYR2", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR3 = new TH1F("CALELYR3", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR4 = new TH1F("CALELYR4", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR5 = new TH1F("CALELYR5", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR6 = new TH1F("CALELYR6", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    TH1F *CALELYR7 = new TH1F("CALELYR7", "Cal Digi Energy - faces/2 - Layer 0",
        200, 0, 1000);
    
    
    TH1F *CALNLYR0 = new TH1F("CALNLYR0", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR1 = new TH1F("CALNLYR1", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR2 = new TH1F("CALNLYR2", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR3 = new TH1F("CALNLYR3", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR4 = new TH1F("CALNLYR4", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR5 = new TH1F("CALNLYR5", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR6 = new TH1F("CALNLYR6", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);
    TH1F *CALNLYR7 = new TH1F("CALNLYR7", "Cal Digi Energy - faces/2 - Layer 0",
        10, 0, 10);

    TH1F *ACDPHATILE41 = new TH1F ("ACDPHATILE41", "PHA for Tile 041", 
        50, 0, 500);

    TH1F *ACDTOTE = new TH1F("ACDTOTE", "ACD Total Energy",
        20, 0, 10);
    
}

void RootTreeAnalysis::ReconHistDefine() {
    // Purpose and Method:  Reconstruction histogram defintions
    // Must cd into the histFile to be sure these new histograms are created
    // within this histogram ROOT file.
    histFile->cd();

    TH1F *TKRNUMFITTRACKS = new TH1F("TKRNUMFITTRACKS", "Number of Fit Tracks",
        30, 0, 30);

    TH1F *TKRNUMHITSPERTRACK = new TH1F("TKRNUMHITSPERTRACK", "Number of Hits/Tracks",
        40, 0, 40);

    TH1F *ACDDOCA = new TH1F("ACDDOCA", "ACD DOCA",
        40, 0, 200);

    TH1F *ACDACTDIST = new TH1F("ACDACTDIST", "ACD Active Distance",
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
    
    const TObjArray* tkrDigiCol = evt->getTkrDigiCol();
    if (!tkrDigiCol) return;

    return;
}

void RootTreeAnalysis::DigiCal() {
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
    
    while (acdDigiItem = (AcdDigi*)acdDigiIter.Next()) {
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


void RootTreeAnalysis::ReconTkr() {
    // Purpose and Method:  Process one TkrRecon event
    
    TkrRecon *tkrRec = rec->getTkrRecon();
    
    // If no TRKRECON data is available then return
    if (!tkrRec)  return;

    ((TH1F*)GetObjectPtr("TKRNUMFITTRACKS"))->Fill(tkrRec->getTrackCol()->GetEntries());

    const TObjArray *trackCol = tkrRec->getTrackCol();
    TIter trackIter(trackCol);
    TkrKalFitTrack *track = 0;
    
    while (track = (TkrKalFitTrack*)trackIter.Next()) {
        ((TH1F*)GetObjectPtr("TKRNUMHITSPERTRACK"))->Fill(track->getNumHits());
    }
}

void RootTreeAnalysis::ReconCal() {
    // Purpose and Method:  Process on CalRecon event
    
    CalRecon *calRec = rec->getCalRecon();
    
    if (!calRec) return;
    
    
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
        mcTree->SetBranchStatus("*", 0);    // disable all branches
        // Activate desired branches...
        mcTree->SetBranchStatus("m_eventId", 1);
        mcTree->SetBranchStatus("m_particleCol", 1);
        mcTree->SetBranchStatus("m_runId", 1);        
        mcTree->SetBranchStatus("m_integratingHitCol", 1);        
        mcTree->SetBranchStatus("m_positionHitCol", 1);        
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
