#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCollection.h"
#include "digiRootData/DigiEvent.h"

/** @file testDigiClasses.cxx
* @brief This defines a test routine for the Digitization ROOT classes.
*
* This program create a new Digi ROOT file, and the opens it up again
* for reading.  The contents are checked against the values known to be stored
* during the original writing of the file.
* The contents of the file are printed to the screen.
* The program returns 0 if the test passed.
* If failure, the program returns -1.
*
* $Header$
*/
const UInt_t runNum = 1;
const Int_t numXtals = 10;
const UInt_t numDigi = 12;
const UInt_t numAcd = 11;
const Int_t numCalDiag = 4;
const Int_t numTkrDiag = 2;
const Bool_t fromMc = true;
const Int_t numTem = 5;
Float_t randNum;

bool floatInRange(Double_t actual, Double_t desired) {
  const Double_t fudge=0.00001;
  if ( (actual >= (desired - fudge)) && (actual <= (desired+fudge)) ){
    return true;
  }
  return false;

}

int checkDigiEvent(DigiEvent *evt, UInt_t ievent) {
    if (evt->getRunId() != runNum) {
        std::cout << "Run Id is wrong: " << evt->getRunId() << std::endl;
        return -1;
    }
    if (evt->getEventId() != ievent) {
        std::cout << "Event Id is wrong: " << evt->getEventId() << std::endl;
        return -1;
    }

    if (!floatInRange(evt->getTimeStamp(), randNum*ievent) ) {
        std::cout << "Time Stamp is wrong: " << evt->getTimeStamp() << std::endl;
        return -1;
    }

    if (!floatInRange(evt->getLiveTime(), randNum*ievent) ) {
        std::cout << "LiveTime is wrong: " << evt->getLiveTime() << std::endl;
        return -1;
    }
    if (evt->getFromMc() != fromMc) {
        std::cout << "From MC flag is wrong" << std::endl;
        return -1;
    }
    return 0;
}

int checkTem(Tem* curTem) {
    const ErrorData& err = curTem->getError();
    if (err.getCal() != 1) {
        std::cout << "Error Data Cal is wrong" << std::endl;
        return -1;
    }

    if (err.getTkr() != 2) {
       std::cout << "Error Data Tkr is wrong " << std::endl;
       return -1;
    }

    if (err.getPhs() != true) {
        std::cout << "Phased Error is wrong " << std::endl;
        return -1;
    }

    if (err.getTmo() != false) {
        std::cout << "Timeout is wrong" << std::endl;
        return -1;
    }

    return 0;
}


int checkCalDiagnostic(const CalDiagnosticData *calDiag) {
    if (calDiag->getDataWord() != 10101) {
        std::cout << "CalDiagnostic data word is wrong" << std::endl;
        return -1;
    }
    if (calDiag->tower() != 2) {
        std::cout << "CalDiagnostic tower is wrong" << std::endl;
        return -1;
    }
    if (calDiag->layer() != 1) {
        std::cout << "CalDiagnostic layer is wrong" << std::endl;
        return -1;
    }

    return 0;
}

int checkTkrDiagnostic(const TkrDiagnosticData *tkrDiag) {
    if (tkrDiag->getDataWord() != 20301) {
        std::cout << "TkrDiagnostic datum is wrong" << std::endl;
        return -1;
    }
    if (tkrDiag->tower() != 4) {
        std::cout << "TkrDiagnostic tower is wrong" << std::endl;
        return -1;
    }
    if (tkrDiag->gtcc() != 2) {
        std::cout << "TkrDiagnostic gtcc is wrong" << std::endl;
        return -1;
    }
   

    return 0;
}

int checkL1T(const L1T &level1) {

    if (level1.getTriggerWord() != 19) {
        std::cout << "Trigger Word is wrong: " << level1.getTriggerWord() 
            << std::endl;
        return -1;
    }

    if (level1.getAcdLow() != true) {
        std::cout << "Acd Low is False when it should be True" << std::endl;
        return -1;
    }

    if (level1.getAcdHigh() != true) {
        std::cout << "Acd High is false when it should be true" << std::endl;
        return -1;
    }

    if (level1.getTkr3InARow() != true) {
        std::cout << "Tkr 3-in-a-row is False when it should be True" << std::endl;
        return -1;
    }

    if (level1.getCalLow() != false) {
        std::cout << "Cal Low is true when it should be false" << std::endl;
        return -1;
    }

    if (level1.getCalHigh() != false) {
        std::cout << "Cal High is True when it shoule be False" << std::endl;
        return -1;
    }


    unsigned int iTower;
    for (iTower = 0; iTower < 16; iTower++) {
        if (level1.getDigiTriRowBits(iTower) != (iTower%2)) {
            std::cout << "DigiTriRowBits is wrong for tower: " << iTower << std::endl;
            return -1;
        }
        if (level1.getTrgReqTriRowBits(iTower) != (iTower%2)) {
            std::cout << "TrgReqTriRowBits is wrong for tower: " << iTower << std::endl;
            return -1;
        }
   
    }

    return 0;
}

int checkCalDigi(CalDigi *digi, UInt_t ievent) {
    // Checks the contents of one CalDigi object

    CalXtalId id = digi->getPackedId();
    if ( (id.getTower() != 5) || (id.getLayer() != 4) || (id.getColumn() != 3) ) {
        std::cout << "CalDigi id is wrong: " << id.getTower() << ", "
        << id.getLayer() << ", " << id.getColumn() << std::endl;
        return -1;
    }
    CalXtalId::CalTrigMode mode = digi->getMode();
    if (mode != CalXtalId::BESTRANGE) {
        std::cout << "Cal Digi Mode is wrong: " << mode << std::endl;
        return -1;
    }

    
    const CalXtalReadout *readoutCol = digi->getReadoutCol();
    UInt_t numEntries = digi->getNumReadouts();
    if (numEntries != 1) {
        std::cout << "Wrong number of readouts in Cal Digi " 
            << numEntries << std::endl;
        return -1;
    }
    const CalXtalReadout *calReadout;

    unsigned int iReadout;
    for (iReadout = 0; iReadout<numEntries; iReadout++){
        calReadout = &(readoutCol[iReadout]);
        calReadout->Print();
        UShort_t adcP = calReadout->getAdc(CalXtalId::POS);
        UShort_t adcM = calReadout->getAdc(CalXtalId::NEG);
        Char_t rangeM = calReadout->getRange(CalXtalId::NEG);
        Char_t rangeP = calReadout->getRange(CalXtalId::POS);
        if ((adcP != 4095) || (adcM != 4095) ) {
            std::cout << "adcP or adcM are wrong: " << adcP << " "
                << adcM << std::endl;
            return -1;
        }
        if ( (rangeM != CalXtalId::LEX8) || (rangeP != CalXtalId::HEX8)) {
            std::cout << "rangeP or rangeM are wrong: " << rangeP << " "
                << rangeM << std::endl;
            return -1;
        }
    }

    // Checking results of calling the getAdcSelectedRange method
    Short_t adcP_lex8 = digi->getAdcSelectedRange(CalXtalId::LEX8, CalXtalId::POS);
    Short_t adcP_hex8 = digi->getAdcSelectedRange(CalXtalId::HEX8, CalXtalId::POS);
    Short_t adcP_lex1 = digi->getAdcSelectedRange(CalXtalId::LEX1, CalXtalId::POS);
    Short_t adcP_hex1 = digi->getAdcSelectedRange(CalXtalId::HEX1, CalXtalId::POS);
    Short_t adcM_lex8 = digi->getAdcSelectedRange(CalXtalId::LEX8, CalXtalId::NEG);
    Short_t adcM_hex8 = digi->getAdcSelectedRange(CalXtalId::HEX8, CalXtalId::NEG);
    Short_t adcM_lex1 = digi->getAdcSelectedRange(CalXtalId::LEX1, CalXtalId::NEG);
    Short_t adcM_hex1 = digi->getAdcSelectedRange(CalXtalId::HEX1, CalXtalId::NEG);

    if ( (adcP_lex8 != (Char_t)-1) || (adcP_lex1 != (Char_t)-1) ||
        (adcP_hex1 != (Char_t)-1) ){
        std::cout << "adcP defined for an incorrect range" << std::endl;
        return -1;
    }

    if ( (adcM_hex8 != (Char_t)-1) || (adcM_lex1 != (Char_t)-1) ||
        (adcM_hex1 != (Char_t)-1) ){
        std::cout << "adcM defined for an incorrect range" << std::endl;
        return -1;
    }


    if ( ((UShort_t)adcP_hex8 != 4095) || ( (UShort_t)adcM_lex8 != 4095) ) {
        std::cout << "adcP_hex8 or adcM_lex8 is not defined correctly: "
            << (UShort_t)adcP_hex8 << " " << (UShort_t)adcM_lex8 << std::endl;
        return -1;
    }
    std::cout << "All adc values over all ranges: " << std::endl;
    std::cout << "Pos: lex8, hex8, lex1, hex1 " << adcP_lex8 << " "
        << (UShort_t)adcP_hex8 << " " << adcP_lex1 << " " 
        << adcP_hex1 << std::endl;
    std::cout << "Neg: lex8, hex8, lex1, hex1 " << (UShort_t)adcM_lex8 << " "
        << adcM_hex8 << " " << adcM_lex1 << " " << adcM_hex1 << std::endl;

    return 0;
}

int checkTkrDigi(TkrDigi *digi, UInt_t ievent, Int_t idigi) {

    TowerId id = digi->getTower();
    if ((id.ix() != 3) || (id.iy() != 2) ) {
        std::cout << "TkrDigi id is wrong (x,y): (" << id.ix() << "," 
            << id.iy() << ")" << std::endl;
        return -1;
    }

    if (id.id() != 3+4*2)  {
        std::cout << "TkrDigi id is computer wrong: " << id.id() << std::endl;
        return -1;
    }

    UShort_t bilayer = digi->getBilayer();
    if (bilayer != idigi) {
        std::cout << "TkrDigi bilayer is wrong:  " << bilayer << std::endl;
        return -1;
    }

    GlastAxis::axis a = digi->getView();
    if (a != GlastAxis::Y) {
        std::cout << "TkrDigi Axis is wrong: " << a << std::endl;
        return -1;
    }

    if (digi->getNumHits() != idigi+1) {
        std::cout << "TkrDigi number of hit strip is wrong: " << digi->getNumHits() << std::endl;
        return -1;
    }
    Int_t istrip;
    for (istrip = 0; istrip < idigi; istrip++) {
        Int_t strip = digi->getHit(istrip);
        if (strip != istrip) {
            std::cout << "TkrDigi strip is wrong: " << strip << std::endl;
            return -1;
        }
    }
    Int_t strip = digi->getHit(idigi);
    if (strip != idigi*2) {
        std::cout << "TkrDigi strip is wrong: " << strip << std::endl;
        return -1;
    }

    strip = digi->getHit(idigi+1);
    if (strip != -1) {
        std::cout << "TkrDigi non-existant strip is defined: " << digi->getHit(idigi+1) << std::endl;
        return -1;
    }
    
    if (digi->getToT(0) != idigi) {
        std::cout << "TkrDigi ToT 0 is wrong: " << digi->getToT(0) << std::endl;
        return -1;
    }
    
    if (digi->getToT(1) != idigi+1) {
        std::cout << "TkrDigi ToT 1 is wrong: " << digi->getToT(1) << std::endl;
        return -1;
    }
    
    if (digi->getToT(2) != -1) {
        std::cout << "TkrDigi non-existant ToT is defined: " << digi->getToT(2) << std::endl;
        return -1;
    }
    
    return 0;
}

int checkGem(const Gem &gem) {
   GemTileList list = gem.getTileList();
   if (list.getXzm() != 20) return -1;
   if (list.getXzp() != 21) return -1;
   if (list.getYzm() != 22) return -1;
   if (list.getYzp() != 23) return -1;
   if (list.getXy() != 24) return -1; 
   if (list.getRbn() != 25) return -1;
   if (list.getNa() != 26) return -1;
   GemOnePpsTime time = gem.getOnePpsTime();
   if (time.getTimebase() != 50) return -1;
   if (time.getSeconds() != 100) return -1;

   if (gem.getTkrVector() != 1) return -1;
   if (gem.getRoiVector() != 2) return -1;
   if (gem.getCalLeVector() != 3) return -1;
   if(gem.getCalHeVector() != 4) return -1;
   if (gem.getCnoVector() != 5) return -1;
   if (gem.getConditionSummary() != 6) return -1;
   if (gem.getMissed() != 1) return -1;
   if (gem.getLiveTime() != 7) return -1;
   if (gem.getPrescaled() != 8) return -1;
   if (gem.getDiscarded() != 9) return -1;
   if (gem.getCondArrTime().condArr() != 10) return -1;
   if (gem.getTriggerTime() != 11) return -1;
   if (gem.getDeltaEventTime() != 12) return -1;
   if (gem.getDeltaWindowOpenTime() != 13) return -1;

   return 0; 
}

int checkAcdDigi(AcdDigi *digi, UInt_t ievent, UInt_t idigi) {

    Float_t f = ievent;

    if (!floatInRange(digi->getEnergy(), ievent*randNum) ) {
        std::cout << "AcdDigi Energy is wrong " << digi->getEnergy() << std::endl;
        return -1;
    }

    AcdId id = digi->getId();
    if (id.getLayer() != 0) {
        std::cout << "AcdId Layer is wrong " << id.getLayer() << std::endl;
        return -1;
    }

    VolumeIdentifier volId = digi->getVolId();
    if ( (volId.getBits0to31() != 0) || (volId.getBits32to63() != 16777216) || (volId.size() != 1) ) {
        std::cout << "AcdDigi VolId is incorrect" << std::endl;
        return -1;
    }

    if (id.getFace() != 2) {
        std::cout << "AcdId Face is wrong " << id.getFace() << std::endl;
        return -1;
    }

    if (id.getRow() != 5) {
        std::cout << "AcdId Row is wrong " << id.getRow() << std::endl;
        return -1;
    }

    if (id.getColumn() != 4) {
        std::cout << "AcdId Col is wrong " << id.getColumn() << std::endl;
        return -1;
    }

    UShort_t pha_A = digi->getPulseHeight(AcdDigi::A);
    if (pha_A != 4095) {
        std::cout << "AcdDigi PHA A is wrong " << pha_A << std::endl;
        return -1;
    }

    Bool_t veto_A = digi->getVeto(AcdDigi::A);
    if (veto_A != kFALSE) {
        std::cout << "Veto A is wrong " << veto_A << std::endl;
        return -1;
    }

    Bool_t low_A = digi->getLowDiscrim(AcdDigi::A);
    if (low_A != kTRUE) {
        std::cout << "Low A is wrong " << low_A << std::endl;
        return -1;
    }

    Bool_t high_A = digi->getHighDiscrim(AcdDigi::A);
    if (high_A != kFALSE) {
        std::cout << "High A is wrong " << high_A << std::endl;
        return -1;
    }

    UShort_t pha_B = digi->getPulseHeight(AcdDigi::B);
    if (pha_B != 1) {    
        std::cout << "AcdDigi PHA B is wrong " << pha_B << std::endl;
        return -1;
    }

    Bool_t veto_B = digi->getVeto(AcdDigi::B);
    if (veto_B != kTRUE) {
        std::cout << "Veto B is wrong " << veto_B << std::endl;
        return -1;
    }

    Bool_t low_B = digi->getLowDiscrim(AcdDigi::B);
    if (low_B != kTRUE) {
        std::cout << "Low B is wrong " << low_B << std::endl;
        return -1;
    }

    Bool_t high_B = digi->getHighDiscrim(AcdDigi::B);
    if (high_B != kTRUE) {
        std::cout << "High B is wrong " << high_B << std::endl;
        return -1;
    }

    return 0;
}


/// Read in the ROOT file just generated via the write method
int read(char* fileName, int numEvents) {
    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("Digi");
    DigiEvent *evt = 0;
    t->SetBranchAddress("DigiEvent", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    Int_t ievent;
    for (ievent = 0; ievent < numEvents; ievent++) {
        t->GetEvent(ievent);
        std::cout << "DigiEvent ievent = " << ievent << std::endl;
        evt->Print();
        if (checkDigiEvent(evt, ievent) < 0) return -1;
        
        L1T level1Trigger = evt->getL1T();
        if (checkL1T(level1Trigger) < 0) return -1;

        Gem gemData = evt->getGem();
        if (checkGem(gemData) < 0) return -1;

        UInt_t numCalDigi = evt->getCalDigiCol()->GetEntries();
        const TObjArray *calDigiCol = evt->getCalDigiCol();
        TIter calDigiIt(calDigiCol);
        CalDigi *digi = 0;
        while((digi = (CalDigi*)calDigiIt.Next())!=0) {
            digi->Print();
            if (checkCalDigi(digi, ievent) < 0) return -1;
        }
        const TObjArray *tkrDigiCol = evt->getTkrDigiCol();
        TIter tkrDigiIt(tkrDigiCol);
        TkrDigi *tDigi = 0;
        UInt_t idigi = 0;
        while ((tDigi = (TkrDigi*)tkrDigiIt.Next())!=0) {
            tDigi->Print();
            if (checkTkrDigi(tDigi, ievent, idigi) < 0) return -1;
            idigi++;
        }

        const TClonesArray *acdDigiCol = evt->getAcdDigiCol();
        TIter acdDigiIt(acdDigiCol);
        AcdDigi *aDigi = 0;
        idigi = 0;
        while ((aDigi = (AcdDigi*)acdDigiIt.Next())!=0) {
            aDigi->Print();
            if (checkAcdDigi(aDigi, ievent, idigi) < 0) return -1;
            idigi++;
        }


        const TClonesArray *calDiagCol = evt->getCalDiagnosticCol();
        if (calDiagCol->GetEntries() != numCalDiag) return -1;
        TIter calDiagIt(calDiagCol);
        CalDiagnosticData *cDiag = 0;
        while ((cDiag=(CalDiagnosticData*)calDiagIt.Next())!=0) {
            cDiag->Print();
            if (checkCalDiagnostic(cDiag) < 0) return -1;
        }

        const TClonesArray *tkrDiagCol = evt->getTkrDiagnosticCol();
        if (tkrDiagCol->GetEntries() != numTkrDiag) return -1;
        TIter tkrDiagIt(tkrDiagCol);
        TkrDiagnosticData *tDiag = 0;
        while ((tDiag=(TkrDiagnosticData*)tkrDiagIt.Next())!=0) {
            tDiag->Print();
            if (checkTkrDiagnostic(tDiag) < 0) return -1;
        }

        const TClonesArray *temCol = evt->getTemCol();
        TIter temIt(temCol);
        Tem* curTem = 0;
        while ((curTem = (Tem*)temIt.Next())) {
            if (checkTem(curTem) < 0) return -1;
        }
 }
    
    f->Close();
    delete f;
    
    return 0;
}

/// Create a new Monte Carlo ROOT file
int write(char* fileName, int numEvents) {
    Int_t buffer = 64000;
    Int_t splitLevel = 1;
    
    TFile *f =  new TFile(fileName, "RECREATE");
    TTree *t = new TTree("Digi", "Digi");
    DigiEvent *ev = new DigiEvent();
    t->Branch("DigiEvent", "DigiEvent", &ev, buffer, splitLevel);
    
    std::cout << "Created new ROOT file" << std::endl;
    
    TRandom randGen;
    Int_t ievent, ixtal;
    randNum = randGen.Rndm();

    for (ievent = 0; ievent < numEvents; ievent++) {
        
        UInt_t digiTriRowBits[16] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
        UInt_t trgReqTriRowBits[16] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
        L1T level1(19, digiTriRowBits,trgReqTriRowBits);
        EventSummaryData summary(0);
        summary.initEventFlags(1);
        ev->initialize(ievent, runNum, randNum*ievent, randNum*ievent, level1, summary, fromMc);
        Gem myGem;
        GemTileList tileList(20, 21, 22, 23, 24, 25, 26);
        myGem.initTrigger(1, 2, 3, 4, 5, 6, 1, tileList);
        GemOnePpsTime ppsTime(50, 100);
        myGem.initSummary(7, 8, 9, 10, 11, ppsTime, 12, 13);
        ev->initGem(myGem);
        
        for (ixtal = 0; ixtal < numXtals; ixtal ++) {
            CalDigi *cal = ev->addCalDigi();
            CalXtalId::CalTrigMode mode = CalXtalId::BESTRANGE;
            Short_t tower = 5;
            Short_t layer = 4;
            Short_t col = 3;
            CalXtalId xtalId(tower, layer, col);
            cal->initialize(mode, xtalId);
            Char_t rangeM = CalXtalId::LEX8;
            Char_t rangeP = CalXtalId::HEX8;
            UShort_t adcM = 4095;
            UShort_t adcP = 4095;
            cal->addReadout(rangeP, adcP, rangeM, adcM);            
        }

        UInt_t idigi;
        for (idigi = 0; idigi < numDigi; idigi++) {
            TkrDigi *tkr = new TkrDigi();
            Int_t tot[2] = {idigi, idigi+1};
            TowerId id(3, 2);
            tkr->initialize(idigi, GlastAxis::Y, id, tot);
            UInt_t istrip;
            for (istrip = 0; istrip < idigi; istrip++) {
                tkr->addC0Hit(istrip);
            }
            tkr->addC1Hit(idigi*2);
            ev->addTkrDigi(tkr);
        }

        for (idigi = 0; idigi < numAcd; idigi++) {
            AcdId id(0, 2, 5, 4);
            Float_t energy = ievent*randNum;
            Bool_t veto[2] = {kFALSE, kTRUE};
            Bool_t low[2] = {kTRUE, kTRUE};
            Bool_t high[2] = {kFALSE, kTRUE};
            UShort_t pha[2] = {4095, 1};
            VolumeIdentifier volId;
            volId.Clear();
            volId.append(1);
            ev->addAcdDigi(id, volId, energy, pha, veto, low, high);
        }


        int idiag;
        for (idiag = 0; idiag < numCalDiag; idiag++) {
            CalDiagnosticData *calDiag = ev->addCalDiagnostic();
            UInt_t dataWord = 10101;
            UInt_t tower = 2;
            UInt_t layer = 1;
            calDiag->initialize(dataWord,tower,layer);
        }
        for (idiag=0; idiag<numTkrDiag; idiag++) {
            TkrDiagnosticData *tkrDiag = ev->addTkrDiagnostic();
            UInt_t dataWord = 20301;
            UInt_t tower = 4;
            UShort_t gtcc = 2;
            tkrDiag->initialize(dataWord,tower,gtcc);
        }

        int item;
        for (item = 0; item < numTem; item++) {
            Tem* curTem = ev->addTem();
            ErrorData err(1, 2, true, false);
            curTem->init(item, err);
        }

        t->Fill();
        ev->Clear();
    }
    
    std::cout << "Filled ROOT file with " << numEvents << " events" << std::endl;
    delete ev;
    
    f->Write();
    f->Close();
    delete f;
    return(0);
}


/// Main program
/// Return 0 for success.
/// Returns -1 for failure.
int main(int argc, char **argv) {
    char *fileName = "digi.root";
    int n = 1;
    int numEvents = 10;
    if (argc > 1) {
        fileName = argv[n++];
    } 
    if (argc > 2) {
        numEvents = atoi(argv[n++]);
    } 
    
    int sc = 0;
    sc = write(fileName, numEvents);
    sc = read(fileName, numEvents);
    
    if (sc == 0) {
        std::cout << "DIGI ROOT file writing and reading succeeded!" << std::endl;
    } else {
        std::cout << "FAILED digi writing and reading" << std::endl;
    }
    
    return(sc);
}


