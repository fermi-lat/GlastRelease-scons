#include <iostream>
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
const UInt_t numXtals = 10;
const Bool_t fromMc = true;
Float_t randNum;

int checkDigiEvent(DigiEvent *evt, UInt_t ievent) {
    if (evt->getRunId() != runNum) {
        std::cout << "Run Id is wrong: " << evt->getRunId() << std::endl;
        return -1;
    }
    if (evt->getEventId() != ievent) {
        std::cout << "Event Id is wrong: " << evt->getEventId() << std::endl;
        return -1;
    }

    if (evt->getFromMc() != fromMc) {
        std::cout << "From MC flag is wrong" << std::endl;
        return -1;
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
    const TClonesArray *readoutCol = digi->getReadoutCol();
    UInt_t numEntries = readoutCol->GetEntries();
    if (numEntries != 1) {
        std::cout << "Wrong number of readouts in Cal Digi " 
            << numEntries << std::endl;
        return -1;
    }
    TIter readoutIt(readoutCol);
    CalXtalReadout *calReadout;
    while(calReadout = (CalXtalReadout*) readoutIt.Next()) {
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

/// Read in the ROOT file just generated via the write method
int read(char* fileName, int numEvents) {
    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("Digi");
    DigiEvent *evt = 0;
    t->SetBranchAddress("DigiEvent", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    UInt_t ievent;
    for (ievent = 0; ievent < numEvents; ievent++) {
        t->GetEvent(ievent);
        std::cout << "DigiEvent ievent = " << ievent << std::endl;
        evt->Print();
        if (checkDigiEvent(evt, ievent) < 0) return -1;
        UInt_t numCalDigi = evt->getCalDigiCol()->GetEntries();
        const TObjArray *calDigiCol = evt->getCalDigiCol();
        TIter calDigiIt(calDigiCol);
        CalDigi *digi = 0;
        while(digi = (CalDigi*)calDigiIt.Next()) {
            digi->Print();
            if (checkCalDigi(digi, ievent) < 0) return -1;
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
        
        ev->initialize(ievent, runNum, fromMc);
        
        for (ixtal = 0; ixtal < numXtals; ixtal ++) {
            CalDigi *cal = new CalDigi();
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
            ev->addCalDigi(cal);
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


