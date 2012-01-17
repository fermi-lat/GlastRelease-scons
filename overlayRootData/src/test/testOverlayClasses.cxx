#include <cstdlib>
#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCollection.h"
#include "overlayRootData/EventOverlay.h"

/** @file testOverlayClasses.cxx
* @brief This defines a test routine for the overlay ROOT classes.
*
* This program create a new overlay ROOT file, and the opens it up again
* for reading.  The contents are checked against the values known to be stored
* during the original writing of the file.
* The contents of the file are printed to the screen.
* The program returns 0 if the test passed.
* If failure, the program returns -1.
*
* $Header$
*/
const UInt_t   runNum = 1;
const Int_t    numCal = 10;
const UInt_t   numTkr = 12;
const UInt_t   numAcd = 11;
const Bool_t   fromMc = true;

const Double_t depEnergy = 100.;
const TVector3 hitPosition(350.,350.,500.);

Float_t        randNum;

bool floatInRange(Double_t actual, Double_t desired) 
{
    const Double_t fudge=0.00001;

    if ( (actual >= (desired - fudge)) && (actual <= (desired+fudge)) ) return true;

    return false;
}

int checkEventOverlay(EventOverlay *evt, UInt_t ievent) 
{
    if (evt->getRunId() != runNum) 
    {
        std::cout << "Run Id is wrong: " << evt->getRunId() << std::endl;
        return -1;
    }

    if (evt->getEventId() != ievent) 
    {
        std::cout << "Event Id is wrong: " << evt->getEventId() << std::endl;
        return -1;
    }

    if (!floatInRange(evt->getTimeStamp(), randNum*ievent) ) 
    {
        std::cout << "Time Stamp is wrong: " << evt->getTimeStamp() << std::endl;
        return -1;
    }

    if (!floatInRange(evt->getLiveTime(), randNum*ievent) ) 
    {
        std::cout << "LiveTime is wrong: " << evt->getLiveTime() << std::endl;
        return -1;
    }

    if (evt->getFromMc() != fromMc) 
    {
        std::cout << "From MC flag is wrong" << std::endl;
        return -1;
    }
    return 0;
}

int checkCalOverlay(CalOverlay *digi, Int_t ievent) 
{
    // Checks the contents of one CalOverlay object

    CalXtalId id = digi->getPackedId();
    if ( (id.getTower() != ievent) || (id.getLayer() != 4) || (id.getColumn() != 3) ) 
    {
        std::cout << "CalOverlay id is wrong: " << id.getTower() << ", "
        << id.getLayer() << ", " << id.getColumn() << std::endl;
        return -1;
    }

    if (digi->getEnergy() != depEnergy)
    {
        std::cout << "CalOverlay energy is wrong " << digi->getEnergy() 
                  << ", " << depEnergy << std::endl;
        return -1;
    }

    if (digi->getPosition() != hitPosition)
    {
        std::cout << "CalOverlay position is wrong" << std::endl;
        return -1;
    }

    return 0;
}

int checkTkrOverlay(TkrOverlay *digi, UInt_t /*ievent*/, UInt_t idigi) 
{
    TowerId id = digi->getTower();
    if ((id.ix() != 3) || (id.iy() != 2) ) 
    {
        std::cout << "TkrOverlay id is wrong (x,y): (" << id.ix() << "," 
            << id.iy() << ")" << std::endl;
        return -1;
    }

    if (id.id() != 3+4*2)  
    {
        std::cout << "TkrOverlay id is computer wrong: " << id.id() << std::endl;
        return -1;
    }

    UShort_t bilayer = digi->getBilayer();
    if (bilayer != idigi) 
    {
        std::cout << "TkrOverlay bilayer is wrong:  " << bilayer << std::endl;
        return -1;
    }

    GlastAxis::axis a = digi->getView();
    if (a != GlastAxis::Y) 
    {
        std::cout << "TkrOverlay Axis is wrong: " << a << std::endl;
        return -1;
    }

    if (digi->getNumHits() != idigi+1) 
    {
        std::cout << "TkrOverlay number of hit strip is wrong: " << digi->getNumHits() << std::endl;
        return -1;
    }

    UInt_t istrip;
    for (istrip = 0; istrip < idigi; istrip++) 
    {
        Int_t strip = digi->getHit(istrip);
        if (strip != (Int_t)istrip) {
            std::cout << "TkrOverlay strip is wrong: " << strip << std::endl;
            return -1;
        }
    }

    Int_t strip = digi->getHit(idigi);
    if (strip != idigi*2) 
    {
        std::cout << "TkrOverlay strip is wrong: " << strip << std::endl;
        return -1;
    }

    strip = digi->getHit(idigi+1);
    if (strip != -1) 
    {
        std::cout << "TkrOverlay non-existant strip is defined: " << digi->getHit(idigi+1) << std::endl;
        return -1;
    }
    
    if (digi->getToT(0) != (Int_t)idigi) 
    {
        std::cout << "TkrOverlay ToT 0 is wrong: " << digi->getToT(0) << std::endl;
        return -1;
    }
    
    if (digi->getToT(1) != (Int_t)idigi+1) 
    {
        std::cout << "TkrOverlay ToT 1 is wrong: " << digi->getToT(1) << std::endl;
        return -1;
    }
    
    if (digi->getToT(2) != -1) 
    {
        std::cout << "TkrOverlay non-existant ToT is defined: " << digi->getToT(2) << std::endl;
        return -1;
    }
    
    return 0;
}

int checkGem(const GemOverlay &gem) 
{
    GemOverlayTileList list = gem.getTileList();

    if (list.getXzm()                  != 20)  return -1;
    if (list.getXzp()                  != 21)  return -1;
    if (list.getYzm()                  != 22)  return -1;
    if (list.getYzp()                  != 23)  return -1;
    if (list.getXy()                   != 24)  return -1; 
    if (list.getRbn()                  != 25)  return -1;
    if (list.getNa()                   != 26)  return -1;

    GemOverlayOnePpsTime time = gem.getOnePpsTime();
    if (time.getTimebase()             != 50)  return -1;
    if (time.getSeconds()              != 100) return -1;

    if (gem.getTkrVector()             != 1)   return -1;
    if (gem.getRoiVector()             != 2)   return -1;
    if (gem.getCalLeVector()           != 3)   return -1;
    if(gem.getCalHeVector()            != 4)   return -1;
    if (gem.getCnoVector()             != 5)   return -1;
    if (gem.getConditionSummary()      != 6)   return -1;
    if (gem.getMissed()                != 1)   return -1;
    if (gem.getLiveTime()              != 7)   return -1;
    if (gem.getPrescaled()             != 8)   return -1;
    if (gem.getDiscarded()             != 9)   return -1;
    if (gem.getCondArrTime().condArr() != 10)  return -1;
    if (gem.getTriggerTime()           != 11)  return -1;
    if (gem.getDeltaEventTime()        != 12)  return -1;
    if (gem.getDeltaWindowOpenTime()   != 13)  return -1;

    return 0; 
}

int checkAcdOverlay(AcdOverlay *digi, UInt_t ievent, UInt_t /*idigi*/) 
{
    // HMK Unused? Float_t f = ievent;

    if (!floatInRange(digi->getEnergyDep(), ievent*randNum) ) 
    {
        std::cout << "AcdOverlay Energy is wrong " << digi->getEnergyDep() << std::endl;
        return -1;
    }

    AcdId id = digi->getAcdId();
    if (id.getLayer() != 0) 
    {
        std::cout << "AcdId Layer is wrong " << id.getLayer() << std::endl;
        return -1;
    }

    VolumeIdentifier volId = digi->getVolId();
    if ( (volId.getBits0to31() != 0) || (volId.getBits32to63() != 16777216) || (volId.size() != 1) ) 
    {
        std::cout << "AcdOverlay VolId is incorrect" << std::endl;
        return -1;
    }

    if (id.getFace() != 2) 
    {
        std::cout << "AcdId Face is wrong " << id.getFace() << std::endl;
        return -1;
    }

    if (id.getRow() != 5) 
    {
        std::cout << "AcdId Row is wrong " << id.getRow() << std::endl;
        return -1;
    }

    if (id.getColumn() != 4) 
    {
        std::cout << "AcdId Col is wrong " << id.getColumn() << std::endl;
        return -1;
    }

    return 0;
}


/// Read in the ROOT file just generated via the write method
int read(char* fileName, int numEvents) 
{
    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("Overlay");
    EventOverlay *evt = 0;
    t->SetBranchAddress("EventOverlay", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    for (Int_t ievent = 0; ievent < numEvents; ievent++) 
    {
        t->GetEvent(ievent);
        std::cout << "EventOverlay ievent = " << ievent << std::endl;
        evt->Print();
        if (checkEventOverlay(evt, ievent) < 0) return -1;

        GemOverlay gemData = evt->getGemOverlay();
        if (checkGem(gemData) < 0) return -1;

        /*UInt_t numCalOverlay =*/ evt->getCalOverlayCol()->GetEntries();
        const TObjArray *CalOverlayCol = evt->getCalOverlayCol();
        TIter CalOverlayIt(CalOverlayCol);
        CalOverlay *digi = 0;
        while((digi = (CalOverlay*)CalOverlayIt.Next())!=0) 
        {
            digi->Print();
            if (checkCalOverlay(digi, ievent) < 0) return -1;
        }

        const TObjArray *TkrOverlayCol = evt->getTkrOverlayCol();
        TIter TkrOverlayIt(TkrOverlayCol);
        TkrOverlay *tDigi = 0;
        UInt_t idigi = 0;
        while ((tDigi = (TkrOverlay*)TkrOverlayIt.Next())!=0) 
        {
            tDigi->Print();
            if (checkTkrOverlay(tDigi, ievent, idigi) < 0) return -1;
            idigi++;
        }

        const TObjArray *AcdOverlayCol = evt->getAcdOverlayCol();
        TIter AcdOverlayIt(AcdOverlayCol);
        AcdOverlay *aDigi = 0;
        idigi = 0;
        while ((aDigi = (AcdOverlay*)AcdOverlayIt.Next())!=0) 
        {
            aDigi->Print();
            if (checkAcdOverlay(aDigi, ievent, idigi) < 0) return -1;
            idigi++;
        }
    }
    
    f->Close();
    delete f;
    
    return 0;
}

/// Create a new Monte Carlo ROOT file
int write(char* fileName, int numEvents) 
{
    Int_t buffer = 64000;
    Int_t splitLevel = 1;
    
    TFile *f =  new TFile(fileName, "RECREATE");
    TTree *t = new TTree("Overlay", "Digi");
    EventOverlay *ev = new EventOverlay();
    t->Branch("EventOverlay", "EventOverlay", &ev, buffer, splitLevel);
    
    std::cout << "Created new ROOT file" << std::endl;
    
    TRandom randGen;
    Int_t ievent;
    randNum = randGen.Rndm();

    for (ievent = 0; ievent < numEvents; ievent++) 
    {
        ev->initialize(ievent, runNum, randNum*ievent, randNum*ievent, fromMc);

        GemOverlay myGem;
        GemOverlayTileList tileList(20, 21, 22, 23, 24, 25, 26);
        myGem.initTrigger(1, 2, 3, 4, 5, 6, 1, tileList);
        GemOverlayOnePpsTime ppsTime(50, 100);
        myGem.initSummary(7, 8, 9, 10, 11, ppsTime, 12, 13);
        ev->setGemOverlay(myGem);
        
        for (int idx = 0; idx < numCal; idx++) 
        {
            CalOverlay *cal = new CalOverlay();

            Short_t tower = idx;
            Short_t layer = 4;
            Short_t col   = 3;
            CalXtalId xtalId(tower, layer, col);

            Double_t depEnergy = 100.;

            cal->initialize(xtalId, hitPosition, depEnergy, 0);
        }

        for (UInt_t idx = 0; idx < numTkr; idx++) 
        {
            TkrOverlay *tkr = new TkrOverlay();
            Int_t tot[2] = {idx, idx+1};
            TowerId id(3, 2);
            tkr->initialize(idx, GlastAxis::Y, id, tot);

            for (UInt_t istrip = 0; istrip < idx; istrip++) 
            {
                tkr->addC0Hit(istrip);
            }

            tkr->addC1Hit(idx*2);
            ev->addTkrOverlay(tkr);
        }

        for (UInt_t idx = 0; idx < numAcd; idx++) 
        {
            AcdId id(0, 2, 5, 4);
            VolumeIdentifier volId;
            volId.Clear();
            volId.append(1);

            Float_t energy = ievent*randNum;

            AcdOverlay* acd = new AcdOverlay(volId, id, energy, hitPosition);

            ev->addAcdOverlay(acd);
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
int main(int argc, char **argv) 
{
    char *fileName = "overlay.root";
    int n = 1;
    int numEvents = 10;
    if (argc > 1) 
    {
        fileName = argv[n++];
    } 
    if (argc > 2) 
    {
        numEvents = atoi(argv[n++]);
    } 
    
    int sc = 0;
    sc = write(fileName, numEvents);
    sc = read(fileName, numEvents);
    
    if (sc == 0) 
    {
        std::cout << "OVERLAY ROOT file writing and reading succeeded!" << std::endl;
    } 
    else 
    {
        std::cout << "FAILED overlay writing and reading" << std::endl;
    }
    
    return(sc);
}


