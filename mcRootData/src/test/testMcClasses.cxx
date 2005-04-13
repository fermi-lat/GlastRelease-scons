#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "mcRootData/McEvent.h"

/** @file testMcClasses.cxx
* @brief This defines a test routine for the Monte Carlo ROOT classes.
*
* This program create a new Monte Carlo ROOT file, and the opens it up again
* for reading.  The contents are checked against the values known to be stored
* during the original writing of the file.
* The contents of the file are printed to the screen.
* The program returns 0 if the test passed.
* If failure, the program returns -1.
*
* $Header$
*/
const UInt_t runNum = 1;
const UInt_t numParticles = 5;
const Int_t sourceId = 7;
const UInt_t sequence = 4;
Float_t randNum;

bool floatInRange(Double_t actual, Double_t desired) {
  const Double_t fudge=0.00001;
  if ( (actual >= (desired - fudge)) && (actual <= (desired+fudge)) ){
    return true;
  }
  return false;

}

int checkMcEvent(const McEvent* mcEvt, UInt_t ievent) {
    // Purpose and Method: Check the contents of the McEvent read in from the 
    //    ROOT file.

    if (mcEvt->getEventId() != ievent) {
        std::cout << "EventId does not match!!" << std::endl;
        return -1;
    }
    if (mcEvt->getRunId() != runNum) {
        std::cout << "Run Number does not match!!" << std::endl;
        return -1;
    }
	if (mcEvt->getSourceId() != sourceId) {
		std::cout << "Source Id does not match!!" << std::endl;
		return -1;
	}
	if (mcEvt->getSequence() != sequence) {
		std::cout << "Sequence does not match!!" << std::endl;
		return -1;
	}
    
    return 0;
}

int checkMcParticle(McParticle* mcPart, Int_t ipart, UInt_t ievent) {
    // Purpose and Method: Check the contents of the McParticle read in from 
    //   the ROOT file
    Float_t f = Float_t(ipart);
    Float_t fr = f*randNum;
    
    if (mcPart->getParticleId() != ipart) {
        std::cout << "McParticle Id does not match!!" << std::endl;
        return -1;
    }
    if (mcPart->getStatusFlags() != 0) {
        std::cout << "McParticle flags do not match!!" << std::endl;
        return -1;
    }
    TLorentzVector initMom = mcPart->getInitialFourMomentum();
    if ( !floatInRange(initMom.X(), f) || !floatInRange(initMom.Y(), f) || 
        !floatInRange(initMom.Z(), f) || !floatInRange(initMom.T(), f) ) {
        std::cout << "McParticle initial momentum does not match" << std::endl;
        return -1;
    }
    TLorentzVector finalMom = mcPart->getFinalFourMomentum();
    if ( !floatInRange(finalMom.X(), fr+ievent) || 
         !floatInRange(finalMom.Y(), fr+ievent) || 
        !floatInRange(finalMom.Z(), fr+ievent) || 
        !floatInRange(finalMom.T(), fr+ievent) ) {
        std::cout << "McParticle final momentum: (" << finalMom.X() << ","
            << finalMom.Y() << "," << finalMom.Z() << "," << finalMom.T() << ")" << std::endl;
        return -1;
    }

    TVector3 initPos = mcPart->getInitialPosition();
    if ( !floatInRange(initPos.X(), randNum) ||
         !floatInRange(initPos.Y(), 2.0*fr) ||
         !floatInRange(initPos.Z(), 4.0*fr) ) {
        std::cout << "McParticle initial position does not match" << std::endl;
        return -1;
    }

    TVector3 finalPos = mcPart->getFinalPosition();
    if ( !floatInRange(finalPos.X(), fr) || 
         !floatInRange(finalPos.Y(), fr) || 
        !floatInRange(finalPos.Z(), fr) ) {
        std::cout << "McParticle final position does not match" << std::endl;
        return -1;
    }
    
    if (mcPart->getMother()) {
        if (mcPart->getMother()->GetUniqueID() != mcPart->GetUniqueID()) {
            std::cout << "McParticle mother is not itself!" << std::endl;
            return -1;
        }
    } else {
        std::cout << "Cannot check TRefs using ROOT 3.02.07" << std::endl;
    }

    if (mcPart->getProcess() != "") {
        std::cout << "Process String is wrong: " << mcPart->getProcess() << std::endl;
        return -1;
    }
    
    return 0;
}

int checkMcPositionHit(const McPositionHit* mcPosHit, Int_t ipart, UInt_t ievent) {
    // Purpose and Method: Check the contents of the McPositionHit read in from 
    //   the ROOT file
    Float_t f = Float_t(ipart);
    Float_t fr = f*randNum;
    
    if (mcPosHit->getMcParticleId() != 7) {
        std::cout << "MC Particle Id is wrong: " << mcPosHit->getMcParticleId() << std::endl;
        return -1;
    }

    if (mcPosHit->getOriginMcParticleId() != -13) {
        std::cout << "MC Origin Id is wrong: " << mcPosHit->getOriginMcParticleId() << std::endl;
        return -1;
    }

    commonRootData::VolumeIdentifier id = mcPosHit->getVolumeId();
    std::cout << "McPosHit Id = " << id.name() << std::endl;
    if ( (id.getBits0to31() != 0) || (id.getBits32to63() != 16777216) || (id.size() != 1) ) {
        std::cout << "McPosHit VolId is incorrect" << std::endl;
        return -1;
    }
    
    TVector3 entry = mcPosHit->getEntryPosition();
    if  ( !floatInRange(entry.X(), 1.) || 
          !floatInRange(entry.Y(), 1.) ||
          !floatInRange(entry.Z(), 1.) ) {
        std::cout << "McPosHit entry is (" << entry.X() << "," << entry.Y() << ","
            << entry.Z() << ")" << std::endl;
        return -1;
    }
    
    TVector3 exit = mcPosHit->getExitPosition();
    if ( !floatInRange(exit.X(), fr) || !floatInRange(exit.Y(), fr) || 
         !floatInRange(exit.Z(), fr) ) {
        std::cout << "McPosHit exit is (" << exit.X() << "," << exit.Y() 
            << "," << exit.Z() << std::endl;
        return -1;
    }

    TVector3 gEntry = mcPosHit->getGlobalEntryPosition();
    if ( !floatInRange(gEntry.X(), 3.) || !floatInRange(gEntry.Y(), 3.) || 
         !floatInRange(gEntry.Z(), 3.) ) {
        std::cout << "McPosHit Global Entry is (" << gEntry.X() << "," << gEntry.Y() 
            << "," << gEntry.Z() << std::endl;
        return -1;
    }

    TVector3 gExit = mcPosHit->getGlobalExitPosition();
    if ( !floatInRange(gExit.X(), fr*2.) || !floatInRange(gExit.Y(), fr*2.) || 
         !floatInRange(gExit.Z(), fr*2.) ) {
        std::cout << "McPosHit Global exit is (" << gExit.X() << "," << gExit.Y() 
            << "," << gExit.Z() << std::endl;
        return -1;
    }
    
    if (!floatInRange(mcPosHit->getDepositedEnergy(), randNum)) {
        std::cout << "Error:  McPosHit dep Energy: " << mcPosHit->getDepositedEnergy() 
            << std::endl;
        return -1;
    }
    
    if (!floatInRange(mcPosHit->getParticleEnergy(), randNum*0.1)) {
        std::cout << "Error: McPosHit particle Energy = " << mcPosHit->getParticleEnergy() 
            << std::endl;
        
        return -1;
    }
    
    if (!floatInRange(mcPosHit->getTimeOfFlight(), randNum*0.4)) {
        std::cout << "Error: McPosHit TOF: " << mcPosHit->getTimeOfFlight() << std::endl;
        return -1;
    }
    
    return 0;
}

int checkMcIntegratingHit(McIntegratingHit* mcIntHit, UInt_t ipart, 
                          UInt_t ievent) {
    // Purpose and Method: Check the contents of the McIntegratingHit read in 
    //    from the ROOT file
    commonRootData::VolumeIdentifier id = mcIntHit->getVolumeId();
    std::cout << "McIntHit volumeId: " << id.name() << std::endl;
    if ( (id.size() != 1) || (id.getBits0to31() != 0) || (id.getBits32to63() != 0) ) {
        std::cout << "McIntHit VolId is incorrect" << std::endl;
        return -1;
    }

    mcIntHit->itemizedEnergyReset();
    const McParticle *myPart;
    double energy;
    if (mcIntHit->itemizedEnergySize() != 1) {
        std::cout << "McIntHit map size: " << mcIntHit->itemizedEnergySize() << std::endl;
        return -1;
    }
    int count = 0;
    while ( (myPart = mcIntHit->itemizedEnergyNext(energy)) ) {
        ++count;
        if (!floatInRange(energy, 1.5)) {
            std::cout << "McInt map energy is incorrect: " << energy << std::endl;
            return -1;
        }
    }
    
    energy = mcIntHit->getMcParticleEnergy(McIntegratingHit::PRIMARY);
    if (!floatInRange(energy, 2.5)) {
        std::cout << "Id energy map is incorrect: " << energy << std::endl;
        return -1;
    }

    energy = mcIntHit->getMcParticleEnergy(McIntegratingHit::ELECTRON);
    if ( !floatInRange(energy, 3.0) ) {
        std::cout << "Id energy map for e- is wrong: " << energy << std::endl;
        return -1;
    }

    energy = mcIntHit->getMcParticleEnergy(McIntegratingHit::POSITRON);
    if ( !floatInRange(energy, 0.) ) {
        std::cout << "Id energy map for e+ is wrong: " << energy << std::endl;
        return -1;
    }

    double totE = mcIntHit->getTotalEnergy();
    if ( !floatInRange(totE, 5.5)) {
        std::cout << "Total E is incorrect: " << totE << std::endl;
        return -1;
    }

    TVector3 moment1 = mcIntHit->getMoment1();

    TVector3 moment2 = mcIntHit->getMoment2();

    if (count == 0) {
        std::cout << "No Error cannot read map which stores TRefs using" 
            << " compiled code and ROOT 3.02.07" << std::endl;
    }
    
    return 0;
}


int read(char* fileName, unsigned int numEvents) {
    // Purpose and Method: Read in the ROOT file just generated via the 
    //   write method

    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("Mc");
    McEvent *evt = 0;
    t->SetBranchAddress("McEvent", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    UInt_t ievent, ipart;
    for (ievent = 0; ievent < numEvents; ievent++) {
        t->GetEvent(ievent);
        std::cout << "McEvent ievent = " << ievent << std::endl;
        evt->Print();
        if (checkMcEvent(evt, ievent) < 0) return -1;
        UInt_t numPart = evt->getMcParticleCol()->GetEntries();
        for (ipart = 0; ipart < numPart; ipart++) {
            evt->getMcParticle(ipart)->Print();
            if (checkMcParticle(evt->getMcParticle(ipart), ipart, ievent) < 0) return -1;
        }
        numPart = evt->getMcPositionHitCol()->GetEntries();
        for (ipart = 0; ipart < numPart; ipart++) {
            evt->getMcPositionHit(ipart)->Print();
            if (checkMcPositionHit(evt->getMcPositionHit(ipart), ipart, ievent) < 0) return -1;
        }
        numPart = evt->getMcIntegratingHitCol()->GetEntries();
        for (ipart = 0; ipart < numPart; ipart++){
            evt->getMcIntegratingHit(ipart)->Print();
            if (checkMcIntegratingHit(evt->getMcIntegratingHit(ipart), ipart, ievent) < 0) return -1;
        }
    }
    
    f->Close();
    delete f;
    
    return(0);
}

int write(char* fileName, UInt_t numEvents) {
    // Purpose and Method: Create a new Monte Carlo ROOT file

    Int_t buffer = 64000;
    Int_t splitLevel = 1;
    
    TFile *f =  new TFile(fileName, "RECREATE");
    TTree *t = new TTree("Mc", "Mc");
    McEvent *ev = new McEvent();
    t->Branch("McEvent", "McEvent", &ev, buffer, splitLevel);
    
    std::cout << "Created new ROOT file" << std::endl;
    
    TRandom randGen;
    commonRootData::VolumeIdentifier id;
    UInt_t ievent, ipart;
    randNum = randGen.Rndm();
    for (ievent = 0; ievent < numEvents; ievent++) {
        
        ev->initialize(ievent, runNum, sourceId, sequence, ievent*1.0);
        
        for (ipart = 0; ipart < numParticles; ipart ++) {
            McParticle *mcPart = new McParticle();
            Float_t f = Float_t(ipart);
            Float_t fr = f*randNum;
            TLorentzVector initMom(f, f, f, f);
            TLorentzVector finalMom(fr+ievent, fr+ievent, fr+ievent, fr+ievent);
            TVector3 initPos(randNum, f*2.0*randNum, f*4.0*randNum);
            TVector3 finalPos(f*randNum, f*randNum, f*randNum);
            mcPart->initialize(mcPart, ipart, 0, initMom, finalMom, initPos, finalPos);
            ev->addMcParticle(mcPart);
            
            McPositionHit *posHit = new McPositionHit();
            TVector3 entry(1., 1., 1.);
            TVector3 exit(fr, fr, fr);
            TVector3 gEntry(3., 3., 3.);
            TVector3 gExit(fr*2., fr*2., fr*2.);
            id.Clear();
            id.append(1);
            Double_t depE = randNum;
            Double_t partE = randNum*0.1;
            Double_t tof = randNum*0.4;
            UInt_t flags = 0;
            McParticle *originMcPart = 0;
            Int_t particleId = 7;
            Int_t originId = -13;
            posHit->initialize(particleId, originId, depE, id, 
                entry, exit, gEntry, gExit,
                mcPart, originMcPart, partE, tof, flags);
            ev->addMcPositionHit(posHit);
            
            McIntegratingHit *intHit = new McIntegratingHit();
            id.Clear();
            id.append(0);
            intHit->initialize(id);
            TVector3 pos = mcPart->getFinalPosition();
            TVector3 pos2(1.3, 0.0, 12.0);
            intHit->addEnergyItem(1.5, mcPart, pos);
            double totE = 5.5;
            double energyArr[3] = { 2.5, 3.0, 0.0 };
            TVector3 moment1(1.0, 2.0, 3.0);
            TVector3 moment2(2.0, 4.0, 6.0);
            intHit->setEnergyItems(totE, energyArr, moment1, moment2);
            ev->addMcIntegratingHit(intHit);
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


int main(int argc, char **argv) {
    // Purpose and Method: Main program
    // Output:  Return 0 for success. Returns -1 for failure.

    char *fileName = "mc.root";
    int n = 1;
    unsigned int numEvents = 10;
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
      std::cout << "MC ROOT file writing and reading suceeded!" << std::endl;
    } else {
      std::cout << "FAILED" << std::endl;
    }
    
    return(sc);
}


