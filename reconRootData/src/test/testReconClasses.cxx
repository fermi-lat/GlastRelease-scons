#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCollection.h"
#include "reconRootData/ReconEvent.h"
#include <iostream>

/** @file testReconClasses.cxx
* @brief This defines a test routine for the Reconstruction ROOT classes.
*
* This program creates a new Recon ROOT file, and the opens it up again
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
const UInt_t numClusters = 20;
Float_t randNum;

bool floatInRange(Double_t actual, Double_t desired) {
    const Double_t fudge=0.00001;
    if ( (actual >= (desired - fudge)) && (actual <= (desired+fudge)) ){
        return true;
    }
    return false;
    
}

int checkReconEvent(ReconEvent *evt, UInt_t ievent) {
    if (evt->getRunId() != runNum) {
        std::cout << "Run Id is wrong: " << evt->getRunId() << std::endl;
        return -1;
    }
    if (evt->getEventId() != ievent) {
        std::cout << "Event Id is wrong: " << evt->getEventId() << std::endl;
        return -1;
    }
    
    return 0;
}

int checkCalCluster(const CalCluster* cluster, UInt_t ievent) {
    Float_t f = (Float_t)ievent;
    Float_t fr = f*randNum;
    
    (*cluster).Print();
    
    if (!floatInRange((*cluster).getEnergyLeak(), randNum)) {
        std::cout << "Energy Leak is : " << (*cluster).getEnergyLeak() << std::endl;
        return -1;
    }
    if (!floatInRange((*cluster).getRmsLong(), randNum*f) ) {
        std::cout << "RMS Long is: " << (*cluster).getRmsLong() << std::endl;
        return -1;
    }
    if (!floatInRange((*cluster).getRmsTrans(), randNum*f*2.) ) {
        std::cout << "RMS Trans is: " << (*cluster).getRmsTrans() << std::endl;
        return -1;
    }
    if (!floatInRange((*cluster).getTransvOffset(), randNum*f*3.) ) {
        std::cout << "TransvOffset is: " << (*cluster).getTransvOffset() << std::endl;
        return -1;
    }
    TVector3 dir = (*cluster).getDirection();
    if ( (!floatInRange(dir.X(), randNum)) || (!floatInRange(dir.Y(), randNum*3.))
        || (!floatInRange(dir.Z(), randNum*5.) ) ){
        std::cout << "Dir: ( " << dir.X() << "," << dir.Y() << "," 
            << dir.Z() << ")" << std::endl;
        return -1;
    }
    
    std::vector<Double_t> eLayer = (*cluster).getEneLayer();
    if (eLayer.size() != 2) {
        std::cout << "E Layer vector is wrong size " 
            << eLayer.size() << std::endl;
        return -1;
    }
    if (!floatInRange(eLayer[0], 15.0)) {
        std::cout << "1st energy Layer item: " << eLayer[0] << std::endl;
        return -1;
        
    }
    if (!floatInRange(eLayer[1], 22.0)) {
        std::cout << "2nd energy Layer item: " << eLayer[1] << std::endl;
        return -1;
    }
    
    
    std::vector<TVector3> pLayer = (*cluster).getPosLayer();
    if (pLayer.size() != 2) {
        std::cout << "POS Layer vector is wrong size " 
            << pLayer.size() << std::endl;
        return -1;
    }
    TVector3 pos0 = pLayer[0];
    TVector3 pos1 = pLayer[1];
    if ( (!floatInRange(pos0.X(), f)) || (!floatInRange(pos0.Y(), f))
        || (!floatInRange(pos0.Z(), f) ) ){
        std::cout << "1st pos Layer item: ( " << pos0.X() << ","
            << pos0.Y() << "," << pos0.Z() << ")" << std::endl;
        return -1;
        
    }
    if ( (!floatInRange(pos1.X(), fr)) || (!floatInRange(pos1.Y(), fr*2.))
        || (!floatInRange(pos1.Z(), fr*3.) ) ) {
        std::cout << "2nd pos Layer item: ( " << pos1.X() << ","
            << pos1.Y() << "," << pos1.Z() << ")" << std::endl;
        return -1;
    }
    
    std::vector<TVector3> rLayer = (*cluster).getRmsLayer();
    if (rLayer.size() != 3) {
        std::cout << "RMS Layer vector is wrong size " 
            << rLayer.size() << std::endl;
        return -1;
    }
    TVector3 rms0 = rLayer[0];
    TVector3 rms1 = rLayer[1];
    TVector3 rms2 = rLayer[2];
    if ( (!floatInRange(rms0.X(), randNum)) || (!floatInRange(rms0.Y(), randNum*5.))
        || (!floatInRange(rms0.Z(), randNum*10.) ) ) {
        std::cout << "1st rms Layer item: ( " << rms0.X() << ","
            << rms0.Y() << "," << rms0.Z() << ")" << std::endl;
        return -1;
        
    }
    if ( (!floatInRange(rms1.X(), f+1.)) || (!floatInRange(rms1.Y(), f+2.))
        || (!floatInRange(rms1.Z(), f+3.) ) ) {
        std::cout << "2nd rms Layer item: ( " << rms1.X() << ","
            << rms1.Y() << "," << rms1.Z() << ")" << std::endl;
        return -1;
    }
    if ( (!floatInRange(rms2.X(), 14.)) || (!floatInRange(rms2.Y(), 22.))
        || (!floatInRange(rms2.Z(), 44.) ) ){
        std::cout << "3rd rms Layer item: ( " << rms2.X() << ","
            << rms2.Y() << "," << rms2.Z() << ")" << std::endl;
        return -1;
    }
    
    
    if(!floatInRange((*cluster).getFitEnergy(), f) ){
        std::cout << "Fit energy: " << (*cluster).getFitEnergy() << std::endl;
        return -1;
    }
    if (!floatInRange((*cluster).getProfChisq(), fr) ) {
        std::cout << "Chi2: " << (*cluster).getProfChisq() << std::endl;
        return -1;
    }
    if( !floatInRange((*cluster).getCsiAlpha(), f*f)) {
        std::cout << "CsiAlpha: " << (*cluster).getCsiAlpha() << std::endl;
        return -1;
    }
    if( !floatInRange((*cluster).getCsiLambda(), f+f)){
        std::cout << "CsiLambda: " << (*cluster).getCsiLambda() << std::endl;
        return -1;
    }
    if( !floatInRange((*cluster).getCsiStart(), randNum*randNum)) {
        std::cout << "CsiAlpha: " << (*cluster).getCsiStart() << std::endl;
        return -1;
    }
    
    return 0;
}

int checkCalXtalRec(const CalXtalRecData *rec, UInt_t ievent) {
    
    Float_t f = Float_t(ievent);
    Float_t fr = f*randNum;
    
    (*rec).Print();
    
    if ((*rec).getMode() != CalXtalId::BESTRANGE) {
        std::cout << "Xtal mode is not BESTRANGE" << std::endl;
        return -1;
    }
    CalXtalId id = (*rec).getPackedId();
    if ( (id.getTower() != 1) || (id.getLayer() != 2) || (id.getColumn() != 3) ) {
        std::cout << "Xtal Id is wrong (tower, Layer, Col): (" << id.getTower
            << "," << id.getLayer() << "," << id.getColumn() << ")" << std::endl;
        return -1;
    }
    
    TVector3 pos = (*rec).getPosition();
    if ( (!floatInRange(pos.X(), 4.5) ) || (!floatInRange(pos.Y(), 7.5) )
        || (!floatInRange(pos.Z(), 8.5) ) ) {
        std::cout << "Xtal pos is (" << pos.X() << "," << pos.Y()
            << "," << pos.Z() << ")" << std::endl;
        return -1;
    }
    
    Char_t rangeP = (*rec).getRange(0, CalXtalId::POS);
    Char_t rangeM = (*rec).getRange(0, CalXtalId::NEG);
    
    if (rangeP != CalXtalId::LEX8) {
        std::cout << "POS range: " << rangeP << std::endl;
        return -1;
    }
    if (rangeM != CalXtalId::HEX8) {
        std::cout << "NEG range: " << rangeM << std::endl;
        return -1;
    }
    Double_t energy = (*rec).getEnergy();
    Double_t energyP = (*rec).getEnergy(0, CalXtalId::POS);
    Double_t energyM = (*rec).getEnergy(0, CalXtalId::NEG);
    if (!floatInRange(energyP, fr) ) {
        std::cout << "Range POS energy: " << energyP << std::endl;
        return -1;
    }
    if (!floatInRange(energyM, randNum*4.) ) {
        std::cout << "Range NEG energy: " << energyM << std::endl;
        return -1;
    }
    
    if ( !floatInRange(energy, (energyP+energyM)/2.) ) {
        std::cout << "energy of RangeRecData is not average: " << energy << std::endl;
        return -1;
    }
    
    const CalRangeRecData* rangeRecData = (*rec).getRangeRecData(0);
    if (!floatInRange(energyP, rangeRecData->getEnergy(CalXtalId::POS)) ) {
        std::cout << "CalRangeRecData POS energy differs " << std::endl;
        return -1;
    }
    if (!floatInRange(energyM, rangeRecData->getEnergy(CalXtalId::NEG)) ) {
        std::cout << "CalRangeRecData NEG energy differs " << std::endl;
        return -1;
    }
    if (CalXtalId::LEX8 != rangeRecData->getRange(CalXtalId::POS)) {
        std::cout << "CalRangeRecData POS range differs" << std::endl;
        return -1;
    }
    if (CalXtalId::HEX8 != rangeRecData->getRange(CalXtalId::NEG)) {
        std::cout << "CalRangeRecData NEG range differs" << std::endl;
        return -1;
    }
    
    Double_t energyRangeP = (*rec).getEnergySelectedRange(CalXtalId::LEX8, CalXtalId::POS);
    if (!floatInRange(energyRangeP, energyP)) {
        std::cout << "get selected POS range differs" << std::endl;
        return -1;
    }
    Double_t energyRangeM = (*rec).getEnergySelectedRange(CalXtalId::HEX8, CalXtalId::NEG);
    if (!floatInRange(energyRangeM, energyM) ) {
        std::cout << "get selected NEG range differs" << std::endl;
        return -1;
    }
    
    // Check range and faces that should return -1, since they are undefined
    if ( !floatInRange((*rec).getEnergySelectedRange(CalXtalId::LEX1, CalXtalId::POS), -1.) ) {
        std::cout << "undefined range face pair does not return -1" << std::endl;
        return -1;
    }
    if ( !floatInRange((*rec).getEnergySelectedRange(CalXtalId::HEX1, CalXtalId::POS), -1.) ) {
        std::cout << "undefined range face pair does not return -1" << std::endl;
        return -1;
    }
    if ( !floatInRange((*rec).getEnergySelectedRange(CalXtalId::HEX8, CalXtalId::POS), -1.) ) {
        std::cout << "undefined range face pair does not return -1" << std::endl;
        return -1;
    }
    
    if ( !floatInRange((*rec).getEnergySelectedRange(CalXtalId::LEX1, CalXtalId::NEG), -1.) ) {
        std::cout << "undefined range face pair does not return -1" << std::endl;
        return -1;
    }
    if ( !floatInRange((*rec).getEnergySelectedRange(CalXtalId::LEX8, CalXtalId::NEG), -1.) ) {
        std::cout << "undefined range face pair does not return -1" << std::endl;
        return -1;
    }
    if ( !floatInRange((*rec).getEnergySelectedRange(CalXtalId::HEX1, CalXtalId::NEG), -1.) ) {
        std::cout << "undefined range face pair does not return -1" << std::endl;
        return -1;
    }
    
    return 0;
}

int checkCalRecon(CalRecon *cal, UInt_t ievent) {
    // Checks the contents of one CalRecon object
    
    Float_t f = (Float_t)ievent;
    Float_t fr = f*randNum;
    
    std::vector<CalCluster> clusterCol = cal->getCalClusterCol();
    if (clusterCol.size() != numClusters) {
        std::cout << "Number of clusters in collection is wrong" << std::endl;
        return -1;
    }
    std::vector<CalCluster>::const_iterator clusterIt;
    for (clusterIt = clusterCol.begin(); clusterIt != clusterCol.end(); clusterIt++) {
        if (checkCalCluster(clusterIt, ievent) < 0) return -1;
    }
    
    std::vector<CalXtalRecData> rec = cal->getCalXtalRecCol();
    if (rec.size() != numXtals) {
        std::cout << "Number of CalXtalRecData objects is wrong" << std::endl;
        return -1;
    }
    
    std::vector<CalXtalRecData> recCol = cal->getCalXtalRecCol();
    std::vector<CalXtalRecData>::const_iterator recIt;
    for (recIt = recCol.begin(); recIt != recCol.end(); recIt++) {
        if (checkCalXtalRec(recIt, ievent) < 0) return -1;
    }
    
    return 0;
}

//int checkTkrRecon(TkrRecon *tkr, UInt_t ievent) {

//    return 0;
//}


/// Read in the ROOT file just generated via the write method
int read(char* fileName, int numEvents) {
    TFile *f = new TFile(fileName, "READ");
    TTree *t = (TTree*)f->Get("Recon");
    ReconEvent *evt = 0;
    t->SetBranchAddress("ReconEvent", &evt);
    
    std::cout << "Opened the ROOT file for reading" << std::endl;
    
    UInt_t ievent;
    for (ievent = 0; ievent < numEvents; ievent++) {
        t->GetEvent(ievent);
        std::cout << "ReconEvent ievent = " << ievent << std::endl;
        evt->Print();
        if (checkReconEvent(evt, ievent) < 0) return -1;
        CalRecon *cal = evt->getCalRecon();
        cal->Print();
        if (checkCalRecon(cal, ievent) < 0) return -1;
        //TkrRecon *tkr = evt->getTkrRecon();
        //if (checkTkrRecon(tkr, ievent) < 0) return -1;
        evt->Clear();
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
    TTree *t = new TTree("Recon", "Recon");
    ReconEvent *ev = new ReconEvent();
    t->Branch("ReconEvent", "ReconEvent", &ev, buffer, splitLevel);
    
    std::cout << "Created new ROOT file" << std::endl;
    
    TRandom randGen;
    Int_t ievent, ixtal;
    randNum = randGen.Rndm();
    for (ievent = 0; ievent < numEvents; ievent++) {
        
        Float_t f = Float_t(ievent);
        CalRecon *calRec = new CalRecon();
        
        UInt_t icluster;
        for (icluster = 0; icluster < numClusters; icluster++ ) {
            CalCluster cluster;
            std::vector<Double_t> eLayer;
            eLayer.push_back(15.0);
            eLayer.push_back(22.0);
            std::vector<TVector3> pLayer;
            pLayer.push_back(TVector3(f, f, f));
            pLayer.push_back(TVector3(randNum*f, randNum*2.*f, randNum*3.*f));
            std::vector<TVector3> rmsLayer;
            rmsLayer.push_back(TVector3(randNum, randNum*5, randNum*10));
            rmsLayer.push_back(TVector3(f+1., f+2., f+3.));
            rmsLayer.push_back(TVector3(14., 22., 44.));
            TVector3 calDir(randNum, randNum*3, randNum*5);
            Double_t eLeak = randNum;
            Double_t rmsLong = randNum*f;
            Double_t rmsTrans = randNum*2.*f;
            Double_t transOffset = randNum*3.*f;
            cluster.initialize(eLeak, eLayer, pLayer, rmsLayer, rmsLong, 
                rmsTrans, calDir, transOffset);
            
            Double_t fitEnergy = f;
            Double_t chi2 = f*randNum;
            Double_t fStart = randNum*randNum;
            Double_t fitAlpha = f*f;
            Double_t fitLambda = f+f;
            cluster.initProfile(fitEnergy, chi2, fStart, fitAlpha, fitLambda);
            
            calRec->addCalCluster(cluster);
        }
        
        for (ixtal = 0; ixtal < numXtals; ixtal ++) {
            CalXtalRecData xtal;
            CalXtalId id;
            id.init(1, 2, 3);
            xtal.initialize(CalXtalId::BESTRANGE, id);
            CalRangeRecData rec(CalXtalId::LEX8, randNum*f, CalXtalId::HEX8, randNum*4.0);
            TVector3 pos(4.5, 7.5, 8.5);
            rec.initialize(pos);
            xtal.addRangeRecData(rec);
            calRec->addXtalRecData(xtal);
        }

        TkrRecon *tkrRec = new TkrRecon();
        ev->initialize(ievent, runNum, tkrRec, calRec);
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
    char *fileName = "recon.root";
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
        std::cout << "RECON ROOT file writing and reading succeeded!" << std::endl;
    } else {
        std::cout << "FAILED recon writing and reading" << std::endl;
    }
    
    return(sc);
}


