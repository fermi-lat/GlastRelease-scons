#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCollection.h"
#include "TVector3.h"
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
const UInt_t numTracks = 15;
const UInt_t numVertices = 11;
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
        std::cout << "Xtal Id is wrong (tower, Layer, Col): (" << id.getTower()
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
    
    TObjArray *clusterCol = cal->getCalClusterCol();
    if (clusterCol->GetEntries() != numClusters) {
        std::cout << "Number of clusters in collection is wrong" << std::endl;
        return -1;
    }
    TIter clusterIt(clusterCol);
    CalCluster *cluster;
    while ( (cluster = (CalCluster*)clusterIt.Next()) ) {
        if (checkCalCluster(cluster, ievent) < 0) return -1;
    }
    
    TObjArray *recCol = cal->getCalXtalRecCol();
    if (recCol->GetEntries() != numXtals) {
        std::cout << "Number of CalXtalRecData objects is wrong" << std::endl;
        return -1;
    }
    
    TIter recIt(recCol);
    CalXtalRecData *rec;
    while ( (rec = (CalXtalRecData*) recIt.Next()) ) {
        if (checkCalXtalRec(rec, ievent) < 0) return -1;
    }
    
    return 0;
}

int checkTkrSiCluster(const TkrSiCluster *cluster, UInt_t ievent, UInt_t icluster) {
    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;
    if ( cluster->getId() != icluster)  {
        std::cout << "SiCluster id is wrong: " << cluster->getId() << std::endl;
        return -1;
    }
    if ( cluster->getLayer() != 5 ) {
        std::cout << "SiCluster layer is wrong: " << cluster->getLayer() << std::endl;
        return -1;
    }
    if ( cluster->getXY() != TkrSiCluster::X) {
        std::cout << "SiCluster axes is wrong: " << cluster->getXY() << std::endl;
        return -1;
    }
    if ( cluster->getCenterStrip() != 7)  {
        std::cout << "SiCluster center is wrong: " << cluster->getCenterStrip() << std::endl;
        return -1;
    }
    if ( cluster->getNumStrips() != 6) {
        std::cout << "SiCluster numStrips is wrong: " << cluster->getNumStrips() << std::endl;
        return -1;
    }
    if ( !floatInRange((*cluster).getPosition(), f)) {
        std::cout << "SiCluster Position: " << cluster->getPosition() << std::endl;
        return -1;
    }
    if ( !floatInRange(cluster->getZPosition(), fr)) {
        std::cout << "SiCluster Z Position: " << cluster->getZPosition() << std::endl;
        return -1;
    }

    return 0;
}

int checkCandTrack(const TkrCandTrack *track, UInt_t ievent, UInt_t iHit) {

    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;

    if ((*track).getId() != iHit) {
        std::cout << "TkrCandTrack id is: " << track->getId() << std::endl;
        return -1;
    }
    if (track->getLayer() != 3) {
        std::cout << "TkrCandTrack Layer: " << track->getLayer() << std::endl;
        return -1;
    }
    if (track->getTower() != 1) {
        std::cout << "TkrCandTrack Tower: " << track->getTower() << std::endl;
        return -1;
    }

    if (!floatInRange(track->getQuality(), f)) {
        std::cout << "TkrCandTrack Quality: " << track->getQuality() << std::endl;
        return -1;
    }

    if (!floatInRange(track->getEnergy(), fr)) {
        std::cout << "TkrCandTrack energy: " << track->getEnergy() << std::endl;
        return -1;
    }

    TVector3 pos = track->getPosition();
    if ( (!floatInRange(pos.X(), f)) || (!floatInRange(pos.Y(), fr))
        || (!floatInRange(pos.Z(), f)) ) {
        std::cout << "TkrCandTrack pos (x,y,z): (" << pos.X() << ","
            << pos.Y() << "," << pos.Z() << ")" << std::endl;
        return -1;
    }
    TVector3 dir = track->getDirection();
    if ( (!floatInRange(dir.X(), fr)) || (!floatInRange(dir.Y(), f))
        || (!floatInRange(dir.Z(), fr)) ) {
        std::cout << "TkrCandTrack dir (x,y,z): (" << dir.X() << ","
            << dir.Y() << "," << dir.Z() << ")" << std::endl;
        return -1;
    }

    return 0;
}

int checkTrack(const TkrTrack *track,  UInt_t ievent, UInt_t itrack) {

    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;
    if (track->getId() != itrack) {
        std::cout << "Track Id: " << track->getId() << std::endl;
        return -1;
    }
    if (track->getXgaps() != itrack*itrack) {
        std::cout << "Track xgaps: " << track->getXgaps() << std::endl;
        return -1;
    }
    if (track->getYgaps() != 2*itrack) {
        std::cout << "Track ygaps: " << track->getYgaps() << std::endl;
        return -1;
    }
    if (track->getXistGaps() != ievent+itrack) {
        std::cout << "Track xistgaps: " << track->getYistGaps() << std::endl;
        return -1;
    }
    if (track->getYistGaps() != 2*ievent) {
        std::cout << "Track yistgaps: " << track->getYistGaps() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getChiSq(), f) ) {
        std::cout << "Track chisq is wrong: " << track->getChiSq() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getChiSqSmooth(), fr) ) {
        std::cout << "Track chisq smooth is wrong: " << track->getChiSqSmooth() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getKalEnergy(), f*2.0) ) {
        std::cout << "Track KalEnergy is wrong: " << track->getKalEnergy() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getKalThetaMS(), f*6.0) ) {
        std::cout << "Track kalThetaMs is wrong: " << track->getKalThetaMS() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getQuality(), f*f) ) {
        std::cout << "Track quality is wrong: " << track->getQuality() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getRmsResid(), randNum*randNum) ) {
        std::cout << "Track rms is wrong: " << track->getRmsResid() << std::endl;
        return -1;
    }
    return 0;
}

int checkTkrVertex(const TkrVertex* vertex, UInt_t ievent, UInt_t ivertex) {
    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;

    if (vertex->getId() != ivertex) {
        std::cout << "Vertex Id is wrong: " << vertex->getId() << std::endl;
        return -1;
    }
    if (!floatInRange(vertex->getEnergy(), fr) ) {
        std::cout << "Vertex energy is wrong: " << vertex->getEnergy() << std::endl;
        return -1;
    }
    if (!floatInRange(vertex->getQuality(), f) ) {
        std::cout << "Vertex quality is wrong: " << vertex->getQuality() << std::endl;
        return -1;
    }

    if (vertex->getLayer() != 2) {
        std::cout << "Vertex layer is wrong: " << vertex->getLayer() << std::endl;
        return -1;
    }

    if (vertex->getTower() != 4) {
        std::cout << "Vertex tower is wrong: " << vertex->getTower() << std::endl;
        return -1;
    }

    TVector3 pos = vertex->getPosition();
    if ( (!floatInRange(pos.X(), f)) || (!floatInRange(pos.Y(), fr))
        || (!floatInRange(pos.Z(), fr*2.)) ) {
        std::cout << "TkrCandTrack pos (x,y,z): (" << pos.X() << ","
            << pos.Y() << "," << pos.Z() << ")" << std::endl;
        return -1;
    }
    TVector3 dir = vertex->getDirection();
    if ( (!floatInRange(dir.X(), randNum*2.)) || (!floatInRange(dir.Y(), randNum*3.))
        || (!floatInRange(dir.Z(), randNum*4.)) ) {
        std::cout << "TkrCandTrack dir (x,y,z): (" << dir.X() << ","
            << dir.Y() << "," << dir.Z() << ")" << std::endl;
        return -1;
    }

    TkrParams params = vertex->getTrackPar();
    if ( (!floatInRange(params.getXPos(), f) ) ||
    (!floatInRange(params.getXSlope(), f*f) ) ){
        std::cout << "Vertex X Param wrong: (" << params.getXPos() << ","
            << params.getXSlope() << ")" << std::endl;
        return -1;
    }
    
    if ( (!floatInRange(params.getYPos(), fr) ) ||
    (!floatInRange(params.getYSlope(), randNum*randNum) ) ) {
        std::cout << "Vertex Y Param wrong: (" << params.getYPos() << ","
            << params.getYSlope() << ")" << std::endl;
        return -1;
    }


    TkrCovMat mat = vertex->getTrackCov();
    if ( (!floatInRange(mat.getCovX0X0(), f) )||
        (!floatInRange(mat.getCovSxSx(), f*2.) )||
        (!floatInRange(mat.getCovX0Sx(), f*3.) ) ) {
        std::cout << "Vertex matrix x vals are wrong: " << mat.getCovX0X0() << " "
            << mat.getCovSxSx() << " " << mat.getCovX0Sx() << std::endl;
        //return -1;
    }

    if ( (!floatInRange(mat.getCovY0Y0(), f*4.) ) ||
        (!floatInRange(mat.getCovY0Sy(), f*5.) ) ||
        (!floatInRange(mat.getCovSySy(), f*6.) ) ) {
        std::cout << "Vertex matrix y vals are wrong: " << mat.getCovY0Y0() << " "
            << mat.getCovY0Sy() << " " << mat.getCovSySy() << std::endl;
        //return -1;

    }


    if (vertex->getNumTracks() != 2) {
        std::cout << "Vertex number of tracks is wrong: " << vertex->getNumTracks() << std::endl;
        return -1;
    }
    if (vertex->getTrackId(0) != ivertex+1) {
        std::cout << "Vertex Track id is wrong: " << vertex->getTrackId(0) << std::endl;
        return -1;
    }

    if (vertex->getTrackId(1) != ivertex+2) {
        std::cout << "Vertex Track id is wrong: " << vertex->getTrackId(1) << std::endl;
        return -1;
    }

    return 0;
}

int checkTkrRecon(TkrRecon *tkr, UInt_t ievent) {

    TObjArray *siClusterCol = tkr->getSiClusterCol();
    if (siClusterCol->GetEntries() != numClusters) {
        std::cout << "Number of TkrRecon clusters is wrong: " << siClusterCol->GetEntries() << std::endl;
        return -1;
    }
    UInt_t icluster = 0;
    TIter siClusIt(siClusterCol);
    TkrSiCluster *cluster;
    while ( (cluster = (TkrSiCluster*)(siClusIt.Next()) ) ){
        if (checkTkrSiCluster(cluster, ievent, icluster) < 0) return -1;
        icluster++;
    }

    TObjArray *candTrackCol = tkr->getTrackCandCol();
    if (candTrackCol->GetEntries() != numTracks) {
        std::cout << "Number of TkrRecon cand tracks is wrong: " << candTrackCol->GetEntries() << std::endl;
        return -1;
    }
    TIter candTrackIt(candTrackCol);
    UInt_t iHit = 0;
    TkrCandTrack *candTrack;
    while ( (candTrack = (TkrCandTrack*)(candTrackIt.Next()) ) ) {
        if (checkCandTrack(candTrack, ievent, iHit) < 0)
            return -1;
        iHit++;
    }

    UInt_t iTrack = 0;
    TObjArray* trackCol = tkr->getTrackCol();
    if (trackCol->GetEntries() != numTracks) {
        std::cout << "Number of TkrRecon tracks is wrong: " << trackCol->GetEntries() << std::endl;
        return -1;
    }
    TIter trackIt(trackCol);
    TkrTrack *track;
    while ((track = (TkrTrack*)(trackIt.Next()) )) {
        if (checkTrack(track, ievent, iTrack) < 0) 
            return -1;
        iTrack++;
    }

    TObjArray* vertexCol = tkr->getVertexCol();
    if (vertexCol->GetEntries() != numVertices) {
        std::cout << "Number of TkrRecon vertices is wrong: " << vertexCol->GetEntries() << std::endl;
        return -1;
    }
    TIter vertexIt(vertexCol);
    TkrVertex *vertex;
    UInt_t ivertex = 0;
    while ((vertex = (TkrVertex*)(vertexIt.Next()) )){
        vertex->Print();
        if (checkTkrVertex(vertex, ievent, ivertex) < 0) 
            return -1;
        ivertex++;
    }

    return 0;
}


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
        if (!cal) return -1;
        cal->Print();
        if (checkCalRecon(cal, ievent) < 0) return -1;
        TkrRecon *tkr = evt->getTkrRecon();
        if (!tkr) return -1;
        tkr->Print();
        if (checkTkrRecon(tkr, ievent) < 0) return -1;
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

        // Create CalRecon object
        CalRecon *calRec = new CalRecon();
        calRec->initialize();
        UInt_t icluster;
        for (icluster = 0; icluster < numClusters; icluster++ ) {
            CalCluster *cluster = new CalCluster();
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
            cluster->initialize(eLeak, eLayer, pLayer, rmsLayer, rmsLong, 
                rmsTrans, calDir, transOffset);
            
            Double_t fitEnergy = f;
            Double_t chi2 = f*randNum;
            Double_t fStart = randNum*randNum;
            Double_t fitAlpha = f*f;
            Double_t fitLambda = f+f;
            cluster->initProfile(fitEnergy, chi2, fStart, fitAlpha, fitLambda);
            
            calRec->addCalCluster(cluster);
        }
        
        for (ixtal = 0; ixtal < numXtals; ixtal ++) {
            CalXtalRecData *xtal = new CalXtalRecData();
            CalXtalId id;
            id.init(1, 2, 3);
            xtal->initialize(CalXtalId::BESTRANGE, id);
            CalRangeRecData rec(CalXtalId::LEX8, randNum*f, CalXtalId::HEX8, randNum*4.0);
            TVector3 pos(4.5, 7.5, 8.5);
            rec.initialize(pos);
            xtal->addRangeRecData(rec);
            calRec->addXtalRecData(xtal);
        }

        // Create TkrRecon object
        TkrRecon *tkrRec = new TkrRecon();
        tkrRec->initialize();

        TkrSiCluster *siCluster;
        for (icluster = 0; icluster < numClusters; icluster++ ) {
            siCluster = new TkrSiCluster();
            UInt_t layer = 5;
            TkrSiCluster::TKRAxes xy = TkrSiCluster::X;
            UShort_t center = 7;
            UShort_t nStrips = 6;
            Float_t pos = f;
            Float_t zpos = f*randNum;
            siCluster->initialize(icluster, layer, xy, center, nStrips, pos, zpos);
            tkrRec->addSiCluster(siCluster);
        }
        UInt_t itrack;
        TkrCandTrack *candTrack;
        for (itrack = 0; itrack < numTracks; itrack++) {
            candTrack = new TkrCandTrack();
            UInt_t lay = 3;
            UInt_t tow = 1;
            TVector3 pos(f, f*randNum, f);
            TVector3 dir(f*randNum, f, f*randNum);
            Double_t quality = f;
            Double_t e = f*randNum;
            candTrack->initialize(itrack, lay, tow, quality, e, pos, dir);
            tkrRec->addTrackCand(candTrack);
        }
        TkrTrack *track;
        for (itrack=0; itrack < numTracks; itrack++) {
            track = new TkrTrack();
            UInt_t xgaps = itrack * itrack;
            UInt_t ygaps = itrack+itrack;
            UInt_t x1st = ievent+itrack;
            UInt_t y1st = ievent+ievent;
            track->initializeInfo(itrack, xgaps, ygaps, x1st, y1st);
            Double_t chiSq = f;
            Double_t chiSqSmooth = f*randNum;
            Double_t rms = randNum*randNum;
            Double_t qual = f*f;
            Double_t e = f*2.0;
            Double_t ms = f*6.0;
            track->initializeQual(chiSq, chiSqSmooth, rms, qual, e, ms);
            tkrRec->addTrack(track);
        }
        
        UInt_t ivertex;
        TkrVertex *vertex;
        for (ivertex=0; ivertex < numVertices; ivertex++) {
            vertex = new TkrVertex();
            UInt_t layer = 2;
            UInt_t tower = 4;
            Double_t qual = f;
            Double_t energy = f*randNum;
            vertex->initializeInfo(ivertex, layer, tower, qual, energy);
            
            TkrParams vtxPar;
            Double_t ax = f;
            Double_t sx = f*f;
            Double_t ay = f*randNum;
            Double_t sy = randNum*randNum;
            vtxPar.initialize(ax, sx, ay, sy);
            
            TkrCovMat vtxCov;
            Double_t a = f;
            Double_t b = f*2.;
            Double_t c = f*3.;
            Double_t d = f*4.;
            Double_t e = f*5.;
            Double_t g = f*6.;
            vtxCov.initialize(a, b, c, d, e, g);

            TVector3 pos(f, f*randNum, f*randNum*2.);
            TVector3 dir(randNum*2., randNum*3., randNum*4.);
            vertex->addTrackId(ivertex+1);
            vertex->addTrackId(ivertex+2);
            vertex->initializeVals(vtxPar, vtxCov, pos, dir);
            tkrRec->addVertex(vertex);
        }


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


