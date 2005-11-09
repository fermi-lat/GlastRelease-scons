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
const UInt_t numEventEnergies = 2 ;
const UInt_t numClusters = 20;
const UInt_t numMipTracks = 5;
const UInt_t numTracks = 15;
const UInt_t numVertices = 11;
const UInt_t numClusLayerData = 5;
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

int checkCalEventEnergy(const CalEventEnergy* energy, UInt_t rank, UInt_t  /*ievent*/) {

    CalEventEnergy energyRef ;
    energyRef.Fake(rank,randNum) ;
    if (energy->CompareInRange(energyRef))
      return 0 ;
    else
      return -1 ;
}

int checkCalCluster(const CalCluster* cluster, UInt_t rank, UInt_t /*ievent*/ ) {

    CalCluster clusterRef ;
    clusterRef.Fake(rank,randNum) ;
    if (cluster->CompareInRange(clusterRef))
      return 0 ;
    else
      return -1 ;
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

int checkCalMipTrack( const CalMipTrack *mip, Int_t /*ievent*/ ) {
    const TVector3 pos = mip->getPoint();

    const TVector3 dir = mip->getDir();

    if (!floatInRange(mip->getD2C(), 1.0) ) {
        std::cout << "D2C is not 1." << std::endl;
        return -1;
    }
    if (!floatInRange(mip->getD2Edge(), 2.) ){
        std::cout << "D2Edge is not 5." << std::endl;
        return -1;
    }
    if (mip->getCalEdge() != 3) {
        std::cout << "CalEdge is not 3" << std::endl;
        return -1;
    }
    if (!floatInRange(mip->getArcLen(), 4.0) ) { 
        std::cout << "ArcLength is not 4" << std::endl;
        return -1;
    }
    if (!floatInRange(mip->getEcor(), 5.0) ) { 
        std::cout << "Ecor is not 5" << std::endl;
        return -1;
    }
    if (!floatInRange(mip->getEcorRms(), 6.0) ) { 
        std::cout << "EcorRms is not 6" << std::endl;
        return -1;
    }
    if (!floatInRange(mip->getChi2(), 7.0) ) { 
        std::cout << "chi2 is not 7" << std::endl;
        return -1;
    }
    if (!floatInRange(mip->getErm(), 8.0) ) { 
        std::cout << "erm is not 8" << std::endl;
        return -1;
    }
    return 0;

}

int checkCalRecon(CalRecon *cal, UInt_t ievent) {
    // Checks the contents of one CalRecon object

    UInt_t rank ;

    TObjArray *eventEnergyCol = cal->getCalEventEnergyCol();
    if (UInt_t(eventEnergyCol->GetEntries()) != numEventEnergies) {
        std::cout << "Number of event energies in collection is wrong" << std::endl;
        return -1;
    }
    rank = 0 ;
    TIter eventEnergyIt(eventEnergyCol);
    CalEventEnergy * eventEnergy ;
    while ( (eventEnergy = (CalEventEnergy*)eventEnergyIt.Next()) ) {
        if (checkCalEventEnergy(eventEnergy,rank++,ievent)<0) return -1 ;
    }

    TObjArray *clusterCol = cal->getCalClusterCol();
    if (UInt_t(clusterCol->GetEntries()) != numClusters) {
        std::cout << "Number of clusters in collection is wrong" << std::endl;
        return -1;
    }
    rank = 0 ;
    TIter clusterIt(clusterCol);
    CalCluster *cluster;
    while ( (cluster = (CalCluster*)clusterIt.Next()) ) {
        if (checkCalCluster(cluster,rank++,ievent) < 0) return -1;
    }

    TObjArray *recCol = cal->getCalXtalRecCol();
    if (UInt_t(recCol->GetEntries()) != numXtals) {
        std::cout << "Number of CalXtalRecData objects is wrong" << std::endl;
        return -1;
    }

    TIter recIt(recCol);
    CalXtalRecData *rec;
    while ( (rec = (CalXtalRecData*) recIt.Next()) ) {
        if (checkCalXtalRec(rec, ievent) < 0) return -1;
    }

    TObjArray *mipCol = cal->getCalMipTrackCol();
    if (UInt_t(mipCol->GetEntries()) != numMipTracks) {
        std::cout << "Number of MipTracks in collection is wrong" << std::endl;
        return -1;
    } 
    TIter mipIt(mipCol);
    CalMipTrack *mip = 0;
    while ( (mip = (CalMipTrack*) mipIt.Next()) ) {
        if (checkCalMipTrack(mip, ievent) < 0) return -1;
    }

    return 0;
}

int checkTkrCluster(const TkrCluster *cluster, UInt_t ievent, UInt_t icluster) {
    Float_t f = Float_t (ievent);
    //Float_t fr = f*randNum;

    commonRootData::TkrId tkrId = cluster->getTkrId();

    // This may go away in a future release
    //if ( cluster->getId() != icluster)  {
    //    std::cout << "TkrCluster id is wrong: " << cluster->getId() << std::endl;
    //    return -1;
    //}

    if ( tkrId.getView() != commonRootData::TkrId::eMeasureX) {
        std::cout << "TkrCluster axes is wrong: " << tkrId.getView() << std::endl;
        return -1;
    }
    if ( !floatInRange(cluster->getStrip(), 7.))  {
        std::cout << "TkrCluster center is wrong: " << cluster->getStrip() << std::endl;
        return -1;
    }
    if ( !floatInRange(cluster->getSize(), 7.)) {
        std::cout << "TkrCluster size is wrong: " << cluster->getSize() << std::endl;
        return -1;
    }

    if (tkrId.getTowerX() != icluster % 4) {
        std::cout << "TkrCluster Tower X number is wrong: " << tkrId.getTowerX() << std::endl;
        return -1;
    }

    if (tkrId.getTowerY() != (icluster / 4) % 4) {
        std::cout << "TkrCluster Tower Y number is wrong: " << tkrId.getTowerY() << std::endl;
        return -1;
    }

    if (cluster->getRawToT() != 50) {
        std::cout << "TkrCluster rawToT is wrong: " << cluster->getRawToT() << std::endl;
        return -1;
    }

    if ( !floatInRange(cluster->getMips(), f) ) {
        std::cout << "TkrCluster ToT is wrong: " << cluster->getMips() << std::endl;
        return -1;
    }
    TVector3 pos = cluster->getPosition();
    if ( (!floatInRange(pos.X(), f)) || (!floatInRange(pos.Y(), f*2.)) 
        || (!floatInRange(pos.Z(), f*3.)) ) {
            std::cout << "TkrCluster Position: (" << pos.X() << ", "
                << pos.Y() << ", " << pos.Z() << ")" << std::endl;
            return -1;
        }
        if ( !cluster->hitFlagged() ) {
            std::cout << "TkrCluster not flagged " << cluster->hitFlagged() << std::endl;
            return -1;
        }

        return 0;
}

int checkTrack(const TkrTrack *track,  UInt_t ievent, UInt_t itrack) {

    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;
    if (track->getNumXGaps() != itrack*itrack) {
        std::cout << "Track xgaps: " << track->getNumXGaps() << std::endl;
        return -1;
    }
    if (track->getNumYGaps() != 2*itrack) {
        std::cout << "Track ygaps: " << track->getNumYGaps() << std::endl;
        return -1;
    }
    if (track->getNumXFirstGaps() != ievent+itrack) {
        std::cout << "Track xistgaps: " << track->getNumXFirstGaps() << std::endl;
        return -1;
    }
    if (track->getNumYFirstGaps() != 2*ievent) {
        std::cout << "Track yistgaps: " << track->getNumYFirstGaps() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getChiSquareFilter(), f) ) {
        std::cout << "Track chisq is wrong: " << track->getChiSquareFilter() << std::endl;
        return -1;
    }
    if ( !floatInRange(track->getChiSquareSmooth(), fr) ) {
        std::cout << "Track chisq smooth is wrong: " << track->getChiSquareSmooth() << std::endl;
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
    if ( !floatInRange(track->getScatter(), randNum*randNum) ) {
        std::cout << "Track rms is wrong: " << track->getScatter() << std::endl;
        return -1;
    }
    return 0;
}


int checkTkrVertex(const TkrVertex* vertex, UInt_t ievent, UInt_t ivertex) {
    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;

    commonRootData::TkrId tkrId = vertex->getTkrId();

    if (!floatInRange(vertex->getEnergy(), fr) ) {
        std::cout << "Vertex energy is wrong: " << vertex->getEnergy() << std::endl;
        return -1;
    }
    if (!floatInRange(vertex->getQuality(), f) ) {
        std::cout << "Vertex quality is wrong: " << vertex->getQuality() << std::endl;
        return -1;
    }

    if (tkrId.getTray() != ivertex % 15) {
        std::cout << "Vertex layer is wrong: " << tkrId.getTray() << std::endl;
        return -1;
    }

    if (tkrId.getTowerX() != ivertex % 4) {
        std::cout << "Vertex tower X is wrong: " << tkrId.getTowerX() << std::endl;
        return -1;
    }

    if (tkrId.getTowerY() != ivertex / 4) {
        std::cout << "Vertex tower Y is wrong: " << tkrId.getTowerY() << std::endl;
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

            TkrTrackParams params = vertex->getVertexParams();
            if ( (!floatInRange(params.getxPosition(), f) ) ||
                (!floatInRange(params.getxSlope(), f*f) ) ){
                    std::cout << "Vertex X Param wrong: (" << params.getxPosition() << ","
                        << params.getxSlope() << ")" << std::endl;
                    std::cout << "f = " << f << " f*f = " << f*f << std::endl;
                    return -1;
                }

                if ( (!floatInRange(params.getyPosition(), fr) ) ||
                    (!floatInRange(params.getySlope(), randNum*randNum) ) ) {
                        std::cout << "Vertex Y Param wrong: (" << params.getyPosition() << ","
                            << params.getySlope() << ")" << std::endl;
                        return -1;
                    }

                    if ( (!floatInRange(params.getxPosxPos(), f) ) ||
                        (!floatInRange(params.getxPosxSlp(), f*2.) ) ||
                        (!floatInRange(params.getxPosyPos(), f*3.) ) || 
                        (!floatInRange(params.getxPosySlp(), f*4.)) ) 
                    {
                        std::cout << "Vertex matrix row 0 vals are wrong: " << params.getxPosxPos() << " "
                            << params.getxPosxSlp() << " " << params.getxPosyPos() 
                            << " " << params.getxPosySlp() << std::endl;
                        return -1;
                    }

                    if ( (!floatInRange(params(2,2), f*5.) ) ||
                        (!floatInRange(params(2,3), f*6.) ) || 
                        (!floatInRange(params(2,4), f*7.)) ) 
                    {
                        std::cout << "Vertex matrix row 1 vals are wrong: " << params(2,1) << " "
                            << params(2,2) << " " << params(2,3) << " " << params(2,4) << std::endl;
                        return -1;
                    }

                    if ( (!floatInRange(params(3,3), f*8.) ) || 
                        (!floatInRange(params(3,4), f*9.)) ) 
                    {
                        std::cout << "Vertex matrix row 2 vals are wrong: " << params(3,1) << " "
                            << params(3,2) << " " << params(3,3) << " " << params(3,4) << std::endl;
                        return -1;
                    }

                    if ( (!floatInRange(params(4,4), f*10.)) ) 
                    {
                        std::cout << "Vertex matrix row 3 vals are wrong: " << params(4,1) << " "
                            << params(4,2) << " " << params(4,3) << " " << params(4,4) << std::endl;
                        return -1;
                    }

                    //if (vertex->getNumTracks() != 2) {
                    //    std::cout << "Vertex number of tracks is wrong: " << vertex->getNumTracks() << std::endl;
                    //    return -1;
                    //}

                    // Hmm... what to do here...
                    //if (vertex->getTrackId(0) != ivertex+1) {
                    //    std::cout << "Vertex Track id is wrong: " << vertex->getTrackId(0) << std::endl;
                    //    return -1;
                    //}

                    //if (vertex->getTrackId(1) != ivertex+2) {
                    //    std::cout << "Vertex Track id is wrong: " << vertex->getTrackId(1) << std::endl;
                    //    return -1;
                    //}

                    return 0;
}

int checkTkrDiagnostics(const TkrDiagnostics* diag) {

    if (diag->getNumClusters() != 1) {
        std::cout << "tkrDiagnostics numClusters is not 1" << std::endl;
        return -1;
    }
    if (diag->getNumVecPoints() != 2) {
        std::cout << "tkrDiagnostics numVecPoints is not 2" << std::endl;
        return -1;
    }
    if (diag->getNumVecLinks() != 3) {
        std::cout << "tkrDiagnostics nuVecLinks is not 3" << std::endl;
        return -1;
    }
    if (diag->getnLinksNonZeroLayers() != 4) {
        std::cout << "tkrDiagnostics nLinksnonZeroLayers is not 4" << std::endl;
        return -1;
    }
    if (diag->getAveNumLinksLayer() != 5) {
        std::cout << "tkrDiagnostics AveNumLinksLayer is not 5 " << std::endl;
        return -1;
    }
    if ( !floatInRange(diag->getNumLinkCombinations(), 6.) ) {
        std::cout << "tkrDiagnostics NumLinknCombinations is not 6" << std::endl;
        return -1;
    }
    if (diag->getNumTrackElements() != 7) {
        std::cout <<"tkrdiagnostics numTrackElements is not 7" << std::endl;
        return -1;
    }
    if (diag->getNumTkrTracks() != 8) {
        std::cout << "tkrDiagnostics NumTkrTracks is not 8" << std::endl;
        return -1;
    }
 

    return 0;
}

int checkTkrRecon(TkrRecon *tkr, UInt_t ievent) {

    TObjArray *clusterCol = tkr->getClusterCol();
    if (UInt_t(clusterCol->GetEntries()) != numClusters) {
        std::cout << "Number of TkrRecon clusters is wrong: " << clusterCol->GetEntries() << std::endl;
        return -1;
    }
    UInt_t icluster = 0;
    TIter siClusIt(clusterCol);
    TkrCluster *cluster;
    while ( (cluster = (TkrCluster*)(siClusIt.Next()) ) ){
        if (checkTkrCluster(cluster, ievent, icluster) < 0) return -1;
        icluster++;
    }

    UInt_t iTrack = 0;
    TObjArray* trackCol = tkr->getTrackCol();
    if (UInt_t(trackCol->GetEntries()) != numTracks) {
        std::cout << "Number of TkrRecon tracks is wrong: " << trackCol->GetEntries() << std::endl;
        return -1;
    }
    TIter trackIt(trackCol);
    TObject *trackObj = 0;
    while ( TkrTrack* track = (TkrTrack*) trackIt.Next() ) 
    {
        if (checkTrack(track, ievent, iTrack) < 0) return -1;
        iTrack++;
    }

    TObjArray* vertexCol = tkr->getVertexCol();
    if (UInt_t(vertexCol->GetEntries()) != numVertices) {
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

    const TkrDiagnostics *diag = tkr->getDiagnostics();
    if (checkTkrDiagnostics(diag) < 0)
         return -1;

    return 0;
}

int checkAcdRecon(AcdRecon *acd, UInt_t ievent) {

    Float_t f = Float_t (ievent);
    Float_t fr = f*randNum;

    if (!floatInRange(acd->getEnergy(), f)) {
        std::cout << "AcdRecon energy is wrong: " << acd->getEnergy() << std::endl;
        return -1;
    }
    if (!floatInRange(acd->getRibbonEnergy(), f)) {
        std::cout << "AcdRecon ribbon energy is wrong: " << acd->getRibbonEnergy() << std::endl;
        return -1;
    }

    if (acd->getTileCount() != 5) {
        std::cout << "AcdRecon tile count is wrong: " << acd->getTileCount() << std::endl;
        return -1;
    }
    if (acd->getRibbonCount() != 2) {
        std::cout << "AcdRecon ribbon count is wrong: " << acd->getRibbonCount() << std::endl;
        return -1;
    }

    if (!floatInRange(acd->getGammaDoca(), fr)) {
        std::cout << "AcdRecon Gamma Doca is wrong: " << acd->getGammaDoca() << std::endl;
        return -1;
    }
    if (!floatInRange(acd->getCornerDoca(), fr)) {
        std::cout << "AcdRecon Corner Doca is wrong: " << acd->getCornerDoca() << std::endl;
        return -1;
    }

    if (!floatInRange(acd->getDoca(), randNum) ){
        std::cout << "AcdRecon Doca is wrong: " << acd->getDoca() << std::endl;
        return -1;
    }

    if (!floatInRange(acd->getActiveDist(), 2.*f) ) {
        std::cout << "AcdRecon Active Dist is wrong: " << acd->getActiveDist() << std::endl;
        return -1;
    }

    if (!floatInRange(acd->getRibbonActiveDist(), 3.*f) ) {
        std::cout << "AcdRecon Ribbon Active Dist is wrong: " << acd->getRibbonActiveDist() << std::endl;
        return -1;
    }

    AcdId ribActDistId = acd->getRibbonActDistId();
    if ((ribActDistId.getLayer() != 0) || (ribActDistId.getFace() != 1) ||
        (ribActDistId.getRow() != 4) || (ribActDistId.getColumn() != 3) ) {
            std::cout << "ribbonActDist AcdID is wrong: " 
                << ribActDistId.getLayer() << " "
                << ribActDistId.getFace() << " " << ribActDistId.getRow() << " "
                << ribActDistId.getColumn() << std::endl;
            return -1;
     }

    AcdId acdMinDocaId = acd->getMinDocaId();
    if ((acdMinDocaId.getLayer() != 0) || (acdMinDocaId.getFace() != 0) ||
        (acdMinDocaId.getRow() != 3) || (acdMinDocaId.getColumn() != 2) ) {
            std::cout << "MinDoca AcdID is wrong: " << acdMinDocaId.getLayer() << " "
                << acdMinDocaId.getFace() << " " << acdMinDocaId.getRow() << " "
                << acdMinDocaId.getColumn() << std::endl;
            return -1;
        }

  AcdId acdMaxActDistId = acd->getMaxActDistId();   
     if ((acdMaxActDistId.getLayer() != 0) || (acdMaxActDistId.getFace() != 0) ||   
         (acdMaxActDistId.getRow() != 2) || (acdMaxActDistId.getColumn() != 2) )   
     {   
             std::cout << "MaxActDist AcdID is wrong: "   
                 << acdMaxActDistId.getLayer() << " "   
                 << acdMaxActDistId.getFace() << " "   
                 << acdMaxActDistId.getRow() << " "   
                 << acdMaxActDistId.getColumn() << std::endl;   
             return -1;   
      } 


        std::vector<Double_t> rowCol = acd->getRowDocaCol();
        if (rowCol.size() != 2) {
            std::cout << "AcdRecon number of row entries is wrong: " << rowCol.size() << std::endl;
            return -1;
        }

        if (!floatInRange(rowCol[0], randNum) ) {
            std::cout << "AcdRecon row doca 0 is wrong: " << rowCol[0] << std::endl;
            return -1;
        }

        if (!floatInRange(rowCol[1], f) ) {
            std::cout << "AcdRecon row doca 1 is wrong: " << rowCol[1] << std::endl;
            return -1;
        }

        std::vector<Double_t> rowActDistCol = acd->getRowActDistCol();
        if (rowActDistCol.size() != 2) {
            std::cout << "AcdRecon number of row actDist entries is wrong: " << rowActDistCol.size() << std::endl;
            return -1;
        }

        if (!floatInRange(rowActDistCol[0], randNum) ) {
            std::cout << "AcdRecon row actdist 0 is wrong: " << rowActDistCol[0] << std::endl;
            return -1;
        }

        if (!floatInRange(rowActDistCol[1], f) ) {
            std::cout << "AcdRecon row actdist 1 is wrong: " << rowActDistCol[1] << std::endl;
            return -1;
        }

        std::vector<AcdId> idCol = acd->getIdCol();
        if (idCol.size() != 2) {
            std::cout << "AcdRecon number of ids is wrong: " << acd->getIdCol().size() << std::endl;
            return -1;
        }
        if ( !(idCol[0] == AcdId(0, 0, 3, 2)) ) {

        }

        if ( !(idCol[1] == AcdId(0, 2, 2, 1))) {

        }

        std::vector<Double_t> energyCol = acd->getEnergyCol();
        if (energyCol.size() != 2) {
            std::cout << "AcdRecon number of energies is wrong: " << acd->getEnergyCol().size() << std::endl;
            return -1;
        }
        if (!floatInRange(energyCol[0], f) ) {
            std::cout << "AcdRecon first energy is wrong: " << energyCol[0] << std::endl;
            return -1;
        }

        if (!floatInRange(energyCol[1], f*2.) ) {
            std::cout << "AcdRecon 2nd energy is wrong: " << energyCol[1] << std::endl;
            return -1;
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
    for (ievent = 0; ievent < (UInt_t) numEvents; ievent++) {
        t->GetEvent(ievent);
        std::cout << "ReconEvent ievent = " << ievent << std::endl;
        evt->Print();
        if (checkReconEvent(evt, ievent) < 0) return -1;
        AcdRecon *acd = evt->getAcdRecon();
        if (!acd) return -1;
        acd->Print();
        if (checkAcdRecon(acd, ievent) < 0) return -1;
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
    Int_t ievent;
    randNum = randGen.Rndm();
    for (ievent = 0; ievent < numEvents; ievent++) {

        Float_t f = Float_t(ievent);

        // Create AcdRecon object
        AcdRecon *acdRec = new AcdRecon();
        Double_t energy = f;
        Double_t ribbonE = energy;
        Int_t count = 5;
        Int_t ribbonCount = 2;
        Double_t gDoca = f*randNum;
        Double_t doca = randNum;
        Double_t actDist = 2.*f;
        Double_t ribActDist = 3.*f;
        AcdId minDocaId(0, 0, 3, 2);
        AcdId maxActDistId(0,0,2,2);
        AcdId ribActDistId(0,1,4,3);
        std::vector<Double_t> rowDocaCol;
        rowDocaCol.push_back(randNum);
        rowDocaCol.push_back(f);
        std::vector<Double_t> rowActDistCol;
        rowActDistCol.push_back(randNum);
        rowActDistCol.push_back(f);
        std::vector<AcdId> idCol;
        idCol.push_back(AcdId(0, 0, 3, 2));
        idCol.push_back(AcdId(0, 2, 2, 1));
        std::vector<Double_t> energyCol;
        energyCol.push_back(f);
        energyCol.push_back(2.*f);
        Double32_t cornerDoca = f*randNum;
        acdRec->initialize(energy, ribbonE, count, ribbonCount, gDoca, doca, 
            minDocaId,
            actDist, maxActDistId, ribActDist, ribActDistId, rowDocaCol, 
            rowActDistCol, idCol, energyCol, cornerDoca);

        // Create CalRecon object
        CalRecon *calRec = new CalRecon();
        calRec->initialize();

        UInt_t ienergy;
        for (ienergy = 0; ienergy < numEventEnergies; ienergy++ ) {
            CalEventEnergy *energy = new CalEventEnergy();
            energy->Fake(ienergy,randNum) ;
            calRec->addCalEventEnergy(energy);
        }

        UInt_t icluster;
        for (icluster = 0; icluster < numClusters; icluster++ ) {
            CalCluster *cluster = new CalCluster();
            cluster->Fake(icluster,randNum) ;
            calRec->addCalCluster(cluster);
        }

        UInt_t ixtal;
        for (ixtal = 0; ixtal < numXtals; ixtal++) {
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

        unsigned int imip;
        for (imip = 0; imip < numMipTracks; imip++) {
            CalMipTrack *mipTrack = new CalMipTrack(TVector3(1.0,1.5, 1.9),
               TVector3(2.0, 2.5, 2.9), 1.0, 2.0, 3, 4.0, 5.0, 6.0, 7.0, 8.0);
           calRec->addCalMipTrack(mipTrack);

        }

        // Create TkrRecon object
        TkrRecon *tkrRec = new TkrRecon();
        tkrRec->initialize();

        for (icluster = 0; icluster < numClusters; icluster++ ) 
        {
            UInt_t towerX = icluster % 4;
            UInt_t towerY = (icluster / 4) % 4;
            UInt_t tray   = icluster % 15;
            commonRootData::TkrId tkrId(towerX, towerY, tray, true, commonRootData::TkrId::eMeasureX);
            UInt_t strip0 = 4;
            UInt_t stripf = 10;
            TVector3 pos(f, 2.*f, 3.*f);
            UInt_t rawTot = 50;
            Double_t tot =  f;
            UInt_t flag = 1;
            UInt_t tower = 8;
            UInt_t nBad  = 0;
            TkrCluster *cluster = new TkrCluster(tkrId, strip0, stripf, pos, rawTot, tot, flag, nBad);
            tkrRec->addCluster(cluster);
        }

        TkrTrack *track;
        UInt_t itrack;
        for (itrack=0; itrack < numTracks; itrack++) 
        {
            track = new TkrTrack();
            UInt_t xgaps = itrack * itrack;
            UInt_t ygaps = itrack+itrack;
            UInt_t x1st = ievent+itrack;
            UInt_t y1st = ievent+ievent;
            track->setNumXGaps(xgaps);
            track->setNumYGaps(ygaps);
            track->setNumXFirstGaps(x1st);
            track->setNumYFirstGaps(y1st);
            Double_t chiSq = f;
            Double_t chiSqSmooth = f*randNum;
            Double_t rms = randNum*randNum;
            Double_t qual = f*f;
            Double_t e = f*2.0;
            Double_t ms = f*6.0;
            track->setChiSquareFilter(chiSq);
            track->setChiSquareSmooth(chiSqSmooth);
            track->setScatter(rms);
            track->setQuality(qual);
            track->setKalEnergy(e);
            track->setKalThetaMS(ms);
            tkrRec->addTrack(track);
        }

        UInt_t ivertex;
        TkrVertex *vertex;
        for (ivertex=0; ivertex < numVertices; ivertex++) {
            vertex = new TkrVertex();
            UInt_t towerX = ivertex % 4;
            UInt_t towerY = ivertex / 4;
            UInt_t tray   = ivertex % 15;
            commonRootData::TkrId tkrId(towerX, towerY, tray, true, commonRootData::TkrId::eMeasureX);
            vertex->setTkrID(tkrId);
            Double_t qual = f;
            vertex->setQuality(qual);
            Double_t energy = f*randNum;
            vertex->setEnergy(energy);

            TkrTrackParams vtxPar;
            vtxPar(1) = f;
            vtxPar(2) = f*f;
            vtxPar(3) = f*randNum;
            vtxPar(4) = randNum*randNum;

            vtxPar(1,1) = f;
            vtxPar(1,2) = f*2.;
            vtxPar(1,3) = f*3.;
            vtxPar(1,4) = f*4.;
            vtxPar(2,2) = f*5.;
            vtxPar(2,3) = f*6.;
            vtxPar(2,4) = f*7.;
            vtxPar(3,3) = f*8.;
            vtxPar(3,4) = f*9.;
            vtxPar(4,4) = f*10.;
            vertex->setParams(vtxPar);

            TVector3 pos(f, f*randNum, f*randNum*2.);
            vertex->setPosition(pos);
            TVector3 dir(randNum*2., randNum*3., randNum*4.);
            vertex->setDirection(dir);
            //vertex->addTrackId(ivertex+1);
            //vertex->addTrackId(ivertex+2);
            tkrRec->addVertex(vertex);
        }

        TkrDiagnostics *diag = new TkrDiagnostics(1, 2, 3, 4, 5, 6.0, 7, 8);
        tkrRec->addDiagnostics(diag);


        ev->initialize(ievent, runNum, tkrRec, calRec, acdRec);
        ev->initEventFlags(1);
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
    sc = write(fileName, numEvents);    sc = read(fileName, numEvents);

    if (sc == 0) {
        std::cout << "RECON ROOT file writing and reading succeeded!" << std::endl;
    } else {
        std::cout << "FAILED recon writing and reading" << std::endl;
    }

    return(sc);
}
