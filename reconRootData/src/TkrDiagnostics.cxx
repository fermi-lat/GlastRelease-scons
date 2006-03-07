#include "reconRootData/TkrDiagnostics.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(TkrDiagnostics)

TkrDiagnostics::TkrDiagnostics() : m_numClusters(0), 
                                   m_numVecPoints(0), 
                                   m_numVecLinks(0),
                                   m_nLinksNonZeroLayers(0), 
                                   m_aveNumLinksLayer(0),
                                   m_numLinkCombinations(0.), 
                                   m_numTrackElements(0),
                                   m_numTkrTracks(0) 
{
    return;
}

void TkrDiagnostics::initializeInfo(int    numClusters, 
                                    int    numVecPoints, 
                                    int    numVecLinks, 
                                    int    nLinksNonZeroLayers,
                                    int    aveNumLinksLayer, 
                                    double numLinkCombinations, 
                                    int    numTrackElements,
                                    int    numTkrTracks)
{
    m_numClusters         = numClusters;
    m_numVecPoints        = numVecPoints;
    m_numVecLinks         = numVecLinks;
    m_nLinksNonZeroLayers = nLinksNonZeroLayers;
    m_aveNumLinksLayer    = aveNumLinksLayer;
    m_numLinkCombinations = numLinkCombinations;
    m_numTrackElements    = numTrackElements;
    m_numTkrTracks        = numTkrTracks;
}


//======================================================
// For Unit Tests
//======================================================

void TkrDiagnostics::Fake( Int_t /*ievent*/, UInt_t /*rank*/, Float_t /*randNum*/ ) {

    initializeInfo(1, 2, 3, 4, 5, 6.0, 7, 8) ;

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),a.get ## att(),#att)

Bool_t TkrDiagnostics::CompareInRange( const TkrDiagnostics & a, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(NumClusters) && result ;
    result = COMPARE_IN_RANGE(NumVecPoints) && result ;
    result = COMPARE_IN_RANGE(NumVecLinks) && result ;
    result = COMPARE_IN_RANGE(nLinksNonZeroLayers) && result ;
    result = COMPARE_IN_RANGE(AveNumLinksLayer) && result ;
    result = COMPARE_IN_RANGE(NumLinkCombinations) && result ;
    result = COMPARE_IN_RANGE(NumTrackElements) && result ;
    result = COMPARE_IN_RANGE(NumTkrTracks) && result ;
    
    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

