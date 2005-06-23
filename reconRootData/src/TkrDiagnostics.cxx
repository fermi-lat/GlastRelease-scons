#include "reconRootData/TkrDiagnostics.h"

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
