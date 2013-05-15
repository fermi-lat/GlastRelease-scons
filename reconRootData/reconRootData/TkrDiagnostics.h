#ifndef ROOT_TkrDiagnostics_H
#define ROOT_TkrDiagnostics_H

#include "TObject.h"
#include <string>

/** 
 * @class TkrDiagnostics
 * @brief Root object containing Tracker Recon Diagnostics
 * 
 * @author Tracy Usher
 *
 */

class TkrDiagnostics : public TObject
{
public:

    /// Default (null) constructor (just in case...)
    /// Default (null) constructor (just in case...)
    TkrDiagnostics();

    /// Construct all but the track parameters, they must be set during recon stage
    TkrDiagnostics(int numClusters, int numVecPoints, int numVecLinks, int nLinksNonZeroLayers,
                   int aveNumLinksLayer, double numLinkCombinations, int numTrackElements,
                   int numTkrTracks) :
                       m_numClusters(numClusters), 
                       m_numVecPoints(numVecPoints), 
                       m_numVecLinks(numVecLinks), 
                       m_nLinksNonZeroLayers(nLinksNonZeroLayers), 
                       m_aveNumLinksLayer(aveNumLinksLayer), 
                       m_numLinkCombinations(numLinkCombinations), 
                       m_numTrackElements(numTrackElements),
                       m_numTkrTracks(numTkrTracks) {};


    virtual ~TkrDiagnostics() {}

    void Clear(Option_t * /* option= "" */)
    {
        initializeInfo(0,0,0,0,0,0.,0,0);
    }
    
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrDiagnostics &, const std::string & name = "" ) const ; // for tests

    inline void Print(Option_t *option="") const {TObject::Print(option);};

    //! Two initialize methods needed for this class
    void initializeInfo(int    numClusters, 
                        int    numVecPoints, 
                        int    numVecLinks, 
                        int    nLinksNonZeroLayers,
                        int    aveNumLinksLayer, 
                        double numLinkCombinations, 
                        int    numTrackElements,
                        int    numTkrTracks);


    inline const UInt_t   getNumClusters()         const {return m_numClusters;}
    inline const UInt_t   getNumVecPoints()        const {return m_numVecPoints;}
    inline const UInt_t   getNumVecLinks()         const {return m_numVecLinks;}
    inline const UInt_t   getnLinksNonZeroLayers() const {return m_nLinksNonZeroLayers;}
    inline const UInt_t   getAveNumLinksLayer()    const {return m_aveNumLinksLayer;}
    inline const Double_t getNumLinkCombinations() const {return m_numLinkCombinations;}
    inline const UInt_t   getNumTrackElements()    const {return m_numTrackElements;}
    inline const UInt_t   getNumTkrTracks()        const {return m_numTkrTracks;}
    
private:
    /// Define here variables to keep diagnostic information for each event
    UInt_t     m_numClusters;          // Number of clusters this event
    UInt_t     m_numVecPoints;         // Resulting number of VecPoints this event
    UInt_t     m_numVecLinks;          // Number of links between VecPoints
    UInt_t     m_nLinksNonZeroLayers;  // Number of layers with links
    UInt_t     m_aveNumLinksLayer;     // Average number of links per layer
    Double32_t m_numLinkCombinations;  // Keep track of expected number of combinations
    UInt_t     m_numTrackElements;     // Number of found TrackElements
    UInt_t     m_numTkrTracks;         // Number of tracks created 

    ClassDef(TkrDiagnostics,1)
};


#endif
