
#ifndef TkrDiagnostics_H
#define TkrDiagnostics_H

#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "idents/TkrId.h"

/** 
* @class TkrDiagnostics
*
* @brief Diagnostics output class
*
* @author Bill Atwood, Leon Rochester, Johann Cohen-Tanugi, Tracy Usher
*
* $Header$
*/

static const CLID& CLID_TkrDiagnostics = InterfaceID("TkrDiagnostics",  1, 1);

namespace Event { // Namespace

class TkrDiagnostics : virtual public DataObject
{
public:
    /// Default (null) constructor (just in case...)
    TkrDiagnostics() : m_numClusters(0), m_numVecPoints(0), m_numVecLinks(0),
                       m_nLinksNonZeroLayers(0), m_aveNumLinksLayer(0),
                       m_numLinkCombinations(0.), m_numTrackElements(0),
                       m_numTkrTracks(0) {};

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

    //! Destructor
    virtual ~TkrDiagnostics() {return;}

    //! Gaudi stuff: Retrieve pointer to class defininition structure
    virtual const CLID& clID()                 const   { return TkrDiagnostics::classID(); }
    static  const CLID& classID()                      { return CLID_TkrDiagnostics;       }

    inline const int    getNumClusters()         const {return m_numClusters;}
    inline const int    getNumVecPoints()        const {return m_numVecPoints;}
    inline const int    getNumVecLinks()         const {return m_numVecLinks;}
    inline const int    getnLinksNonZeroLayers() const {return m_nLinksNonZeroLayers;}
    inline const int    getAveNumLinksLayer()    const {return m_aveNumLinksLayer;}
    inline const double getNumLinkCombinations() const {return m_numLinkCombinations;}
    inline const int    getNumTrackElements()    const {return m_numTrackElements;}
    inline const int    getNumTkrTracks()        const {return m_numTkrTracks;}

    inline std::ostream& fillStream( std::ostream& s ) const;
    
private:
    /// Define here variables to keep diagnostic information for each event
    int     m_numClusters;          // Number of clusters this event
    int     m_numVecPoints;         // Resulting number of VecPoints this event
    int     m_numVecLinks;          // Number of links between VecPoints
    int     m_nLinksNonZeroLayers;  // Number of layers with links
    int     m_aveNumLinksLayer;     // Average number of links per layer
    double  m_numLinkCombinations;  // Keep track of expected number of combinations
    int     m_numTrackElements;     // Number of found TrackElements
    int     m_numTkrTracks;         // Number of tracks created 
};

inline std::ostream& Event::TkrDiagnostics::fillStream( std::ostream& s ) const 
{ 
  s << "Num Clusters              : " << getNumClusters()         << "\n"
    << "Num Vec Points            : " << getNumVecPoints()        << "\n"
    << "Num Vec Links             : " << getNumVecLinks()         << "\n"
    << "Num links non zero layers : " << getnLinksNonZeroLayers() << "\n"
    << "Ave Num links per layer   : " << getAveNumLinksLayer()    << "\n"
    << "Num Link Combinations     : " << getNumLinkCombinations() << "\n"
    << "Num  Track Elements       : " << getNumTrackElements()    << "\n"
    << "Num TkrTracks             : " << getNumTkrTracks()        << "\n";

  return s; 
}

}; //Namespace

#endif
