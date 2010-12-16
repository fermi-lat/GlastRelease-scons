/** @file TkrVecPointInfo.h

* @class TkrVecPointInfo
*
* @brief This is a companion utility class for TkrVecPoints. It is meant to be 
*        be created when TkrVecPoints are created and to contained pieces of information
*        which are useful to various downstream processes - e.g. building links
*
* last modified 12/08/2010
*
* @authors Tracker Folks
*
* $Header$
*/

#ifndef __TkrVecPointInfo_H
#define __TkrVecPointInfo_H

#include "Event/Recon/TkrRecon/TkrVecPoint.h"
#include <vector>

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrVecPointInfo = InterfaceID("TkrVecPointInfo",  1, 0);

namespace Event {  // NameSpace

/// typedefs to define a mapping between bilayers and the start/end points in the TkrVecPointCol
typedef std::pair<Event::TkrVecPointColPtr, Event::TkrVecPointColPtr> TkrVecPointItrPair;
typedef std::map<int, TkrVecPointItrPair >                            TkrLyrToVecPointItrMap;
typedef TkrLyrToVecPointItrMap::iterator                              TkrLyrToVecPointItrMapItr;
typedef TkrLyrToVecPointItrMap::const_iterator                        TkrLyrToVecPointItrMapConsItr;

class TkrVecPointInfo: virtual public DataObject
{
public:
    // constructors
    TkrVecPointInfo() : m_maxNumSkippedLayers(0),
                        m_numVecPoints(0),
                        m_numBiLayersWVecPoints(0),
                        m_maxNumLinkCombinations(0),
                        m_tkrLyrToVecPointItrMap(0)
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrVecPointInfo::classID(); }
    static  const CLID& classID()         { return CLID_TkrVecPointInfo; }

    // destructor
    virtual ~TkrVecPointInfo()
    {
        // We own this map so we need to take care of it
        if (m_tkrLyrToVecPointItrMap) delete m_tkrLyrToVecPointItrMap;

        return;
    }
    /// @name Data set methods
    //@{
    /// Set the maximum number of skipped layers
    void setMaxNumSkippedLayers(int maxNumSkippedLayers)          {m_maxNumSkippedLayers    = maxNumSkippedLayers;}
    void setNumTkrVecPoints(int numVecPoints)                     {m_numVecPoints           = numVecPoints;}
    void setNumBiLayersWVecPoints(int numBiLayersWVecPoints)      {m_numBiLayersWVecPoints  = numBiLayersWVecPoints;}
    void setMaxNumLinkCombinations(double maxNumLinkCombinations) {m_maxNumLinkCombinations = maxNumLinkCombinations;}
    void setLyrToVecPointItrMap(TkrLyrToVecPointItrMap* map)      {m_tkrLyrToVecPointItrMap = map;}
    //@}

    /// @name Data access methods
    //@{
    /// How many layers are we allowed to skip when building links? 
    int                     getMaxNumSkippedLayers()    {return m_maxNumSkippedLayers;}
    /// Total number of TkrVecPoints in the collection? 
    int                     getNumTkrVecPoints()        {return m_numVecPoints;}
    /// Number of bilayers with TkrVecPoints
    int                     getNumBiLayersWVecPoints()  {return m_numBiLayersWVecPoints;}
    /// Get the maximum number of link combinations which are possible
    double                  getMaxNumLinkCombinations() {return m_maxNumLinkCombinations;}

    /// Access the mapping between bilayers and start/end iterators
    TkrLyrToVecPointItrMap* getLyrToVecPointItrMap()    {return m_tkrLyrToVecPointItrMap;}
    //@}

private:

    // data members
    /// How many layers are we allowed to skip?
    int                     m_maxNumSkippedLayers;

    /// Also keep track of the total number of TkrVecPoints
    int                     m_numVecPoints;

    /// Keep count of the number of bilayers with TkrVecPoints
    int                     m_numBiLayersWVecPoints;
 
    /// Finally, keep track of max possible link combinations
    double                  m_maxNumLinkCombinations;

    /// Pointer to our map
    TkrLyrToVecPointItrMap* m_tkrLyrToVecPointItrMap;
};

}; // Namespace

#endif // __TkrVecPointInfo_H
