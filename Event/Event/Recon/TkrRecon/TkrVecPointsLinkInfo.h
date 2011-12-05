/** @file TkrVecPointsLinkInfo.h

* @class TkrVecPointsLinkInfo
*
* @brief This is a companion utility class for TkrVecPointsLink. It is meant to be 
*        be created when TkrVecPointsLink are created and to contained pieces of information
*        which are useful to various downstream processes - e.g. connecting links to form nodes
*
* last modified 11/29/11
*
* @authors Tracker Folks
*
* $Header$
*/

#ifndef __TkrVecPointsLinkInfo_H
#define __TkrVecPointsLinkInfo_H

#include "Event/Recon/TkrRecon/TkrVecPointsLink.h"
#include <vector>

// Declare Gaudi object interface ID
static const CLID& CLID_TkrVecPointsLinkInfo = InterfaceID("TkrVecPointsLinkInfo",  1, 0);

namespace Event {  // NameSpace

class TkrVecPointsLinkInfo: virtual public DataObject
{
public:
    // constructors
    TkrVecPointsLinkInfo() : m_tkrVecPointsLinkCol(0),
                             m_pointToLinksTab(0)
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrVecPointsLinkInfo::classID(); }
    static  const CLID& classID()         { return CLID_TkrVecPointsLinkInfo; }

    // destructor
    virtual ~TkrVecPointsLinkInfo()
    {
        // We own the actual table, but not the actual relations
        // delete the table but remember that the TDS will take care of the relations
        if (m_pointToLinksTab) delete m_pointToLinksTab;

        return;
    }
    /// @name Data set methods
    //@{
    /// Set the pointer to the collection of links
    void setTkrVecPointsLinkCol(Event::TkrVecPointsLinkCol* linksCol)     {m_tkrVecPointsLinkCol = linksCol;}
    /// Set the pointer to the table of link relations
    void setTkrVecPointToLinksTab(Event::TkrVecPointToLinksTab* linksTab) {m_pointToLinksTab     = linksTab;}
    //@}

    /// @name Data access methods
    //@{
    /// Recover a pointer to the collection of links
    const Event::TkrVecPointsLinkCol*   getTkrVecPointsLinkCol()   const {return m_tkrVecPointsLinkCol;}
    /// Recover a pointer to the link relations table
    const Event::TkrVecPointToLinksTab* getTkrVecPointToLinksTab() const {return m_pointToLinksTab;}
    //@}

private:

    // data members

    /// We will provide a pointer to the link collection but we don't own it and leave it
    /// to the TDS to manage
    Event::TkrVecPointsLinkCol*   m_tkrVecPointsLinkCol;

    /// Here we maintain a link to a pre-built table so that it does not need to be 
    /// reconstructed by downstream tools/algorithms that may need/use it.
    Event::TkrVecPointToLinksTab* m_pointToLinksTab;
};

}; // Namespace

#endif // __TkrVecPointsLinkInfo_H
