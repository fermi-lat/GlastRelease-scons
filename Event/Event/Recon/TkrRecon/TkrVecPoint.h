/** @file TkrVecPoint.h

* @class TkrVecPoint
*
* @brief Container class for the XY hit pairs which are produced by TkrVecPoints.
*
* last modified 11/01/2004
*
* @authors b. allgood, w. atwood and l. rochester
*
* $Header$
*/

#ifndef __TkrVecPoint_H
#define __TkrVecPoint_H

#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/Ray.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include <vector>

#include "GaudiKernel/ObjectList.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrVecPoint = InterfaceID("TkrVecPoint",  1, 0);

namespace Event {  // NameSpace

class TkrVecPoint: virtual public ContainedObject
{
public:
    // Status bit definitions    
    enum StatusBits {ASSOCIATED       = 0x00000001,
                     LINKTOPHIT       = 0x00000010,
                     LINKBOTHIT       = 0x00000020,
                     PRNTLNKBOTHIT    = 0x00000040,
                     ASSOCIATEDTONODE = 0x00000100,
                     DONOTUSE         = 0x80000000};

    // constructors
    TkrVecPoint() : m_layer(-1), m_pXCluster(0), m_pYCluster(0), m_position(0.,0.,0.) 
    {}

    TkrVecPoint(int layer, 
        const Event::TkrCluster* xClus, const Event::TkrCluster* yClus):
         m_status(0), m_layer(layer), m_pXCluster(xClus), m_pYCluster(yClus),
         m_position(xClus->position().x(), yClus->position().y(),
            0.5*(xClus->position().z() + yClus->position().z()))
    {}

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrVecPoint::classID(); }
    static  const CLID& classID()         { return CLID_TkrVecPoint; }

    // destructor
    virtual ~TkrVecPoint() 
    {
        return;
    }

    // Initializer
    void initialize(int layer, unsigned int status, const Event::TkrCluster* xClus, const Event::TkrCluster* yClus)
    {
        m_layer     = layer;
        m_status    = status;
        m_pXCluster = xClus;
        m_pYCluster = yClus;
    }

    /// Set status bits
    void setAssociated()         {m_status |=  ASSOCIATED;}
    void clearAssociated()       {m_status &= ~ASSOCIATED;}
    void setLinkTopHit()         {m_status |=  LINKTOPHIT;}
    void clearLinkTopHit()       {m_status &= ~LINKTOPHIT;}
    void setLinkBotHit()         {m_status |=  LINKBOTHIT;}
    void clearLinkBotHit()       {m_status &= ~LINKBOTHIT;}
    void setPrntLinkBotHit()     {m_status |=  PRNTLNKBOTHIT;}
    void clearPrntLinkBotHit()   {m_status &= ~PRNTLNKBOTHIT;}
    void setAssociatedToNode()   {m_status |=  ASSOCIATEDTONODE;}
    void clearAssociatedToNode() {m_status &= ~ASSOCIATEDTONODE;}
    void setDoNotUse()           {m_status |=  DONOTUSE;}
    void clearDoNotUse()         {m_status &= ~DONOTUSE;}

    /// @name access methods
    //@{
    /// Is this hit associated to a link?
    const bool isAssociated()       const {return (m_status & ASSOCIATED)       != 0;}
    /// Is this hit used as a top hit in a link?
    const bool isLinkTopHit()       const {return (m_status & LINKTOPHIT)       != 0;}
    /// Is this hit used as a bottom hit in a link?
    const bool isLinkBotHit()       const {return (m_status & LINKBOTHIT)       != 0;}
    /// Is the top hit of the link making this a bottom link also a bottom link?
    const bool isPrntLinkBotHit()   const {return (m_status & PRNTLNKBOTHIT)    != 0;}
    /// Is this hit associated to a node?
    const bool isAssociatedToNode() const {return (m_status & ASSOCIATEDTONODE) != 0;}
    /// Is a usable point?
    const bool isUsablePoint()      const {return (m_status & DONOTUSE) == 0;}
    /// Pointer to the cluster in the x plane of this layer
    const Event::TkrCluster*   getXCluster()   const {return m_pXCluster;}
    /// Pointer to the cluster of the y plane of this layer
    const Event::TkrCluster*   getYCluster()   const {return m_pYCluster;}
    /// returns the ray between 2 TkrVecPoints, with corrections for slopes
    Ray getRayTo(const TkrVecPoint* point) const;
    /// at least one of the clusters in this point is flagged
    bool flagged() const { return m_pXCluster->hitFlagged() || m_pYCluster->hitFlagged(); }
    /// distance in layers between these two TkrVecPoints (always positive-definite)
    int  getLayerSeparationFrom(const TkrVecPoint* point) const;
    /// the layer number of this TkrVecPoint (for Neural Net, could go away)
    int  getLayer() const { return m_layer; }
    /// position of this TkrVecPoint, using info from x and y clusters
    const Point& getPosition() const {return m_position;}
//        return Point(m_pXCluster->position().x(), m_pYCluster->position().y(),
//            0.5*(m_pXCluster->position().z() + m_pYCluster->position().z())); }
    /// Tower of this point... (x and y clusters are guaranteed to be in the same tower)
    int getTower() const { return m_pXCluster->tower(); }
    /// x/y distance to a reference point
    double getDistanceSquaredTo(Point refPoint) const {
        Vector diff = refPoint - getPosition();
        return diff.x()*diff.x() + diff.y()*diff.y() + diff.z()*diff.z();    
    }

    // Allow to retrieve the full status word (mainly for root storage)
    unsigned int getStatusWord() {return m_status;}

    //@}

    /// @name other methods
    //@{
    /// equality operator - requires both clusters to match
    bool operator==(const TkrVecPoint& point) const;
    /// equality operator - requires only one cluster to match
    bool operator|=(const TkrVecPoint& point) const;
    /// inequality operator
    bool operator!=(const TkrVecPoint& point) const;
    //@}

private:

    // data members
    /// Status word
    unsigned int             m_status;
    /// layer number
    int                      m_layer;
    /// pointer to x cluster
    const Event::TkrCluster* m_pXCluster;
    /// pointer to y cluster
    const Event::TkrCluster* m_pYCluster;
    /// Position from the above two clusters
    Point                    m_position;
};

inline int TkrVecPoint::getLayerSeparationFrom(const TkrVecPoint* point) const
{
    return abs(m_layer - point->m_layer);
}

inline bool TkrVecPoint::operator==(const TkrVecPoint& point) const
{
    // same point if both clusters matche (note &&)
    return( (m_pXCluster == point.m_pXCluster) &&
        (m_pYCluster == point.m_pYCluster) );
}

inline bool TkrVecPoint::operator|=(const TkrVecPoint& point) const
{
    // same point if either cluster matches (note ||)
    return( (m_pXCluster == point.m_pXCluster) ||
        (m_pYCluster == point.m_pYCluster) );
}

inline bool TkrVecPoint::operator!=(const TkrVecPoint& point) const
{
    // different point if different clusters
    return( (m_pXCluster != point.m_pXCluster) ||
            (m_pYCluster != point.m_pYCluster) );
}

typedef std::list<TkrVecPoint*>         TkrVecPointList;
typedef TkrVecPointList::const_iterator TkrVecPointListConItr;

// Typedefs for gaudi container for these objects
typedef ObjectList<TkrVecPoint>        TkrVecPointCol;
typedef TkrVecPointCol::iterator       TkrVecPointColPtr;
typedef TkrVecPointCol::const_iterator TkrVecPointColConPtr;

}; // Namespace

#endif // __TkrVecPoint_H
