/** @file TkrVecPoint.h

* @class TkrVecPoint
*
* @brief This is the root version of the TkrVecPoint object. The purpose of the TkrVecPoint object is
*        to represent a 3D "point" in the Tracker definied by the crossing of two TkrClusters
*
* last modified 01/29/2013
*
* @authors Usher
*
* $Header: 
*/

#ifndef ROOT__TkrVecPoint_H
#define ROOT__TkrVecPoint_H

#include "TObject.h"
#include "commonRootData/idents/TkrId.h"
#include "TkrCluster.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

class TkrVecPoint: public TObject
{
public:
    // Status bit definitions    
    enum StatusBits {ASSOCIATED       = 0x00000001,   // Has been associated to a link/tree
                     LINKTOPHIT       = 0x00000010,   // Represents a top point in a link
                     LINKBOTHIT       = 0x00000020,   // Represents a bottom point in a link
                     PRNTLNKBOTHIT    = 0x00000040,   // 
                     ASSOCIATEDTONODE = 0x00000100,   // Point has been associated to a node
                     DONOTUSE         = 0x80000000};  // Not to be used in further Tree building

    // constructors - null case
    TkrVecPoint() {Clear();}

    // constructurs - have the information case
    TkrVecPoint(const TkrCluster* xClus, const TkrCluster* yClus):
         m_status(0), m_layer(xClus->getLayer()), m_pXCluster(xClus), m_pYCluster(yClus),
         m_position(xClus->getPosition().x(), yClus->getPosition().y(),
            0.5*(xClus->getPosition().z() + yClus->getPosition().z()))
    {}

    // destructor
         virtual ~TkrVecPoint() {Clear();}

    // Root mandated methods
    void Clear(Option_t *option="");
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const TkrVecPoint &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    // Initializer
    void initializeInfo(int layer, unsigned int status, const TkrCluster* xClus, const TkrCluster* yClus)
    {
        m_layer     = layer;
        m_status    = status;
        m_pXCluster = xClus;
        m_pYCluster = yClus;
        m_position  = TVector3(xClus->getPosition().x(), 
                               yClus->getPosition().y(),
                               0.5*(xClus->getPosition().z() + yClus->getPosition().z()));
    }

    /// Set status bits
    inline void setAssociated()         {m_status |=  ASSOCIATED;}
    inline void clearAssociated()       {m_status &= ~ASSOCIATED;}
    inline void setLinkTopHit()         {m_status |=  LINKTOPHIT;}
    inline void clearLinkTopHit()       {m_status &= ~LINKTOPHIT;}
    inline void setLinkBotHit()         {m_status |=  LINKBOTHIT;}
    inline void clearLinkBotHit()       {m_status &= ~LINKBOTHIT;}
    inline void setPrntLinkBotHit()     {m_status |=  PRNTLNKBOTHIT;}
    inline void clearPrntLinkBotHit()   {m_status &= ~PRNTLNKBOTHIT;}
    inline void setAssociatedToNode()   {m_status |=  ASSOCIATEDTONODE;}
    inline void clearAssociatedToNode() {m_status &= ~ASSOCIATEDTONODE;}
    inline void setDoNotUse()           {m_status |=  DONOTUSE;}
    inline void clearDoNotUse()         {m_status &= ~DONOTUSE;}

    /// @name access methods
    //@{
    /// Is this hit associated to a link?
    inline const bool isAssociated()               const {return (m_status & ASSOCIATED)       != 0;}
    /// Is this hit used as a top hit in a link?
    inline const bool isLinkTopHit()               const {return (m_status & LINKTOPHIT)       != 0;}
    /// Is this hit used as a bottom hit in a link?
    inline const bool isLinkBotHit()               const {return (m_status & LINKBOTHIT)       != 0;}
    /// Is the top hit of the link making this a bottom link also a bottom link?
    inline const bool isPrntLinkBotHit()           const {return (m_status & PRNTLNKBOTHIT)    != 0;}
    /// Is this hit associated to a node?
    inline const bool isAssociatedToNode()         const {return (m_status & ASSOCIATEDTONODE) != 0;}
    /// Is a usable point?
    inline const bool isUsablePoint()              const {return (m_status & DONOTUSE) == 0;}
    /// Pointer to the cluster in the x plane of this layer
    inline const TkrCluster* getXCluster()         const {return m_pXCluster;}
    /// Pointer to the cluster of the y plane of this layer
    inline const TkrCluster* getYCluster()         const {return m_pYCluster;}
    /// returns the ray between 2 TkrVecPoints, with corrections for slopes
//    Ray getRayTo(const TkrVecPoint* point)         const;
    /// at least one of the clusters in this point is flagged
    inline bool flagged()                          const { return m_pXCluster->hitFlagged() || m_pYCluster->hitFlagged(); }
    /// distance in layers between these two TkrVecPoints (always positive-definite)
    Int_t  getLayerSeparationFrom(const TkrVecPoint* point) const;
    /// the layer number of this TkrVecPoint (for Neural Net, could go away)
    inline Int_t  getLayer()                       const { return m_layer; }
    /// position of this TkrVecPoint, using info from x and y clusters
    inline const TVector3& getPosition()           const {return m_position;}
    /// Tower of this point... (x and y clusters are guaranteed to be in the same tower)
//    inline Int_t getTower()                        const { return m_pXCluster->tower(); }
    /// x/y distance to a reference point
    double getDistanceSquaredTo(TVector3 refPoint) const {
        TVector3 diff = refPoint - getPosition();
        return diff.x()*diff.x() + diff.y()*diff.y() + diff.z()*diff.z();    
    }

    // Allow access to full status word, mainly for converting to TDS
    UInt_t getStatusWord() {return m_status;}
    //@}

    /// @name other methods
    //@{
    /// equality operator - requires both clusters to match
    inline bool operator==(const TkrVecPoint& point) const;
    /// equality operator - requires only one cluster to match
    inline bool operator|=(const TkrVecPoint& point) const;
    /// inequality operator
    inline bool operator!=(const TkrVecPoint& point) const;
    //@}

    /// Override the "new" operator in order to use the ReconObjectManager pool management
    void* operator new(size_t size);

    void* operator new(size_t size, void* vp);

    /// Override the "delete" operator in order to use the ReconObjectManager pool management
    void  operator delete(void* p);

private:

    // data members
    /// Status word
    UInt_t            m_status;
    /// layer number
    Int_t             m_layer;
    /// pointer to x cluster
    const TkrCluster* m_pXCluster;
    /// pointer to y cluster
    const TkrCluster* m_pYCluster;
    /// Position from the above two clusters
    TVector3           m_position;

    ClassDef(TkrVecPoint,1)
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

//typedef std::list<TkrVecPoint*>         TkrVecPointList;
//typedef TkrVecPointList::const_iterator TkrVecPointListConItr;

// Typedefs for gaudi container for these objects
//typedef ObjectList<TkrVecPoint>        TkrVecPointCol;
//typedef TkrVecPointCol::iterator       TkrVecPointColPtr;
//typedef TkrVecPointCol::const_iterator TkrVecPointColConPtr;


#endif // __TkrVecPoint_H
