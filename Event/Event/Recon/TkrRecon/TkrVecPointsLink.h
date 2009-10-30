/**
 * @class TkrVecPointLinks
 *
 * @brief This defines a class to associate TkrVecPoints into "links"
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef TkrVecPointsLink_h
#define TkrVecPointsLink_h

#include "Event/Recon/TkrRecon/TkrVecPoint.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/RelTable/RelTable.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrVecPointsLink = InterfaceID("TkrVecPointsLink",  1, 0);

namespace Event {  // NameSpace


class TkrVecPointsLink: virtual public ContainedObject
{
public:
    
    enum StatusBits {ASSOCIATED = 0x0001,
                     FIRSTLINK  = 0x0002,
                     SKIP1LAYER = 0x0010,
                     SKIP2LAYER = 0x0020,
                     SAMETOWER  = 0x0100};

    // Constructors
    TkrVecPointsLink(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, double ang);

    virtual ~TkrVecPointsLink() 
    {
        return;
    }

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrVecPointsLink::classID(); }
    static  const CLID& classID()         { return CLID_TkrVecPointsLink; }

    void setUnAssociated()              {m_statusBits      &= ~ASSOCIATED;}
    void setAssociated()                {m_statusBits      |=  ASSOCIATED;}
    void setNotFirstLink()              {m_statusBits      &= ~FIRSTLINK;}
    void setFirstLink()                 {m_statusBits      |=  FIRSTLINK;}
    void setSameTower()                 {m_statusBits      |=  SAMETOWER;}
    void setSkip1Layer()                {m_statusBits      |=  SKIP1LAYER;}
    void setSkip2Layer()                {m_statusBits      |=  SKIP2LAYER;}
    void setMaxScatAngle(double ang)    {m_maxScatAngle     =  ang;}
    void setAngleToNextLink(double ang) {m_angleToNextLink  =  ang;}
    void setRmsAngleSum(double ang)     {m_rmsAngleSum      =  ang;}
    void setNumAnglesInSum(int num)     {m_numAnglesInSum   =  num;}
    void addAngleToRmsSum(double ang);

    const bool         associated()         const {return (m_statusBits & ASSOCIATED) == ASSOCIATED;}
    const bool         firstLink()          const {return (m_statusBits & FIRSTLINK)  == FIRSTLINK;}
    const bool         sameTower()          const {return (m_statusBits & SAMETOWER)  == SAMETOWER;}
    const bool         skip1Layer()         const {return (m_statusBits & SKIP1LAYER) == SKIP1LAYER;}
    const bool         skip2Layer()         const {return (m_statusBits & SKIP2LAYER) == SKIP2LAYER;}
    const unsigned int getStatusBits()      const {return m_statusBits;}
    const Point&       getPosition()        const {return m_position;}
    const Vector&      getVector()          const {return m_vector;}
    const TkrVecPoint* getFirstVecPoint()   const {return m_firstVecPoint;}
    const TkrVecPoint* getSecondVecPoint()  const {return m_secondVecPoint;}
    const double       getMaxScatAngle()    const {return m_maxScatAngle;}
    const double       getAngleToNextLink() const {return m_angleToNextLink;}
    const double       getRmsAngleSum()     const {return m_rmsAngleSum;}
    const int          getNumAnglesInSum()  const {return m_numAnglesInSum;}

    // This used to match links
    const bool         matchFirst( const TkrVecPointsLink& linkToMatch) {return linkToMatch.getSecondVecPoint() == m_firstVecPoint;}
    const bool         matchSecond(const TkrVecPointsLink& linkToMatch) {return linkToMatch.getFirstVecPoint()  == m_secondVecPoint;}

    double             angleToNextLink(const TkrVecPointsLink& linkToMatch);

    
    const bool operator<(const TkrVecPointsLink* right) const {return m_position.z() < right->getPosition().z();}

private:
    const TkrVecPoint* m_firstVecPoint;
    const TkrVecPoint* m_secondVecPoint;
    
    unsigned int       m_statusBits;
    Point              m_position;
    Vector             m_vector;

    // Calculated expected maximum scattering angle
    double             m_maxScatAngle;

    // warning! a very volatile variables for specific use in the link attachment stage!
    double             m_angleToNextLink;
    double             m_rmsAngleSum;
    int                m_numAnglesInSum;
};

inline TkrVecPointsLink::TkrVecPointsLink(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, double ang) : 
                                   m_firstVecPoint(firstPoint), 
                                   m_secondVecPoint(secondPoint), 
                                   m_statusBits(0), 
                                   m_maxScatAngle(ang),
                                   m_angleToNextLink(0.),
                                   m_rmsAngleSum(0.),
                                   m_numAnglesInSum(0)
{
    // Are start and end points in the same tower?
    if (firstPoint->getXCluster()->tower() == secondPoint->getXCluster()->tower())
        m_statusBits = SAMETOWER;

    // Get the slopes in both X and Y to make the direction of this link
    double deltaX  = firstPoint->getXCluster()->position().x() - secondPoint->getXCluster()->position().x();
    double deltaXZ = firstPoint->getXCluster()->position().z() - secondPoint->getXCluster()->position().z();
    double slopeX  = deltaX / deltaXZ;

    double deltaY  = firstPoint->getYCluster()->position().y() - secondPoint->getYCluster()->position().y();
    double deltaYZ = firstPoint->getYCluster()->position().z() - secondPoint->getYCluster()->position().z();
    double slopeY  = deltaY / deltaYZ;

    Vector init_dir(-slopeX, -slopeY, -1.);
    m_vector       = init_dir.unit();

    // Now get the position which we'll set to be the bilayer z, between the two clusters
    double linkZ   = 0.5 * (firstPoint->getXCluster()->position().z() + firstPoint->getYCluster()->position().z());
    double linkX   = firstPoint->getXCluster()->position().x() + slopeX * (linkZ - firstPoint->getXCluster()->position().z());
    double linkY   = firstPoint->getYCluster()->position().y() + slopeY * (linkZ - firstPoint->getYCluster()->position().z());

    m_position     = Point(linkX, linkY, linkZ);
}

inline double TkrVecPointsLink::angleToNextLink(const TkrVecPointsLink& linkToMatch)
{
    double cosAng = m_vector.dot(linkToMatch.getVector());

    m_angleToNextLink = acos(cosAng);

    return m_angleToNextLink;
}

inline void TkrVecPointsLink::addAngleToRmsSum(double ang)   
{
    m_rmsAngleSum +=  ang * ang;
    m_numAnglesInSum++;

    return;
}

// These typedefs useful for point/link association
typedef std::vector<TkrVecPointsLink>               TkrVecPointsLinkVec;
typedef std::vector<TkrVecPointsLink*>              TkrVecPointsLinkPtrVec;

// Typedefs for gaudi container for these objects
typedef ObjectVector<TkrVecPointsLink>              TkrVecPointsLinkCol;
typedef TkrVecPointsLinkCol::iterator               TkrVecPointsLinkColPtr;
typedef TkrVecPointsLinkCol::const_iterator         TkrVecPointsLinkColConPtr;

typedef RelTable<TkrVecPoint, TkrVecPointsLink>     TkrVecPointToLinksTab;
typedef Relation<TkrVecPoint, TkrVecPointsLink>     TkrVecPointToLinksRel;
typedef RelationList<TkrVecPoint, TkrVecPointsLink> TkrVecPointToLinksTabList;

}; // Namespace

#endif

