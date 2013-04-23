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
    
    enum StatusBits {ASSOCIATED   = 0x00000001,
                     FIRSTLINK    = 0x00000002,
                     SKIP1LAYER   = 0x00000010,
                     SKIP2LAYER   = 0x00000020,
                     SKIP3LAYER   = 0x00000040,
                     SKIPNLAYER   = 0x00000080,
                     SAMETOWER    = 0x00000100,
                     VERIFIED     = 0x00000800,
                     INTERTOWER   = 0x00010000,
                     WAFERGAP     = 0x00020000,
                     WAFERGAPPLUS = 0x00040000,
                     GAPANDCLUS   = 0x00080000,
                     TRUNCATED    = 0x00100000,
                     BADSTRIPS    = 0x00080000
                                         };

    enum skipsBits  {SKIPSLAYERS = SKIP1LAYER | SKIP2LAYER | SKIP3LAYER | SKIPNLAYER};

    // Constructors
    TkrVecPointsLink(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, double ang);

    virtual ~TkrVecPointsLink() 
    {
        return;
    }

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrVecPointsLink::classID(); }
    static  const CLID& classID()         { return CLID_TkrVecPointsLink; }

    void setUnAssociated()                {m_statusBits      &= ~ASSOCIATED;}
    void setAssociated()                  {m_statusBits      |=  ASSOCIATED;}
    void setNotFirstLink()                {m_statusBits      &= ~FIRSTLINK;}
    void setFirstLink()                   {m_statusBits      |=  FIRSTLINK;}
    void setSameTower()                   {m_statusBits      |=  SAMETOWER;}
    void setSkip1Layer()                  {m_statusBits      |=  SKIP1LAYER;}
    void setSkip2Layer()                  {m_statusBits      |=  SKIP2LAYER;}
    void setSkip3Layer()                  {m_statusBits      |=  SKIP3LAYER;}
    void setSkipNLayer()                  {m_statusBits      |=  SKIPNLAYER;}
    void setVerified()                    {m_statusBits      |=  VERIFIED;}
    void updateStatusBits(unsigned bits)  {m_statusBits      |=  bits;}
    void clearStatusBits(unsigned bits)   {m_statusBits      &= ~bits;}
    void setMaxScatAngle(double ang)      {m_maxScatAngle     =  ang;}
    void setAngleToNextLink(double ang)   {m_angleToNextLink  =  ang;}
    void setDistToNextLink(double dist)   {m_distToNextLink   =  dist;}
    void setRmsAngleSum(double ang)       {m_rmsAngleSum      =  ang;}
    void setRmsDistSum(double dist)       {m_rmsDistSum       =  dist;}
    void setNumAnglesInSum(int num)       {m_numAnglesInSum   =  num;}
    void setNumDistsInSum(int num)        {m_numDistsInSum    =  num;}
    void addAngleToRmsSum(double ang);
    void addDistToRmsSum(double dist);

    const bool         associated()          const {return (m_statusBits & ASSOCIATED)  == ASSOCIATED;}
    const bool         firstLink()           const {return (m_statusBits & FIRSTLINK)   == FIRSTLINK;}
    const bool         sameTower()           const {return (m_statusBits & SAMETOWER)   == SAMETOWER;}
    const bool         skipsLayers()         const {return (m_statusBits & SKIPSLAYERS) != 0;}
    const bool         skip1Layer()          const {return (m_statusBits & SKIP1LAYER)  == SKIP1LAYER;}
    const bool         skip2Layer()          const {return (m_statusBits & SKIP2LAYER)  == SKIP2LAYER;}
    const bool         skip3Layer()          const {return (m_statusBits & SKIP3LAYER)  == SKIP3LAYER;}
    const bool         skipNLayer()          const {return (m_statusBits & SKIPNLAYER)  == SKIPNLAYER;}
    const bool         verified()            const {return (m_statusBits & VERIFIED)    == VERIFIED;}
    const unsigned int getStatusBits()       const {return m_statusBits;}
    const Point&       getPosition()         const {return m_position;}
    const Point        getPosition(double z) const;
    const Point        getBotPosition()      const {return getPosition(m_secondVecPoint->getPosition().z());}
    const Vector&      getVector()           const {return m_vector;}
    const TkrVecPoint* getFirstVecPoint()    const {return m_firstVecPoint;}
    const TkrVecPoint* getSecondVecPoint()   const {return m_secondVecPoint;}
    const double       getMaxScatAngle()     const {return m_maxScatAngle;}
    const double       getAngleToNextLink()  const {return m_angleToNextLink;}
    const double       getDistToNextLink()   const {return m_distToNextLink;}
    const double       getRmsAngleSum()      const {return m_rmsAngleSum;}
    const double       getRmsDistSum()       const {return m_rmsDistSum;}
    const int          getNumAnglesInSum()   const {return m_numAnglesInSum;}
    const int          getNumDistsInSum()    const {return m_numDistsInSum;}

    // This used to match links
    const bool         matchFirst( const TkrVecPointsLink& linkToMatch) {return linkToMatch.getSecondVecPoint() == m_firstVecPoint;}
    const bool         matchSecond(const TkrVecPointsLink& linkToMatch) {return linkToMatch.getFirstVecPoint()  == m_secondVecPoint;}

    const double       angleToNextLink(const TkrVecPointsLink& linkToMatch);
    const double       distToNextLink( const TkrVecPointsLink& linkToMatch);

    
    const bool operator<(const TkrVecPointsLink* right) const {return m_position.z() < right->getPosition().z();}

	void setPosition(Point pos)           {m_position = pos; }
	void setVector(Vector dir)            {m_vector = dir; }

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
    double             m_distToNextLink;
    double             m_rmsAngleSum;
    double             m_rmsDistSum;
    int                m_numAnglesInSum;
    int                m_numDistsInSum;
};

inline TkrVecPointsLink::TkrVecPointsLink(const TkrVecPoint* firstPoint, const TkrVecPoint* secondPoint, double ang) : 
                                   m_firstVecPoint(firstPoint), 
                                   m_secondVecPoint(secondPoint), 
                                   m_statusBits(0), 
                                   m_maxScatAngle(ang),
                                   m_angleToNextLink(0.),
                                   m_distToNextLink(0.),
                                   m_rmsAngleSum(0.),
                                   m_rmsDistSum(0.),
                                   m_numAnglesInSum(0),
                                   m_numDistsInSum(0)
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
    
inline const Point TkrVecPointsLink::getPosition(double z) const
{
    // Set the position to the links position
    Point position = m_position;

    // First get the arclength of the step
    double arcLen = (z - position.z()) / m_vector.z();

    // Advance the link position to the new z
    position = position + arcLen * m_vector;

    return position;
}

inline const double TkrVecPointsLink::angleToNextLink(const TkrVecPointsLink& linkToMatch) 
{
    double cosAng = m_vector.dot(linkToMatch.getVector());

    m_angleToNextLink = acos(cosAng);

    return m_angleToNextLink;
}

inline const double TkrVecPointsLink::distToNextLink(const TkrVecPointsLink& linkToMatch) 
{
    // This assumes we are determining distance between point at end of THIS link to
    // the point at the start of the NEXT link
    Vector pointDiff = getPosition(m_secondVecPoint->getPosition().z()) - linkToMatch.getPosition();

//    m_distToNextLink = pointDiff.mag();

    double xDiffNorm = fabs(pointDiff.x()) / (0.03291 * m_secondVecPoint->getXCluster()->size());
    double yDiffNorm = fabs(pointDiff.y()) / (0.03291 * m_secondVecPoint->getYCluster()->size());
    
    m_distToNextLink = sqrt(xDiffNorm*xDiffNorm + yDiffNorm*yDiffNorm);

    return m_distToNextLink;
}

inline void TkrVecPointsLink::addAngleToRmsSum(double ang)   
{
    m_rmsAngleSum +=  ang * ang;
    m_numAnglesInSum++;

    return;
}

inline void TkrVecPointsLink::addDistToRmsSum(double dist)   
{
    m_rmsDistSum +=  dist * dist;
    m_numDistsInSum++;

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

