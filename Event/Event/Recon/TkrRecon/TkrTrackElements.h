/**
 * @class TkrTrackElements
 *
 * @brief Last level of the Vector Links pattern recognition. This 
 *        class contains the links associated together to form a track candidate
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef TkrTrackElements_h
#define TkrTrackElements_h

#include "Event/Recon/TkrRecon/TkrVecPointsLink.h"

#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/RelTable/RelTable.h"

// Declare Gaudi object interface ID
static const CLID& CLID_TkrTrackElements = InterfaceID("TkrTrackElements",  1, 0);

namespace Event {  // NameSpace

class TkrTrackElements: virtual public ContainedObject
{
public:
    enum StatusBits {UNUSED      = 0x0001,
                     USED        = 0x0002,
                     BESTTRK     = 0x0004,
                     NOTBEST     = 0x0008,
                     SHARE1STHIT = 0x0010,
                     SHARE2NDHIT = 0x0020,
                     SHARE3RDHIT = 0x0040};

    // Constructors
    TkrTrackElements(int numLinks, int numBiLayers, double rmsAngle, TkrVecPointsLink* firstLink) : 
                  m_statusBits(UNUSED), 
                  m_numLinks(numLinks), 
                  m_numBiLayers(numBiLayers),
                  m_rmsAngle(rmsAngle),
                  m_firstLink(firstLink)
                  {}

   virtual ~TkrTrackElements() 
   {
       return;
   }

    //! Retrieve pointer to class defininition structure
    virtual const CLID& clID()    const   { return TkrTrackElements::classID(); }
    static  const CLID& classID()         { return CLID_TkrTrackElements; }

    unsigned int      getStatusBits()  const {return m_statusBits;}
    int               getNumLinks()    const {return m_numLinks;}
    int               getNumBiLayers() const {return m_numBiLayers;}
    double            getRmsAngle()    const {return m_rmsAngle;}
    TkrVecPointsLink* getFirstLink()   const {return m_firstLink;}

    bool              isUnUsed()       const {return (m_statusBits & UNUSED) == UNUSED;}

    // Allow to reset number of bi layers encountered by track element
    void              setNumBiLayers(int numBiLayers) {m_numBiLayers = numBiLayers;}

    const bool operator<(const TkrTrackElements& right)  const;
    const bool operator==(const TkrTrackElements& right) const;

    void setStatusBits(unsigned int bits) {m_statusBits = bits;}
    void setUsed()                        {m_statusBits |= USED; m_statusBits &= ~UNUSED;}

private:
    unsigned int      m_statusBits;      // The status bits for this track element
    int               m_numLinks;        // The actual number of TkrVecPointsLinks in track element
    int               m_numBiLayers;     // The number of bilayers traversed by the track element
    double            m_rmsAngle;        // "rms" angle between links
    TkrVecPointsLink* m_firstLink;       // First TkrVecPointsLink on trackk element
};

inline const bool TkrTrackElements::operator<(const TkrTrackElements& right) const
{
//    if (m_numLinks > 16 && right.getNumLinks() > 16)     // 18 is max
//    {
//        // calculate weight factors
//        double leftWeight  = m_rmsAngle / (10. * m_numLinks);
//        double rightWeight = right.getRmsAngle() / (10. * right.getNumLinks());
//
//        if (leftWeight < rightWeight) return true;
//    }
//    else
//    {
//        int lenDiff = abs(m_numLinks - right.getNumLinks());
//
//        // Divide between close in length and not close
//        if (lenDiff > 0)
//        {
//            // Longer track wins
//            if (m_numLinks > right.getNumLinks()) return true;
//        }
//        else
//        {
//            // straighter wins
//            //if (m_numLinks == right.getNumLinks()) 
//            //{
//                if (m_rmsAngle < right.getRmsAngle()) return true;
//            //}
//        }
//    }
    // if big difference in length then take longer track
    if (m_numBiLayers - right.getNumBiLayers() > 3) return true;
    if (right.getNumBiLayers() - m_numBiLayers > 3) return false;

    double myRmsAngle = std::max(0.000001, m_rmsAngle);
    double myWeight   = m_numBiLayers * m_numBiLayers / myRmsAngle;

    double rtRmsAngle = std::max(0.000001, right.getRmsAngle());
    double rtWeight   = right.getNumBiLayers() * right.getNumBiLayers() / rtRmsAngle;

    double wghtRatio  = rtWeight / myWeight;

    // if the straightness is ~ the same, prefer the one starting higher up
    if (wghtRatio > 0.9 && wghtRatio < 1.1)
    {
        // Get start z coordinate for this track element
        double myPosZ = m_firstLink->getPosition().z();
        double rtPosZ = right.getFirstLink()->getPosition().z();

        if      (myPosZ > rtPosZ) return true;
        else if (rtPosZ > myPosZ) return false;
    }

    if (myWeight > rtWeight) return true;

    return false;
}

inline const bool TkrTrackElements::operator==(const TkrTrackElements& right) const
{
    // same length
    if (m_numBiLayers != right.getNumBiLayers()) return false;
    
    // and same angle
    if (m_rmsAngle != right.getRmsAngle()) return false;

    return true;
}

// These typedefs useful for associating tracks/hits
typedef std::vector<TkrTrackElements>                         TkrTrackElementsVec;
typedef std::vector<TkrTrackElements*>                        TkrTrackElementsPtrVec;

// Typedefs for gaudi container for these objects
typedef ObjectVector<TkrTrackElements>                        TkrTrackElementsCol;
typedef TkrTrackElementsCol::iterator                         TkrTrackElementsColPtr;
typedef TkrTrackElementsCol::const_iterator                   TkrTrackElementsColConPtr;

typedef RelTable<TkrTrackElements, TkrVecPoint>               TkrTrackElemToPointsTab;
typedef Relation<TkrTrackElements, TkrVecPoint>               TkrTrackElemToPointsRel;
typedef ObjectList< Relation<TkrTrackElements, TkrVecPoint> > TkrTrackElemToPointsTabList;

typedef RelTable<TkrTrackElements, TkrVecPointsLink>          TkrTrackElemToLinksTab;
typedef Relation<TkrTrackElements, TkrVecPointsLink>          TkrTrackElemToLinksRel;
typedef RelationList<TkrTrackElements, TkrVecPointsLink>      TkrTrackElemToLinksTabList;

typedef RelTable<TkrTrackElements, TkrCluster>                TkrTrackElemToClustersTab;
typedef Relation<TkrTrackElements, TkrCluster>                TkrTrackElemToClustersRel;
}; //Namespace


#endif

