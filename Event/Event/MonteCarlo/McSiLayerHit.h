/**
 * @class McSiLayerHit
 *
 * @brief Represents a Monte Carlo hit in a tracker layer
 *
 * @author The Tracking Software Group
 *
 * $Header$
 */
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRefVector.h"
#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/MonteCarlo/McParticle.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "idents/VolumeIdentifier.h"
#include "Event/RelTable/RelTable.h"

#ifndef McSiLayerHit_h
#define McSiLayerHit_h

namespace Event {
class McSiLayerHit : virtual public ContainedObject
{
public:
    //! Define some bits to help classify the event
    enum StatusBits{  
        CLUSTERHIT = 1 ,    //! This McSiLayerHit is associated with a valid TkrCluster
        SHAREDCLUS = 1<<1,  //! The TkrCluster is shared with another McSiLayerHit
        OWNCLUSTER = 1<<2,  //! The TkrCluster is "owned" by this McSiLayerHit
    };

    /// Standard Gaudi Tool interface constructor
    McSiLayerHit(const Event::McParticle* particle);
   ~McSiLayerHit();

    void addMcPositionHit(const Event::McPositionHit* posHit);
    void setTkrCluster(const Event::TkrCluster* cluster);
    void setStatusBit(StatusBits bitToSet) {m_statusBits |= bitToSet;}
    void clearStatusBit(StatusBits bitToClear) {m_statusBits &= ~bitToClear;}
    
    const Event::McParticle*                    getMcParticle()        const {return  m_McParticle;}
    const Event::TkrCluster*                    getTkrCluster()        const {return  m_cluster;}
    const SmartRefVector<Event::McPositionHit>* getMcPositionHitsVec() const {return &m_PositionHitsVec;}
    const idents::VolumeIdentifier              getVolumeIdent()       const {return  m_volIdent;}
    const unsigned long                         getStatusBits()        const {return  m_statusBits;}
    const HepPoint3D                            getHitPosition()       const;

private:
    const Event::McParticle*             m_McParticle;
    SmartRefVector<Event::McPositionHit> m_PositionHitsVec;
    const Event::TkrCluster*             m_cluster;
    idents::VolumeIdentifier             m_volIdent;
    unsigned long                        m_statusBits;
};

// typedefs for the cluster to McPositionHits 
typedef Event::RelTable<Event::TkrCluster, Event::McPositionHit>   ClusMcPosHitTab;
typedef Event::Relation<Event::TkrCluster, Event::McPositionHit>   ClusMcPosHitRel;
typedef ObjectList<ClusMcPosHitRel>                                ClusMcPosHitTabList;

// typedefs for the mc particles to McSiLayerHits
typedef Event::RelTable<Event::McParticle, Event::McSiLayerHit>    McPartToHitTab;
typedef Event::Relation<Event::McParticle, Event::McSiLayerHit>    McPartToHitRel;
typedef ObjectList<McPartToHitRel>                                 McPartToHitTabList;
typedef std::vector<Event::McPartToHitRel*>                        McPartToHitVec;

// typedes for MC "tracks" (McParticles) and the hits associated with them
//typedef std::vector<const Event::McParticle*> McPartTracks;
//typedef SmartRefVector<Event::McParticle>   McPartTracks;
typedef std::vector<Event::McPartToHitRel*> McPartTrack;

// typedefs for relating TkrClusters to McSiLayerHits
typedef Event::RelTable<Event::TkrCluster, Event::McSiLayerHit>    ClusToLyrHitTab;
typedef Event::Relation<Event::TkrCluster, Event::McSiLayerHit>    ClusToLyrHitRel;
typedef ObjectList<ClusToLyrHitRel>                                ClusToLyrHitTabList;
typedef std::vector<Event::ClusToLyrHitRel*>                       ClusToLyrHitVec;

// typedefs for relating McSiLayerHits to McPositionHits
typedef Event::RelTable<Event::McSiLayerHit, Event::McPositionHit> McLyrToHitTab;
typedef Event::Relation<Event::McSiLayerHit, Event::McPositionHit> McLyrToHitRel;
typedef ObjectList<McLyrToHitRel>                                  McLyrToHitTabList;
typedef std::vector<Event::McLyrToHitRel*>                         McLyrToHitVec;


typedef ObjectVector<McSiLayerHit> McSiLayerHitCol;
};

#endif