/**
 * @class McRelTableDefs
 *
 * @brief This header file serves to define the various relational tables used with 
 *        the Monte Carlo information
 *
 * @author Tracy Usher
 *
 * $Header$
 */

#ifndef McRelTableDefs_h
#define McRelTableDefs_h

// Necessary Gaudi stuff
#include "GaudiKernel/SmartRefVector.h"
#include "Event/RelTable/RelTable.h"

// Start the defintions
namespace Event {

// Monte Carlo objects with the "truth"
class McParticle;
class McPositionHit;
class McSiLayerHit;
class McCandTrack;

// Forward declarations for recon information
class TkrCluster;
class TkrTrack;
class TkrTrackHit;
class TkrPatCandHit;
class TkrPatCand;

// typedefs for relating McParticles to associated McPositionHits
typedef Event::RelTable<Event::McParticle, Event::McPositionHit>   McPartToPosHitTab;
typedef Event::Relation<Event::McParticle, Event::McPositionHit>   McPartToPosHitRel;
typedef ObjectList<McPartToPosHitRel>                              McPartToPosHitTabList;
typedef std::vector<Event::McPartToPosHitRel*>                     McPartToPosHitVec;

// typedefs for relating TkrClusters to McPositionHits 
typedef Event::RelTable<Event::TkrCluster, Event::McPositionHit>   ClusMcPosHitTab;
typedef Event::Relation<Event::TkrCluster, Event::McPositionHit>   ClusMcPosHitRel;
typedef ObjectList<ClusMcPosHitRel>                                ClusMcPosHitTabList;
typedef std::vector<Event::ClusMcPosHitRel*>                       ClusMcPosHitVec;

// typedefs for relating McParticles to TkrClusters
typedef Event::RelTable<Event::McParticle, Event::TkrCluster>      McPartToClusTab;
typedef Event::Relation<Event::McParticle, Event::TkrCluster>      McPartToClusRel;
typedef ObjectList<McPartToClusRel>                                McPartToClusTabList;
typedef std::vector<Event::McPartToClusRel*>                       McPartToClusVec;

// typedefs for relating McParticles to ClusMcPosHitRels
typedef Event::RelTable<Event::McParticle, ClusMcPosHitRel>        McPartToClusPosHitTab;
typedef Event::Relation<Event::McParticle, ClusMcPosHitRel>        McPartToClusPosHitRel;
typedef ObjectList<McPartToClusPosHitRel>                          McPartToClusPosHitTabList;
typedef std::vector<Event::McPartToClusPosHitRel*>                 McPartToClusPosHitVec;

// typedefs for relating TkrPatCandHits to McParticles
typedef Event::RelTable<Event::McParticle, Event::TkrPatCandHit>   McPartToTkrCandHitTab;
typedef Event::Relation<Event::McParticle, Event::TkrPatCandHit>   McPartToTkrCandHitRel;
typedef ObjectList<McPartToTkrCandHitRel>                          McPartToTkrCandHitTabList;
typedef std::vector<McPartToTkrCandHitRel*>                        McPartToTkrCandHitVec;

// typedefs for relating TkrPatCands to McParticles
typedef Event::RelTable<Event::McParticle, Event::TkrPatCand>      McPartToTkrPatCandTab;
typedef Event::Relation<Event::McParticle, Event::TkrPatCand>      McPartToTkrPatCandRel;
typedef ObjectList<Event::McPartToTkrPatCandRel>                   McPartToTkrPatCandTabList;
typedef std::vector<Event::McPartToTkrPatCandRel*>                 McPartToTkrPatCandVec;

// typedefs for relating TkrTrackHits to McParticles
typedef Event::RelTable<Event::McParticle, Event::TkrTrackHit>     McPartToTkrTrackHitTab;
typedef Event::Relation<Event::McParticle, Event::TkrTrackHit>     McPartToTkrTrackHitRel;
typedef ObjectList<McPartToTkrTrackHitRel>                         McPartToTkrTrackHitTabList;
typedef std::vector<McPartToTkrTrackHitRel*>                       McPartToTkrTrackHitVec;

// typedefs for relating TkrTracks to McParticles
typedef Event::RelTable<Event::McParticle, Event::TkrTrack>        McPartToTkrTrackTab;
typedef Event::Relation<Event::McParticle, Event::TkrTrack>        McPartToTkrTrackRel;
typedef ObjectList<Event::McPartToTkrTrackRel>                     McPartToTkrTrackTabList;
typedef std::vector<Event::McPartToTkrTrackRel*>                   McPartToTkrTrackVec;
};

#endif
