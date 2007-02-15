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
class McTrajectory;
class McTrajectoryPoint;
class McPositionHit;
class McIntegratingHit;

// Forward declarations for recon information
class TkrCluster;
class TkrTrack;
class TkrTrackHit;

// typedefs for relating McParticles to associated McTrajectorys
typedef Event::RelTable<Event::McParticle, Event::McTrajectory>            McPartToTrajectoryTab;
typedef Event::Relation<Event::McParticle, Event::McTrajectory>            McPartToTrajectoryRel;
typedef ObjectList<McPartToTrajectoryRel>                                  McPartToTrajectoryTabList;
typedef std::vector<Event::McPartToTrajectoryRel*>                         McPartToTrajectoryVec;

// typedefs for relating McTrajectoryPoints to associated McPositionHits
typedef Event::RelTable<Event::McTrajectoryPoint, Event::McPositionHit>    McPointToPosHitTab;
typedef Event::Relation<Event::McTrajectoryPoint, Event::McPositionHit>    McPointToPosHitRel;
typedef ObjectList<McPointToPosHitRel>                                     McPointToPosHitTabList;
typedef std::vector<Event::McPointToPosHitRel*>                            McPointToPosHitVec;

// typedefs for relating McTrajectoryPoints to associated McIntegratingHits
typedef Event::RelTable<Event::McTrajectoryPoint, Event::McIntegratingHit> McPointToIntHitTab;
typedef Event::Relation<Event::McTrajectoryPoint, Event::McIntegratingHit> McPointToIntHitRel;
typedef ObjectList<McPointToIntHitRel>                                     McPointToIntHitTabList;
typedef std::vector<Event::McPointToIntHitRel*>                            McPointToIntHitVec;

// typedefs for relating McParticles to associated McPositionHits
typedef Event::RelTable<Event::McParticle, Event::McPositionHit>           McPartToPosHitTab;
typedef Event::Relation<Event::McParticle, Event::McPositionHit>           McPartToPosHitRel;
typedef ObjectList<McPartToPosHitRel>                                      McPartToPosHitTabList;
typedef std::vector<Event::McPartToPosHitRel*>                             McPartToPosHitVec;

// typedefs for relating McParticles to associated McIntegratingHits
typedef Event::RelTable<Event::McParticle, Event::McIntegratingHit>        McPartToIntHitTab;
typedef Event::Relation<Event::McParticle, Event::McIntegratingHit>        McPartToIntHitRel;
typedef ObjectList<McPartToIntHitRel>                                      McPartToIntHitTabList;
typedef std::vector<Event::McPartToIntHitRel*>                             McPartToIntHitVec;

// Below here relatest MC to TkrRecon objects
// typedefs for relating TkrClusters to McPositionHits 
typedef Event::RelTable<Event::TkrCluster, Event::McPositionHit>           ClusMcPosHitTab;
typedef Event::Relation<Event::TkrCluster, Event::McPositionHit>           ClusMcPosHitRel;
typedef ObjectList<ClusMcPosHitRel>                                        ClusMcPosHitTabList;
typedef std::vector<Event::ClusMcPosHitRel*>                               ClusMcPosHitVec;

// typedefs for relating McParticles to TkrClusters
typedef Event::RelTable<Event::McParticle, Event::TkrCluster>              McPartToClusTab;
typedef Event::Relation<Event::McParticle, Event::TkrCluster>              McPartToClusRel;
typedef ObjectList<McPartToClusRel>                                        McPartToClusTabList;
typedef std::vector<Event::McPartToClusRel*>                               McPartToClusVec;

// typedefs for relating McParticles to ClusMcPosHitRels
typedef Event::RelTable<Event::McParticle, ClusMcPosHitRel>                McPartToClusPosHitTab;
typedef Event::Relation<Event::McParticle, ClusMcPosHitRel>                McPartToClusPosHitRel;
typedef ObjectList<McPartToClusPosHitRel>                                  McPartToClusPosHitTabList;
typedef std::vector<Event::McPartToClusPosHitRel*>                         McPartToClusPosHitVec;

// typedefs for relating TkrTrackHits to McParticles
typedef Event::RelTable<Event::McParticle, Event::TkrTrackHit>             McPartToTkrTrackHitTab;
typedef Event::Relation<Event::McParticle, Event::TkrTrackHit>             McPartToTkrTrackHitRel;
typedef ObjectList<McPartToTkrTrackHitRel>                                 McPartToTkrTrackHitTabList;
typedef std::vector<McPartToTkrTrackHitRel*>                               McPartToTkrTrackHitVec;

// typedefs for relating TkrTracks to McParticles
typedef Event::RelTable<Event::McParticle, Event::TkrTrack>                McPartToTkrTrackTab;
typedef Event::Relation<Event::McParticle, Event::TkrTrack>                McPartToTkrTrackRel;
typedef ObjectList<Event::McPartToTkrTrackRel>                             McPartToTkrTrackTabList;
typedef std::vector<Event::McPartToTkrTrackRel*>                           McPartToTkrTrackVec;
};

#endif
