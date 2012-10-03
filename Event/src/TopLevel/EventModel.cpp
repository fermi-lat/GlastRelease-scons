// File and Version Information:
// $Header$

#define _Event_EventModel_CPP_

// Define this in order to export symbols
#define EVT_DLL_EXPORTS

#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/ClassID.h"

std::string EventModel::EventHeader = "/Event";

/** 
 *  @class MC
 *
*/
std::string EventModel::MC::Event                     = EventModel::EventHeader + "/MC";
std::string EventModel::MC::McParticleCol             = EventModel::MC::Event  + "/McParticleCol";

std::string EventModel::MC::McPositionHitCol          = EventModel::MC::Event  + "/PositionHitsCol";
std::string EventModel::MC::McIntegratingHitCol       = EventModel::MC::Event  + "/IntegratingHitsCol";
std::string EventModel::MC::McTrajectoryCol           = EventModel::MC::Event  + "/TrajectoryCol";
std::string EventModel::MC::McTkrStripCol             = EventModel::MC::Event  + "/StripCol";
std::string EventModel::MC::D2EntryCol                = EventModel::MC::Event  + "/D2EntryCol";
std::string EventModel::MC::ExposureCol               = EventModel::MC::Event  + "/ExposureCol";

std::string EventModel::MC::McEventStructure          = EventModel::MC::Event  + "/McEventStructure";
std::string EventModel::MC::McPartToPosHitTab         = EventModel::MC::Event  + "/McPartToPosHitTab";
std::string EventModel::MC::McPartToClusTab           = EventModel::MC::Event  + "/McPartToClusTab";
std::string EventModel::MC::McPartToClusHitTab        = EventModel::MC::Event  + "/McPartToClusHitTab";

std::string EventModel::MC::McPartToTkrCandHitTab     = EventModel::MC::Event  + "/McPartToTkrCandHitTab";
std::string EventModel::MC::McPartToTkrPatCandTab     = EventModel::MC::Event  + "/McPartToTkrPatCandTab";
std::string EventModel::MC::McPartToTkrTrackHitTab    = EventModel::MC::Event  + "/McPartToTkrTrackHitTab";
std::string EventModel::MC::McPartToTkrTrackTab       = EventModel::MC::Event  + "/McPartToTkrTrackTab";

std::string EventModel::MC::McAcdTkrPointCol          = EventModel::MC::Event  + "/AcdTkrPointCol";
std::string EventModel::MC::McAcdTkrHitPocaCol        = EventModel::MC::Event  + "/AcdTkrHitPocaCol";
std::string EventModel::MC::McAcdTkrAssocCol          = EventModel::MC::Event  + "/AcdTkrAssocCol";


// Digi event
std::string EventModel::Digi::Event                   = EventModel::EventHeader + "/Digi";
std::string EventModel::Digi::AcdDigiCol              = EventModel::Digi::Event + "/AcdDigiCol";
std::string EventModel::Digi::TkrDigiCol              = EventModel::Digi::Event + "/TkrDigiCol";
std::string EventModel::Digi::CalDigiCol              = EventModel::Digi::Event + "/CalDigiCol";
std::string EventModel::Digi::CalDigiHitTab           = EventModel::Digi::Event + "/CalDigiHitTab";
std::string EventModel::Digi::TkrDigiHitTab           = EventModel::Digi::Event + "/TkrDigiHitTab";
std::string EventModel::Digi::TkrClusterHitTab        = EventModel::Digi::Event + "/TkrClusterHitTab";

// Digi Overlay 
std::string EventModel::Overlay::Event                = EventModel::EventHeader + "/Overlay";
std::string EventModel::Overlay::EventHeader          = EventModel::Overlay::Event + "/EventHeader";
std::string EventModel::Overlay::TriRowBits           = EventModel::Overlay::Event + "/TriRowBits";
std::string EventModel::Overlay::Time                 = EventModel::Overlay::Event + "/Time";
std::string EventModel::Overlay::EventSummary         = EventModel::Overlay::Event + "/EventSummary";
std::string EventModel::Overlay::Gem                  = EventModel::Overlay::Event + "/Gem";
std::string EventModel::Overlay::Error                = EventModel::Overlay::Event + "/Error";
std::string EventModel::Overlay::Diagnostic           = EventModel::Overlay::Event + "/Diagnostic";
std::string EventModel::Overlay::ObfFilterStatus      = EventModel::Overlay::Event + "/ObfFilterStatus";
std::string EventModel::Overlay::ObfFilterTrack       = EventModel::Overlay::Event + "/ObfFilterTrack";
std::string EventModel::Overlay::MetaEvent            = EventModel::Overlay::Event + "/MetaEvent";
std::string EventModel::Overlay::Ccsds                = EventModel::Overlay::Event + "/Ccsds";
std::string EventModel::Overlay::AncillaryEvent       = EventModel::Overlay::Event + "/AncillaryEvent";
std::string EventModel::Overlay::AncillaryEventDigi   = EventModel::Overlay::Event + "/AncillaryEventDigi";
std::string EventModel::Overlay::AcdDigiCol           = EventModel::Overlay::Event + "/AcdDigiCol";
std::string EventModel::Overlay::TkrDigiCol           = EventModel::Overlay::Event + "/TkrDigiCol";
std::string EventModel::Overlay::CalDigiCol           = EventModel::Overlay::Event + "/CalDigiCol";
std::string EventModel::Overlay::CalDigiHitTab        = EventModel::Overlay::Event + "/CalDigiHitTab";
std::string EventModel::Overlay::TkrDigiHitTab        = EventModel::Overlay::Event + "/TkrDigiHitTab";
std::string EventModel::Overlay::TkrClusterHitTab     = EventModel::Overlay::Event + "/TkrClusterHitTab";
      
// reconstructed data (Overall)
std::string EventModel::Recon::Event                         = EventModel::EventHeader + "/Recon";
std::string EventModel::Recon::TreeClusterRelationCol        = EventModel::Recon::Event + "/TreeClusterRelationCol";
std::string EventModel::Recon::TreeToRelationMap             = EventModel::Recon::Event + "/TreeToRelationMap";
std::string EventModel::Recon::ClusterToRelationMap          = EventModel::Recon::Event + "/ClusterToRelationMap";
      
// reconstructed data (Tracker)
std::string EventModel::TkrRecon::Event                      = EventModel::EventHeader + "/TkrRecon";
std::string EventModel::TkrRecon::TkrClusterCol              = EventModel::TkrRecon::Event + "/TkrClusterCol";
std::string EventModel::TkrRecon::TkrIdClusterMMap           = EventModel::TkrRecon::Event + "/TkrIdClusterMMap";
std::string EventModel::TkrRecon::TkrIdClusterMap            = EventModel::TkrRecon::Event + "/TkrIdClusterMap";
std::string EventModel::TkrRecon::TkrTreeCol                 = EventModel::TkrRecon::Event + "/TkrTreeCol";
std::string EventModel::TkrRecon::TkrTrackCol                = EventModel::TkrRecon::Event + "/TkrTrackCol";
std::string EventModel::TkrRecon::TkrCRTrackCol              = EventModel::TkrRecon::Event + "/TkrCRTrackCol";
std::string EventModel::TkrRecon::TkrTrackMap                = EventModel::TkrRecon::Event + "/TkrTrackMap";
std::string EventModel::TkrRecon::TkrTrackHitCol             = EventModel::TkrRecon::Event + "/TkrTrackHitCol";
std::string EventModel::TkrRecon::TkrVertexCol               = EventModel::TkrRecon::Event + "/TkrVertexCol";
std::string EventModel::TkrRecon::TkrDiagnostics             = EventModel::TkrRecon::Event + "/TkrDiagnostics";
std::string EventModel::TkrRecon::TkrEventParams             = EventModel::TkrRecon::Event + "/TkrEventParams";
std::string EventModel::TkrRecon::TkrFilterParamsToBoxTab    = EventModel::TkrRecon::Event + "/TkrEventParamsToBoxTab";
std::string EventModel::TkrRecon::TkrFilterParamsToLinksTab  = EventModel::TkrRecon::Event + "/TkrEventParamsToLinksTab";
std::string EventModel::TkrRecon::TkrFilterParamsToPointsTab = EventModel::TkrRecon::Event + "/TkrEventParamsToPointsTab";
std::string EventModel::TkrRecon::TkrFilterParamsCol         = EventModel::TkrRecon::Event + "/TkrFilterParamsCol";
std::string EventModel::TkrRecon::TkrTruncatedPlane          = EventModel::TkrRecon::Event + "/TkrTruncatedPlane";
std::string EventModel::TkrRecon::TkrTruncationInfo          = EventModel::TkrRecon::Event + "/TkrTruncationInfo";
std::string EventModel::TkrRecon::TkrBoundBoxCol             = EventModel::TkrRecon::Event + "/TkrBoundBoxCol";
std::string EventModel::TkrRecon::TkrBoundBoxLinksCol        = EventModel::TkrRecon::Event + "/TkrBoundBoxLinksCol";
std::string EventModel::TkrRecon::TkrBoundBoxPointsCol       = EventModel::TkrRecon::Event + "/TkrBoundBoxPointsCol";
std::string EventModel::TkrRecon::TkrBoundBoxPointsToBoxTab  = EventModel::TkrRecon::Event + "/TkrBoundBoxPointsToBoxTab";
std::string EventModel::TkrRecon::TkrVecPointCol             = EventModel::TkrRecon::Event + "/TkrVecPointCol";
std::string EventModel::TkrRecon::TkrVecPointInfo            = EventModel::TkrRecon::Event + "/TkrVecPointInfo";
std::string EventModel::TkrRecon::TkrVecPointsLinkCol        = EventModel::TkrRecon::Event + "/TkrVecPointsLinkCol";
std::string EventModel::TkrRecon::TkrVecPointsLinkInfo       = EventModel::TkrRecon::Event + "/TkrVecPointsLinkInfo";
std::string EventModel::TkrRecon::TkrTrackElementsCol        = EventModel::TkrRecon::Event + "/TkrTrackElementsCol";
std::string EventModel::TkrRecon::TkrTrackElemsToLinksTab    = EventModel::TkrRecon::Event + "/TkrTrackElemsToLinksTab";
std::string EventModel::TkrRecon::TkrDiagnosticFlag          = EventModel::TkrRecon::Event + "/TkrDiagnosticFlag";

//reconstructed Cal data
std::string EventModel::CalRecon::Event               = EventModel::EventHeader + "/CalRecon";
std::string EventModel::CalRecon::CalXtalRecCol       = EventModel::CalRecon::Event + "/CalXtalRecCol";
std::string EventModel::CalRecon::CalMipTrackCol      = EventModel::CalRecon::Event + "/CalMipTrackCol";
std::string EventModel::CalRecon::CalClusterCol       = EventModel::CalRecon::Event + "/CalClusterCol";
std::string EventModel::CalRecon::CalClusterMap       = EventModel::CalRecon::Event + "/CalClusterMap";
std::string EventModel::CalRecon::CalRawClusterVec    = EventModel::CalRecon::Event + "/CalRawClusterVec";
std::string EventModel::CalRecon::CalUberCluster      = EventModel::CalRecon::Event + "/CalUberCluster";
std::string EventModel::CalRecon::CalUber2Cluster     = EventModel::CalRecon::Event + "/CalUber2Cluster";
std::string EventModel::CalRecon::CalEventEnergyCol   = EventModel::CalRecon::Event + "/CalEventEnergyCol";
std::string EventModel::CalRecon::CalEventEnergyMap   = EventModel::CalRecon::Event + "/CalEventEnergyMap";
std::string EventModel::CalRecon::CalClusterHitTab    = EventModel::CalRecon::Event + "/CalClusterHitTab";
std::string EventModel::CalRecon::GcrXtalCol          = EventModel::CalRecon::Event + "/GcrXtalCol";
std::string EventModel::CalRecon::GcrReconVals        = EventModel::CalRecon::Event + "/GcrReconVals";
std::string EventModel::CalRecon::GcrSelectedXtalsCol = EventModel::CalRecon::Event + "/GcrSelectedXtalsCol";
std::string EventModel::CalRecon::GcrTrack            = EventModel::CalRecon::Event + "/GcrTrack";
std::string EventModel::CalRecon::GcrSelectVals       = EventModel::CalRecon::Event + "/GcrSelectVals";

// reconstructed ACD data
std::string EventModel::AcdRecon::Event               = EventModel::EventHeader + "/AcdRecon";

std::string EventModel::AcdReconV2::Event             = EventModel::EventHeader + "/AcdReconV2";

