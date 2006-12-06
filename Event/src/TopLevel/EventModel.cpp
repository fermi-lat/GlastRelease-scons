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


// Digi event
std::string EventModel::Digi::Event                   = EventModel::EventHeader + "/Digi";
std::string EventModel::Digi::AcdDigiCol              = EventModel::Digi::Event + "/AcdDigiCol";
std::string EventModel::Digi::TkrDigiCol              = EventModel::Digi::Event + "/TkrDigiCol";
std::string EventModel::Digi::CalDigiCol              = EventModel::Digi::Event + "/CalDigiCol";
std::string EventModel::Digi::CalDigiHitTab           = EventModel::Digi::Event + "/CalDigiHitTab";
std::string EventModel::Digi::TkrDigiHitTab           = EventModel::Digi::Event + "/TkrDigiHitTab";
std::string EventModel::Digi::TkrClusterHitTab        = EventModel::Digi::Event + "/TkrClusterHitTab";
      
// reconstructed data (Tracker)
std::string EventModel::TkrRecon::Event               = EventModel::EventHeader + "/TkrRecon";
std::string EventModel::TkrRecon::TkrClusterCol       = EventModel::TkrRecon::Event + "/TkrClusterCol";
std::string EventModel::TkrRecon::TkrIdClusterMMap    = EventModel::TkrRecon::Event + "/TkrIdClusterMMap";
std::string EventModel::TkrRecon::TkrIdClusterMap     = EventModel::TkrRecon::Event + "/TkrIdClusterMap";
std::string EventModel::TkrRecon::TkrTrackCol         = EventModel::TkrRecon::Event + "/TkrTrackCol";
std::string EventModel::TkrRecon::TkrTrackHitCol      = EventModel::TkrRecon::Event + "/TkrTrackHitCol";
std::string EventModel::TkrRecon::TkrVertexCol        = EventModel::TkrRecon::Event + "/TkrVertexCol";
std::string EventModel::TkrRecon::TkrDiagnostics      = EventModel::TkrRecon::Event + "/TkrDiagnostics";
std::string EventModel::TkrRecon::TkrEventParams      = EventModel::TkrRecon::Event + "/TkrEventParams";
std::string EventModel::TkrRecon::TkrTruncatedPlane   = EventModel::TkrRecon::Event + "/TkrTruncatedPlane";
std::string EventModel::TkrRecon::TkrTruncationInfo   = EventModel::TkrRecon::Event + "/TkrTruncationInfo";

//reconstructed Cal data
std::string EventModel::CalRecon::Event               = EventModel::EventHeader + "/CalRecon";
std::string EventModel::CalRecon::CalXtalRecCol       = EventModel::CalRecon::Event + "/CalXtalRecCol";
//@@@FP 07/09/05
//std::string EventModel::CalRecon::CalMIPsCol          = EventModel::CalRecon::Event + "/CalMIPsCol";
//@@@FP 07/09/05
std::string EventModel::CalRecon::CalMipTrackCol      = EventModel::CalRecon::Event + "/CalMipTrackCol";
std::string EventModel::CalRecon::CalClusterCol       = EventModel::CalRecon::Event + "/CalClusterCol";
std::string EventModel::CalRecon::CalEventEnergyCol   = EventModel::CalRecon::Event + "/CalEventEnergyCol";
std::string EventModel::CalRecon::CalClusterHitTab    = EventModel::CalRecon::Event + "/CalClusterHitTab";
//@@@FP 07/09/05
//std::string EventModel::CalRecon::CalXtalMIPsTab      = EventModel::CalRecon::Event + "/CalXtalMIPsTab";
//@@@FP 07/09/05
//@@@CL 01/06/06 BEGIN:
std::string EventModel::CalRecon::GcrXtalCol      = EventModel::CalRecon::Event + "/GcrXtalCol";
//@@@CL 01/06/06 END
//@@@CL 06/27/06 BEGIN:
std::string EventModel::CalRecon::GcrSelectedXtalsCol      = EventModel::CalRecon::Event + "/GcrSelectedXtalsCol";
//@@@CL 06/27/06 END
std::string EventModel::CalRecon::GcrTrack      = EventModel::CalRecon::Event + "/GcrTrack";
std::string EventModel::CalRecon::GcrSelectVals     = EventModel::CalRecon::Event + "/GcrSelectVals";

// reconstructed ACD data
std::string EventModel::AcdRecon::Event               = EventModel::EventHeader + "/AcdRecon";

