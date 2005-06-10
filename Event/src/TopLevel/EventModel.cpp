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

//reconstructed Cal data
std::string EventModel::CalRecon::Event               = EventModel::EventHeader + "/CalRecon";
std::string EventModel::CalRecon::CalXtalRecCol       = EventModel::CalRecon::Event + "/CalXtalRecCol";
std::string EventModel::CalRecon::CalMIPsCol          = EventModel::CalRecon::Event + "/CalMIPsCol";
std::string EventModel::CalRecon::CalMipTrackCol      = EventModel::CalRecon::Event + "/CalMipTrackCol";
std::string EventModel::CalRecon::CalClusterCol       = EventModel::CalRecon::Event + "/CalClusterCol";
std::string EventModel::CalRecon::CalEventEnergy      = EventModel::CalRecon::Event + "/CalEventEnergy";
std::string EventModel::CalRecon::CalClusterHitTab    = EventModel::CalRecon::Event + "/CalClusterHitTab";
std::string EventModel::CalRecon::CalXtalMIPsTab      = EventModel::CalRecon::Event + "/CalXtalMIPsTab";

// reconstructed ACD data
std::string EventModel::AcdRecon::Event               = EventModel::EventHeader + "/AcdRecon";

