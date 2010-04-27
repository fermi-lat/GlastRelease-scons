// File and Version Information:
// $Header$

#define _OverlayEventModel_CPP_

// Define this in order to export symbols
#define OEVT_DLL_EXPORTS

#include "OverlayEvent/OverlayEventModel.h"
#include "GaudiKernel/ClassID.h"

std::string OverlayEventModel::OverlayEventHeader = "/Event";

/** 
 *  @class MC
 *
*/

// Digi Overlay 
std::string OverlayEventModel::Overlay::EventOverlay    = OverlayEventModel::OverlayEventHeader    + "/Overlay";
std::string OverlayEventModel::Overlay::SrcOverlay      = OverlayEventModel::Overlay::EventOverlay + "/Src";
std::string OverlayEventModel::Overlay::AcdOverlayCol   = OverlayEventModel::Overlay::EventOverlay + "/AcdOverlayCol";
std::string OverlayEventModel::Overlay::TkrOverlayCol   = OverlayEventModel::Overlay::EventOverlay + "/TkrOverlayCol";
std::string OverlayEventModel::Overlay::CalOverlayCol   = OverlayEventModel::Overlay::EventOverlay + "/CalOverlayCol";
std::string OverlayEventModel::Overlay::GemOverlay      = OverlayEventModel::Overlay::EventOverlay + "/GemOverlay";
std::string OverlayEventModel::Overlay::DiagDataOverlay = OverlayEventModel::Overlay::EventOverlay + "/DiagDataOverlay";
std::string OverlayEventModel::Overlay::PtOverlay       = OverlayEventModel::Overlay::EventOverlay + "/PtOverlay";

