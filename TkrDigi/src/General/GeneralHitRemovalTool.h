/*
 * @class GeneralHitRemovalTool
 *
 * @brief Truncates the strip lists to account for finite controller and cable buffers.
 * Called from TkrDigiAlg
 *
 * @authors Leon Rochester
 *
 * $Header$
 */

#ifndef __GENERALHitRemovalTOOL_H__
#define __GENERALHitRemovalTOOL_H__

#include "../IHitRemovalTool.h"
#include "../SiLayerList.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <string>

class GeneralHitRemovalTool : public AlgTool, virtual public IHitRemovalTool {

 public:

    /// Standard Gaudi Tool interface constructor
    GeneralHitRemovalTool(const std::string&, const std::string&, const IInterface*);
    /// Initializes the tool
    StatusCode initialize();
    /// runs the tool
    StatusCode execute();

 private:

    /// Pointer to the event data service (aka "eventSvc")
    IDataProviderSvc* m_edSvc;
    /// Pointer to the Glast detector service
    IGlastDetSvc*     m_gdSvc;
    /// Pointer to TkrGeometrySvc
    ITkrGeometrySvc*  m_tkrGeom;
    /// flag to kill bad strips
    bool m_killBadStrips;
    /// flag to kill failed layers
    bool m_killFailed;

};

#endif
