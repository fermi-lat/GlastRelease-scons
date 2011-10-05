/** @file CalRawEnergyTool.cxx
@brief implementation of the class CalRawEnergyTool

$Header$

*/

#include <CalRecon/ICalEnergyCorr.h>
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "CLHEP/Matrix/Vector.h"

/**   
* @class CalRawEnergyTool
* @author CalRecon Rewrite Group
*
* This sets the "raw" energy for an event (after clustering)
*
* $Header$
*/

class CalRawEnergyTool : public AlgTool, virtual public ICalEnergyCorr 
{
public:

    //! Constructor
    CalRawEnergyTool( const std::string& type, const std::string& name, const IInterface* parent);
    //! destructor
    virtual ~CalRawEnergyTool() {}; 

    StatusCode initialize();

    // worker function to get the corrected energy      
    Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrTree* );

private:

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalRawEnergyTool) ;

CalRawEnergyTool::CalRawEnergyTool( const std::string & type, 
                                    const std::string & name, 
                                    const IInterface * parent )
                                    : AlgTool(type,name,parent)
{ 
    declareInterface<ICalEnergyCorr>(this) ; 
}

StatusCode CalRawEnergyTool::initialize()
{
    // This function does following initialization actions:
    //    - Sets up a pointer to the data service as a convenience for the tool
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing CalRawEnergyTool" <<endreq;

    //Locate and store a pointer to the data service which allows access to the TDS
    if ((sc = service("EventDataSvc", m_dataSvc)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }

    return sc;
}


Event::CalCorToolResult* CalRawEnergyTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrTree* )
{
    //Purpose and method:
    //
    //   This function determines the overal "raw" energy for an event
    //   Note: This sums the energy contribution from all clusters and
    //         does a weighted average for centroid position and axis
    // 
    // TDS input:  CalClusters
    // TDS output: CalEventEnergy
    
    MsgStream lm(msgSvc(), name());

    // Set up to loop over all clusters to get total raw energy
    const Event::CalParams&  momParams  = cluster->getMomParams();

    // From those, create a copy for the corrected energy object
    Event::CalParams params = momParams;

    // Create a CalCorToolResult object to hold the information
    Event::CalCorToolResult* corResult = new Event::CalCorToolResult();

    corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
    corResult->setCorrectionName(type());
    corResult->setParams(params);
    corResult->setChiSquare(1.);

    return corResult;
}

