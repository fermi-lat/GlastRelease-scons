
#include <CalRecon/ICalEnergyCorr.h>
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

/**   
* @class CalTransvOffsetTool
*
* @brief Correction of TransvOffset when tracker info available
*
*/


class CalTransvOffsetTool : public AlgTool, virtual public ICalEnergyCorr 
 {
  public :
    
    //! destructor
    CalTransvOffsetTool( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~CalTransvOffsetTool() {} ; 

    StatusCode initialize();
    
    Event::CalCorToolResult* doEnergyCorr(Event::CalClusterCol*, Event::TkrVertex* );

 private:

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc* m_dataSvc;
 } ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalTransvOffsetTool) ;

CalTransvOffsetTool::CalTransvOffsetTool(const std::string & type,
                                         const std::string & name,
                                         const IInterface * parent )
                                        : AlgTool(type,name,parent)
{
    declareInterface<ICalEnergyCorr>(this) ;
}

StatusCode CalTransvOffsetTool::initialize()
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

Event::CalCorToolResult* CalTransvOffsetTool::doEnergyCorr(Event::CalClusterCol* clusters, Event::TkrVertex* vertex)
{
    // calculating the transverse offset of average position in the calorimeter
    // with respect to the position predicted from tracker information
    Event::CalCorToolResult* corResult = 0;
    MsgStream log(msgSvc(), "CalTransvOffsetTool::doEnergyCorr");
    
    if (clusters->empty())
    {
        log << MSG::DEBUG << "Ending doEnergyCorr: No Cluster" 
            << endreq;
        return corResult;
    }
    Event::CalCluster * cluster = clusters->front() ;

    if (vertex != 0)
    {
        const Vector& trackDirection = vertex->getDirection();
        const Point&  trackPosition  = vertex->getPosition();

        // TU: I don't understand this calculation... 
        Vector calOffset     = (cluster->getPosition()) - trackPosition;
        double calLongOffset = trackDirection * calOffset;

        double calTransvOffset = sqrt(calOffset.mag2() - calLongOffset * calLongOffset);

        // Create a CalCorToolResult object to hold the information
        corResult = new Event::CalCorToolResult();

        corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
        corResult->setCorrectionName(type());
        corResult->setParams(cluster->getMomParams());
        corResult->setChiSquare(1.);
        (*corResult)["TransvOffset"] = calTransvOffset ;
    }
    
    return corResult;
}
 

