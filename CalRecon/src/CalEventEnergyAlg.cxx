
#include "src/Utilities/CalException.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/TopLevel/EventModel.h"

#include <CalRecon/ICalEnergyCorr.h>
#include <CalRecon/ICalReconSvc.h>

/**   
* @class CalEventEnergyAlg
*
* @brief An algorithm for controlling and applying the various energy correction tools
*        used to determine the final event energy for GLAST
* 
* $Header$
*/


class CalEventEnergyAlg : public Algorithm
{
public:
    //! constructor
    CalEventEnergyAlg( const std::string & name, ISvcLocator * pSvcLocator ); 
    
    //! destructor
    virtual ~CalEventEnergyAlg() 
    {
        m_corrTools.clear();
    }; 
    
    virtual StatusCode initialize();

    StatusCode execute();

    StatusCode finalize() ;

private:

    //! correction tool names
    StringArrayProperty          m_corrToolNames ;
    
    //! correction tools
    std::vector<ICalEnergyCorr*> m_corrTools ;

    //! Type of correction for primary corrected parameters
    std::string                  m_corrType;

    //! Set status bits depending on which iteration of algorithm
    unsigned int                 m_passBits;
    
    //! package service
    ICalReconSvc *      m_calReconSvc ;
} ;

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_ALGORITHM_FACTORY(CalEventEnergyAlg) ;


CalEventEnergyAlg::CalEventEnergyAlg( const std::string & name, ISvcLocator * pSvcLocator )
 : Algorithm(name,pSvcLocator)
{   
    std::vector<std::string> corrToolVec;
    std::string              corrType;

    // Initial/default parameters depend on whether this would be the first pass
    // through CalRecon or the second. The first pass through should have the 
    // algorithm name "RawEnergy" assigned to it, the second (or subsequent) name
    // doesn't matter for initialization.
    // This branch if first pass
    if (name == "RawEnergy")
    {
        corrToolVec.push_back("CalRawEnergyTool");

        //corrType   = Event::CalCorToolResult::RAWENERGY;
        corrType   = "CalRawEnergyTool";
        m_passBits = Event::CalEventEnergy::PASS_ONE;
    }
    // This the default for the second iteration of CalRecon
    else
    {
		corrToolVec.push_back("CalValsCorrTool");
        corrToolVec.push_back("CalProfileTool");
        corrToolVec.push_back("CalLastLayerLikelihoodTool");
        corrToolVec.push_back("CalFullProfileTool");
        corrToolVec.push_back("CalTkrLikelihoodTool");
        corrToolVec.push_back("CalTransvOffsetTool");

        //corrType   = Event::CalCorToolResult::CALVALS;
        corrType   = "CalValsCorrTool";
        m_passBits = Event::CalEventEnergy::PASS_TWO;
    }

    // Declare the properties with these defaults
    declareProperty("corrToolNames", m_corrToolNames = corrToolVec);
    declareProperty("finalCorrType", m_corrType      = corrType);
}



StatusCode CalEventEnergyAlg::initialize()
{
    // Purpose and Method: Initialize the algorithm:
    //                     - Initializes the vector of pointers to the various
    //                       energy correction tools to all for this particular
    //                       iteration of reconstruction

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Basic initialization first
    log << MSG::INFO << "CalEventEnergyAlg Initialization";
    if( (sc = setProperties()).isFailure()) 
    {
        log << " didn't work!" << endreq;
        return sc;
    }
    log << endreq;
        
    if (service("CalReconSvc",m_calReconSvc,true).isFailure())
    {
        log<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
        return StatusCode::FAILURE ;
    }

    // Now build the list of correction tools to apply during execution
    const std::vector< std::string >& corrToolNames = m_corrToolNames ;
    std::vector< std::string >::const_iterator toolName ;
    ICalEnergyCorr* tool ;
    for (toolName = corrToolNames.begin(); toolName != corrToolNames.end(); toolName++) 
    {
        // Attempt to retrieve the tool
        sc = toolSvc()->retrieveTool(*toolName,tool);

        // Did we fail to find the tool?
        if (sc.isFailure() ) 
        {
            log << MSG::ERROR << " Unable to create " << *toolName << endreq ;
            return sc ;
        }
        // Otherwise add to our list
        else 
        {
            m_corrTools.push_back(tool) ;
        }
    }

    return sc;
}


StatusCode CalEventEnergyAlg::execute()
{
    //Purpose and method: Primary driver for running the various energy correction tools
    //                    Also creates and registers the output TDS object for the Event Energy
    //                    and retrieves the "best" energy as the one to use this event.
    // 
    // TDS input:  CalClusterCol
    // TDS output: CalEventEnergy (with a vector of CalCorToolResult objects for each tool run)
    //
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    log<<MSG::DEBUG<<"Begin execute()"<<endreq ;

    // Retrieve our TDS objects, we use Clusters to output corrected energy in CalEventEnergy
    Event::CalClusterCol*  calClusters = SmartDataPtr<Event::CalClusterCol>(eventSvc(),EventModel::CalRecon::CalClusterCol);
    Event::CalEventEnergy* calEnergy   = SmartDataPtr<Event::CalEventEnergy>(eventSvc(),EventModel::CalRecon::CalEventEnergy);
    Event::TkrVertexCol*   tkrVertices = SmartDataPtr<Event::TkrVertexCol>(eventSvc(),EventModel::TkrRecon::TkrVertexCol);

    // If no CalEnergy object (yet) then create one and register in TDS
    if (calEnergy == 0) {
        calEnergy = new Event::CalEventEnergy();

        if ((eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergy, calEnergy)).isFailure())
        {
            log<<MSG::ERROR<<"Cannot register CalEventEnergy"<<endreq ;
            return StatusCode::FAILURE ;
        }
    // Else reset CalEventEnergy
    } else {
        calEnergy->clear() ;
        calEnergy->clearStatusBit((Event::CalEventEnergy::StatusBits)calEnergy->getStatusBits()) ;
    }
    
    // No clusters no work
    if (calClusters->size() > 0)
    {
        // Set pointer to first vertex if it exists
        Event::TkrVertex* vertex = 0;

        if (tkrVertices != 0 && !tkrVertices->empty()) vertex = tkrVertices->front();

        // apply corrections according to vector of tools
        std::vector<ICalEnergyCorr *>::const_iterator tool ;
        for ( tool = m_corrTools.begin(); tool != m_corrTools.end(); ++tool ) {

            try {

                log<<MSG::DEBUG<<(*tool)->type()<<endreq ;
    
                // Loop over clusters 	 
                for ( Event::CalClusterCol::const_iterator cluster = calClusters->begin(); 	 
                      cluster != calClusters->end(); 	 
                      cluster++) { 	 
                    Event::CalCorToolResult* corResult = (*tool)->doEnergyCorr(*cluster, vertex); 	 
                    if (corResult != 0) {
                        calEnergy->push_back(corResult);
                        if(m_passBits != Event::CalEventEnergy::PASS_ONE) {
                            // Need set the status bit in the CalCluster  
                            (*cluster)->setStatusBit(Event::CalCluster::ENERGYCORR);
                        }
                    }
                }
                
            } catch( CalException & e ) {
                sc = sc && m_calReconSvc->handleError(name()+" CalException",(*tool)->type()+", "+e.what()) ;
            } catch( std::exception & e) {
                sc = sc && m_calReconSvc->handleError(name()+" std::exception",(*tool)->type()+", "+e.what()) ;
            } catch(...) {
                sc = sc && m_calReconSvc->handleError(name(),(*tool)->type()+", "+"unknown exception") ;
            }
        }
        
        // Set the pass number bits
        calEnergy->setStatusBit(m_passBits);

        // Go through and pick the "best" energy
        // For now this consists of choosing the value from CalValsCorrTool
        // It is envisaged that this will be replaced by a call to a tool/method which examines all 
        // possibilites and chooses the best for the given event
        for(Event::CalCorToolResultCol::iterator corIter = calEnergy->begin(); corIter != calEnergy->end(); corIter++)
        {
            Event::CalCorToolResult* corResult = *corIter;

            //if ((corResult->getStatusBits() & m_corrType) == m_corrType)
            if (corResult->getCorrectionName() == m_corrType)
            {
                // Set the main energy parameters to the selected correction tool values
                // Should also be setting a bit in CalEventEnergy to describe this?
                calEnergy->setParams(corResult->getParams());
                break;
            }
        }
    }
    
    log<<MSG::DEBUG<<"End execute()"<<endreq ;
    return sc;
}

StatusCode CalEventEnergyAlg::finalize()
{ 
    return StatusCode::SUCCESS ; 
}




