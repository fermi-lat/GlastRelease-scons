
#include "src/Utilities/CalException.h"

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/TreeClusterRelation.h"
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

// DC: useless ?
//    //! Type of correction for primary corrected parameters
//    std::string                  m_corrType;
//
//    //! Set status bits depending on which iteration of algorithm
//    Event::CalEventEnergy::StatusBits m_passBits ;
//    
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
    if (name == "CalRawEnergy")
    {
        corrToolVec.push_back("CalRawEnergyTool");

// DC: redundant with correctionName
//        corrType   = Event::CalCorToolResult::RAWENERGY;
// DC: useless ?
//        corrType   = "CalRawEnergyTool";
//        m_passBits = Event::CalEventEnergy::PASS_ONE;
    }
    // This the default for the second iteration of CalRecon
    else
    {
      corrToolVec.push_back("CalValsCorrTool");
      // Ph. Bruel : old profile method, not used anymore 
      //        corrToolVec.push_back("CalProfileTool");
      // Philippe Bruel reminds us this isn't used anymore, so remove 9/22/09
      //corrToolVec.push_back("CalLastLayerLikelihoodTool");
      corrToolVec.push_back("CalFullProfileTool");
      // Philippe Bruel reminds us this isn't used anymore, so remove 9/22/09
      //corrToolVec.push_back("CalTkrLikelihoodTool");
      corrToolVec.push_back("CalLikelihoodManagerTool"); //POL
// DC: probably useless
//        corrToolVec.push_back("CalTransvOffsetTool");

        //corrType   = Event::CalCorToolResult::CALVALS;
// DC: useless ?
//        corrType   = "CalValsCorrTool";
//        m_passBits = Event::CalEventEnergy::PASS_TWO;
    }

    // Declare the properties with these defaults
    declareProperty("corrToolNames", m_corrToolNames = corrToolVec);
//    declareProperty("finalCorrType", m_corrType      = corrType);
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
    Event::CalClusterCol*     calClusterCol = SmartDataPtr<Event::CalClusterCol>(eventSvc(),EventModel::CalRecon::CalClusterCol);

    // See if we already have a calEnergyMap in the TDS
    Event::CalEventEnergyMap* calEnergyMap = SmartDataPtr<Event::CalEventEnergyMap>(eventSvc(),EventModel::CalRecon::CalEventEnergyMap);

    // If it is not there, then set up the CalEventEnergyMap in the TDS to relate clusters to energy algorithm
    if (!calEnergyMap)
    {
        calEnergyMap = new Event::CalEventEnergyMap();

        if ((eventSvc()->registerObject(EventModel::CalRecon::CalEventEnergyMap, calEnergyMap)).isFailure())
        {
            log << MSG::ERROR << "Cannot register CalEventEnergyMap in the TDS" << endreq;
            return StatusCode::FAILURE;
        }
    }

    // Now recover the mapping between clusters and trees from the TDS
    Event::ClusterToRelationMap* clusterToRelationMap = 
        SmartDataPtr<Event::ClusterToRelationMap>(eventSvc(), EventModel::Recon::ClusterToRelationMap);

    // there could be no clusters collection, if there was no hits
    if (calClusterCol && !calClusterCol->empty()) 
    {
        // Loop over all clusters in the list (including the uber cluster)
        for(Event::CalClusterCol::iterator clusItr = calClusterCol->begin(); clusItr != calClusterCol->end(); clusItr++)
        {
            // Recover the cluster pointer
            Event::CalCluster* cluster = *clusItr;

            // Search for the existing CalEventEnergy object
            Event::CalEventEnergyMap::iterator clusEnergyItr = calEnergyMap->find(cluster);

            // What we want...
            Event::CalEventEnergy * calEnergy;

            // Is there already one available?
            if (clusEnergyItr == calEnergyMap->end()) 
            {       
                calEnergy = new Event::CalEventEnergy;
                (*calEnergyMap)[cluster].push_back(calEnergy);
            } 
            else  calEnergy = clusEnergyItr->second.front();

            // Recover the pointer to the TkrTree, if there is one available
            Event::TkrTree* tree = 0;

            if (clusterToRelationMap)
            {
                // Recover the tree associated to this cluster, if there is one
                Event::ClusterToRelationMap::iterator clusTreeItr = clusterToRelationMap->find(cluster);

                if (clusTreeItr != clusterToRelationMap->end()) tree = clusTreeItr->second.front()->getTree();
            }
    
            // apply corrections according to vector of tools
            std::vector<ICalEnergyCorr *>::const_iterator tool ;
            for ( tool = m_corrTools.begin(); tool != m_corrTools.end(); ++tool ) 
            {
                try 
                {
                    log<<MSG::DEBUG<<(*tool)->type()<<endreq ;
    
                    // Call the correction tool         
                    Event::CalCorToolResult* corResult = (*tool)->doEnergyCorr(cluster, tree);          
                    if (corResult != 0) 
                    {
                        calEnergy->push_back(corResult);
                    }
                    
                } catch( CalException & e ) {
                    sc = sc && m_calReconSvc->handleError(name()+" CalException",(*tool)->type()+", "+e.what()) ;
                } catch( std::exception & e) {
                    sc = sc && m_calReconSvc->handleError(name()+" std::exception",(*tool)->type()+", "+e.what()) ;
                } catch(...) {
                    sc = sc && m_calReconSvc->handleError(name(),(*tool)->type()+", "+"unknown exception") ;
                }
            }
            
            // Go through and pick the "best" energy
            // For now this consists of choosing the value from CalValsCorrTool
            // It is envisaged that this will be replaced by a call to a tool/method which examines all 
            // possibilites and chooses the best for the given event
            
            const Event::CalCorToolResult * bestCorResult = calEnergy->findLast("CalValsCorrTool") ;
            if (bestCorResult==0) {
                bestCorResult = calEnergy->findLast("CalRawEnergyTool") ;
            }    
            if (bestCorResult!=0) {
                calEnergy->setParams(bestCorResult->getParams()) ;
                calEnergy->setStatusBits(Event::CalEventEnergy::VALIDPARAMS) ;
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




