#include "CalLikelihoodTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Event/TopLevel/EventModel.h"
#include "Event/Digi/TkrDigi.h"

/**   
* @class CalTkrLikelihoodTool
*
* Algorithm for correction of energy degradation in the tracker by correlating
* the energy in the CAL with the number of hit TKR strips.  
*
*/

class CalTkrLikelihoodTool : public CalLikelihoodTool 
{
public:
    //! constructor
    CalTkrLikelihoodTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent);
    virtual ~CalTkrLikelihoodTool(){}
    
    //! Tracker energy degradation correction using number of TKR hit strips
    /*! This method uses the correlation between the energy \"lost\" in the 
    * tracker and the energy deposited in the calorimeter.
    * We used the Monte Carlo simulation of the LAT to determine this 
    * correlation at several energies, from 50 MeV up to 3 GeV, 
    * and angles from 0 to 32\deg. 
    * See CalLikelihoodTool.calculateEvent for more information
    *
    * \par The method takes 2 arguments:
    * \param CalCluster
    * \param TkrVertex
    *
    *\return CalCorToolResult with an energy and error estimate.
    *
    *\author
    */
           
    Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrTree* );
};

#include <GaudiKernel/DeclareFactoryEntries.h>
DECLARE_TOOL_FACTORY(CalTkrLikelihoodTool) ;

CalTkrLikelihoodTool::CalTkrLikelihoodTool( const std::string& type, 
                                            const std::string& name, 
                                            const IInterface* parent)
                                          : CalLikelihoodTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<ICalEnergyCorr>(this);
    declareProperty("dataFile",
                    m_dataFile="$(CALRECONXMLPATH)/CalTkrLikelihood.data");
};

Event::CalCorToolResult* CalTkrLikelihoodTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrTree* tree)
//Purpose and method:
//
//   This function performs:
//   The main actions are:
//      - check wheter the event meets basic requirements (CUTS)
//      - calculate energy by TKR correction method using LikelihoodTool 
// 
// TDS input: CalCluster
// TDS output: CalClusters
{
    MsgStream log(msgSvc(), "CalTkrLikelihoodTool::doEnergyCorr");
    Event::CalCorToolResult* corResult = 0;

    // CALCULATE
    //    - get number of hits in TKR
    int nHits= 0;
    if( hasCorrelation() )
    {
      Event::TkrDigiCol *tkrDigiData =
                  SmartDataPtr<Event::TkrDigiCol>(m_dataSvc, 
                                                  EventModel::Digi::TkrDigiCol);

      for ( Event::TkrDigiCol::iterator it= tkrDigiData->begin(); 
            it!=tkrDigiData->end(); ++it )
        nHits+= (*it)->getNumHits();
    }

    // if reconstructed tracker data doesn't exist or number of tracks is 0:
    if (!tree)
    {
        log << MSG::DEBUG << "Ending doEnergyCorr: No TKR Reconstruction" 
            << endreq;
        return corResult;
    }

    if (!cluster)
    {
        log << MSG::DEBUG << "Ending doEnergyCorr: No Cluster" 
            << endreq;
        return corResult;
    }

    const Vector& trackDirection = tree->getAxisParams()->getEventAxis();
    const Point&  trackPosition  = tree->getAxisParams()->getEventPosition();
    // CUTS
    // this checks whether a set of PDF parameters exist for this event's
    // energy, direction and point of impact.
  
    // Energy:
    if( cluster->getMomParams().getEnergy() < minTrialEnergy()*.1 ) 
    {
      log << MSG::DEBUG << "Ending doEnergyCorr: "
                           "CAL Energy below Method Minimum"
          << endreq;
        return corResult;
    }

    if( cluster->getMomParams().getEnergy() > maxTrialEnergy() ) 
    {
       log << MSG::DEBUG << "Ending doEnergyCorr: "
                            "CAL Energy above Method Maximum"
           << endreq;
        return corResult;
    }
  
    // direction: slope must be above \f$cos(32\circ)$\f
    if( fabs(trackDirection.z()) < minSlope() )
    { 
        log << MSG::DEBUG << "Ending doEnergyCorr: Slope is too Small."
            << endreq;
        return corResult; 
    }


    int vertexPos = findTkrVertex(trackPosition);
    if( vertexPos < 0 || vertexPos > 15 ) 
    { 
        log << MSG::DEBUG << "Ending doEnergyCorr: "
                             "Vertex is out of Method Reach." 
        << endreq;
        return corResult; 
    }

    double geometricCut= findGeometricCut(trackPosition, trackDirection, 
                                          cluster);
    if( geometricCut<.15 )
    {
        log << MSG::DEBUG << "Ending doEnergyCorr: "
                             "Geometric Cut too low." 
        << endreq;
        return corResult;
    }

    if(!m_flight_geom)
      {
        nHits += 40; // 45.3*14/16
      }

    int iData= vertexPos;
    if( geometricCut>.4 ) iData+= 16*(geometricCut>.60)+32;
    else                  iData+= 16*(geometricCut>.25);

    setEventPDFdata(iData);
    setEventPDFparameters(fabs(trackDirection.z()),
                          cluster->getMomParams().getEnergy(),
                          nHits);
    log << MSG::VERBOSE 
        << "PDF Index: " << iData << endreq
        << "Parameters: " << fabs(trackDirection.z()) << ", " 
        << cluster->getMomParams().getEnergy() << ", " << nHits << endreq;
    
    corResult= calculateEvent(cluster, log);

    // Both next values are not calculated in the MeritTuple
    // in GlastRelease v6r8 MeritTuple.TkrTotalHits is actually only those
    // inside a radius around the best track   
    // These values are really of no interest for the rest of reconstruction and
    // the lines could be commented.
    if( corResult ) {
      (*corResult)["GeometricCut"] = geometricCut ;
      (*corResult)["TkrHitStrips"] = nHits ;
    }

    log << MSG::DEBUG << "Ending doEnergyCorr: Reconstruction Done" << endreq;
    return corResult;
}

