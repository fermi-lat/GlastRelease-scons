#include "CalLikelihoodTool.h"
#include "GaudiKernel/GaudiException.h" 
#include "src/Utilities/CalException.h" 

/**   
* @class CalLastLayerLikelihoodTool
*
* Algorithm for correction of energy leak through the bottom of the CAL using
* the correlation between that vallue and the energy deposit in the last layer.
*
*/

class CalLastLayerLikelihoodTool : public CalLikelihoodTool 
{
public:
    //! constructor
    CalLastLayerLikelihoodTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent);
    virtual ~CalLastLayerLikelihoodTool(){}
    StatusCode initialize();
    
    //! Energy leak correction using correlation to the energy in the last layer
    /*! This method uses the correlation between the energy lost in through the 
    * botton of the CAL and the energy deposited in the last layer of the 
    * calorimeter.
    * We used the Monte Carlo simulation of the LAT to determine this 
    * correlation at several energies, from 200 MeV up to 50 GeV,
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
private:
    int m_calNLayers;
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalLastLayerLikelihoodTool) ;
CalLastLayerLikelihoodTool::CalLastLayerLikelihoodTool( const std::string& type, 
                                            const std::string& name, 
                                            const IInterface* parent)
                                          : CalLikelihoodTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<ICalEnergyCorr>(this);
    declareProperty("dataFile",
                    m_dataFile="$(CALRECONXMLPATH)/CalLastLayerLikelihood.data");
};

StatusCode CalLastLayerLikelihoodTool::initialize()
{
  StatusCode sc= CalLikelihoodTool::initialize();
  if( sc==StatusCode::SUCCESS )
  {
    if (!m_detSvc->getNumericConstByName(std::string("CALnLayer"), &m_calNLayers)) 
    { 
        throw GaudiException("GlastDetSvc cannot find [CALnLayer]", name(), StatusCode::FAILURE);
    }
  }
  return sc;
}


Event::CalCorToolResult* CalLastLayerLikelihoodTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrTree* tree)
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
    MsgStream log(msgSvc(), "CalLastLayerLikelihoodTool::doEnergyCorr");
    Event::CalCorToolResult* corResult = 0;

    // if reconstructed tracker data doesn't exist or number of tracks is 0:
    if (tree == 0)
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
    // CALCULATE
    //    - get number of hits in TKR
    double calE7= 0.;
    if( hasCorrelation() )
    {
      calE7= (*cluster)[m_calNLayers-1].getEnergy();
    }


    setEventPDFdata(vertexPos+(geometricCut>.5)*16);
    setEventPDFparameters(fabs(trackDirection.z()),
                          cluster->getMomParams().getEnergy(),
                          calE7);
    log << MSG::VERBOSE 
        << "PDF Index: " << vertexPos+(geometricCut>.5)*16 
        << endreq
        << "Parameters: " << fabs(trackDirection.z()) << ", " 
        << cluster->getMomParams().getEnergy() << ", " << calE7 << endreq;
    
    corResult= calculateEvent(cluster, log);
    if( corResult ) 
      (*corResult)["GeometricCut"] = geometricCut ;


    log << MSG::DEBUG << "Ending doEnergyCorr: Reconstruction Done" << endreq;
    return corResult;
}

