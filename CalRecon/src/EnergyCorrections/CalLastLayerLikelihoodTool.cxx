#include "CalLikelihoodTool.h"
#include "GaudiKernel/GaudiException.h" 

/**   
* @class CalLastLayerLikelihoodTool
*
* Algorithm for correction of energy degradation in the tracker by correlating
* the energy in the CAL with the number of hit TKR strips.  
*
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
    
    //! Tracker energy degradation correction using number of TKR hit strips
    /*! This method uses the correlation between the energy \"lost\" in the 
    * tracker and the energy deposited in the calorimeter.
    * We used the Monte Carlo simulation of the LAT to determine this correlation
    * at several energies, from 50 MeV up to 1 GeV, and angles from 0 to 32\deg. 
    * For one particular incident energy and angle, the bidimensionnal
    * distribution of  the  number of hit strips and the energy deposited in the 
    * CAL can be characterised by the 1D distribution:
    * \f[ E_{CAL} + \alpha TkrHits \f]
    * where  \f$\alpha\f$ is been optimised so as to obtain the narrowest such
    * distribution, normalised to a probability and with its MPV at the incident
    * energy.
    * These distributions can be used to defined a probability density function.
    * The reconstructed energy for a given event then becomes the one maximising
    * the probability, for a reconstruced direction, CAL raw energy,...
    *
    * \par The method takes 4 arguments:
    * \param eTotal Total energy measured in the calorimeter in MeV
    * \param nHits  Total number of hit strips in the CAL
    * \param vertex[3] reconstructed vertex position  
    * \param dir[3] reconstructed direction
    *
    *\return Corrected energy in MeV
    *
    *\warning needs TKR reconstruction
    *
    *\author
    */
           
    Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrVertex* );
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
                    m_dataFile="$(CALRECONROOT)/xml/CalLastLayerLikelihood.data");
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


Event::CalCorToolResult* CalLastLayerLikelihoodTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrVertex* vertex)
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
        log << MSG::DEBUG << "Ending doEnergyCorr: " 
            << endreq;

    // if reconstructed tracker data doesn't exist or number of tracks is 0:
    if (vertex == 0)
    {
        log << MSG::DEBUG << "Ending doEnergyCorr: No TKR Reconstruction" 
            << endreq;
        return corResult;
    }

    const Vector& trackDirection = vertex->getDirection();
    const Point&  trackPosition  = vertex->getPosition();
    // CUTS
    // this checks whether a set of PDF parameters exist for this event's
    // energy, direction and point of impact.
  
    // Energy:
    if( cluster->getCalParams().getEnergy() < 80. ) 
    {
      log << MSG::DEBUG << "Ending doEnergyCorr: "
                           "CAL Energy below Method Minimum"
          << endreq;
        return corResult;
    }

    if( cluster->getCalParams().getEnergy() > 49000. ) 
    {
       log << MSG::DEBUG << "Ending doEnergyCorr: "
                            "CAL Energy above Method Maximum"
           << endreq;
        return corResult;
    }
  
    // direction: slope must be above \f$cos(32\circ)$\f
    if( fabs(trackDirection.z()) < .8480481 )
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
                          cluster->getCalParams().getEnergy(),
                          calE7);
    log << MSG::VERBOSE 
        << "PDF Index: " << vertexPos+(geometricCut>.5)*16 
        << endreq
        << "Parameters: " << fabs(trackDirection.z()) << ", " 
        << cluster->getCalParams().getEnergy() << ", " << calE7 << endreq;
    
    corResult = calculateEvent( 200., 50000., cluster, log );

    log << MSG::DEBUG << "Ending doEnergyCorr: Reconstruction Done" << endreq;
    return corResult;
}

