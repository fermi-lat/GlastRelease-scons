
#include "LastLayerCorrTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
// to access an XML containing Digi parameters file
#include "xmlBase/IFile.h"

DECLARE_TOOL_FACTORY(LastLayerCorrTool) ;

LastLayerCorrTool::LastLayerCorrTool( const std::string& type, const std::string& name, const IInterface* parent)
:EnergyCorr(type,name,parent){

    // declare base interface for all consecutive concrete classes
    declareInterface<IEnergyCorr>(this);
    declareProperty ("xmlFile", m_xmlFile="$(CALRECONROOT)/xml/CalLayer.xml");
};


StatusCode LastLayerCorrTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc

{
    if (EnergyCorr::initialize().isFailure())
     { return StatusCode::FAILURE ; }

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "Initializing LastLayerCorrTool" <<endreq;
    
    // extracting detector geometry constants from xml file
       
 // Read in the parameters from the XML file
    xmlBase::IFile m_ifile(m_xmlFile.c_str());
    if (m_ifile.contains("lastlayer","c0") && m_ifile.contains("lastlayer","c1") && m_ifile.contains("lastlayer","c2") && m_ifile.contains("lastlayer","c3") && m_ifile.contains("lastlayer","b0") && m_ifile.contains("lastlayer","b1") ) {
        m_c0 = m_ifile.getDouble("lastlayer", "c0");
	    log << MSG::INFO << " value for c0 " << m_c0 << endreq;
        m_c1 = m_ifile.getDouble("lastlayer", "c1");
	    log << MSG::INFO << " value for c1 " << m_c1 << endreq;
        m_c2 = m_ifile.getDouble("lastlayer", "c2");
	    log << MSG::INFO << " value for c0 " << m_c2 << endreq;
        m_c3 = m_ifile.getDouble("lastlayer", "c3");
	    log << MSG::INFO << " value for c3 " << m_c3 << endreq;
        m_b0 = m_ifile.getDouble("lastlayer", "b0");
	    log << MSG::INFO << " value for b0 " << m_b0 << endreq;
        m_b1 = m_ifile.getDouble("lastlayer", "b1");
	    log << MSG::INFO << " value for b1 " << m_b1 << endreq;
        m_k0 = m_ifile.getDouble("lastlayer", "k0");
	    log << MSG::INFO << " value for k0 " << m_k0 << endreq;
        m_k1 = m_ifile.getDouble("lastlayer", "k1");
	    log << MSG::INFO << " value for k1 " << m_k1 << endreq;
        m_k2 = m_ifile.getDouble("lastlayer", "k2");
	    log << MSG::INFO << " value for k2 " << m_k2 << endreq;
    }
    else return StatusCode::FAILURE;

//    sc        = serviceLocator()->getService("EventDataSvc", iService);
//    m_dataSvc = dynamic_cast<DataSvc*>(iService);   

     
    return sc;
}


//Purpose and method:
//
//   This function performs 
//   The main actions are:
//      - calculate particle energy by last layer correction method
//          using Leak() function
// 
// TDS input: CalCluster
// TDS output: CalClusters

StatusCode LastLayerCorrTool::doEnergyCorr( Event::CalCluster * cluster )
 {

  MsgStream lm(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
   
  double eTotal = cluster->getEnergySum() ;
  double eCorrected = 0 ;
  
  if (eTotal<500.) {

    eCorrected = -1. ;

  }
  else {

    Event::TkrVertexCol * tkrRecData = SmartDataPtr<Event::TkrVertexCol>(getKernel()->getEventSvc(),
      EventModel::TkrRecon::TkrVertexCol);    
           
    // if reconstructed tracker data doesn't exist - put the debugging message
    if (tkrRecData == 0) {

      eCorrected = -2 ;
      //log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
      //return sc;

    }
    // if data exist and number of tracks not zero 
    // - get information of track position and direction 
    else {

	  // First get reconstructed direction from tracker
	  int ntracks = getKernel()->getTkrNVertices() ;
	  const Vector & trackDirection = getKernel()->getTkrFrontVertexDirection() ;
	  const Point & trackPosition = getKernel()->getTkrFrontVertexPosition() ;
      
      double thetarec=-1;//reconstructed theta
      double phirec=0;  //reconstructed phi
      double px=0;      //projected position on (0,0,0 plane) 
      double py=0;
        
      if ( ntracks > 0 ) {
      
	    thetarec=trackDirection[2];
	    phirec=atan2(trackDirection[1],trackDirection[0]);
	    
	    //find projected positions on the 0,0,0 plane
	    px = trackPosition[0]+sqrt((-trackPosition[2]/thetarec)*(-trackPosition[2]/thetarec)-trackPosition[2]*trackPosition[2])*cos(phirec) ; 
	    py = trackPosition[1]+sqrt((-trackPosition[2]/thetarec)*(-trackPosition[2]/thetarec)-trackPosition[2]*trackPosition[2])*sin(phirec) ; 
	    px = px+(m_k0-m_k1*thetarec+m_k2*thetarec*thetarec)*cos(phirec) ;
	    py = py+(m_k0-m_k1*thetarec+m_k2*thetarec*thetarec)*sin(phirec) ;

	    //log << MSG::INFO << "track direction = " << thetarec << endreq;
            
      }
      else {

	  //log << MSG::DEBUG << "No reconstructed tracks " << endreq;

      }	

      //for now valid up to 26 degrees.. 
      if (-thetarec < 0.898 ) {
          eCorrected = -3. ;
          //return sc ;
      }
      else {
          //make cuts for the 0-26 degree limit. Will later make the (theta,phi) dependence available
    
          if ( fabs(fabs(fabs(fabs(px)-375)-187.5)-187.5)>55 &&
               fabs(fabs(fabs(fabs(py)-375)-187.5)-187.5)>55) {
    
            // Evaluation of energy using correlation method
            // Coefficients fitted using GlastSim.
            /*
            double p0 = -1.49 + 1.72*getKernel()->getSlope();
            double p1 = 0.28 + 0.434 * getKernel()->getSlope();
            double p2 = -15.16 + 11.55 * getKernel()->getSlope();
            double p3 = 13.88 - 10.18 * getKernel()->getSlope();
            double lnE = log(eTotal/1000.);
            double funcoef = (p0 + p1 * lnE )/(1+exp(-p2*(lnE - p3)));
            */
          
            double acoef = m_c0 - m_c1*thetarec + m_c2*thetarec*thetarec;
    
            double lnE = log(eTotal/1000.) ;
            double funcoef = (acoef + m_c3 * lnE ) ;
            double E1 = eTotal + funcoef* cluster->getEneLayer()[getKernel()->getCalNLayers()-1];
            double biascoef = m_b0 + m_b1*log(E1/1000);
            double E2 = E1/biascoef;
            funcoef = (acoef + m_c3 * log(E2/1000) ) ;
            ///first iteration
            double E3 = eTotal + funcoef* cluster->getEneLayer()[getKernel()->getCalNLayers()-1];
            biascoef= m_b0 + m_b1*log(E3/1000);
            double E4 = E3/biascoef;
            funcoef = (acoef + m_c3 * log(E4/1000) );
            ///second iteration.. probably overkill
            eCorrected = E4 ;
    
            // eCorrected = funcoef * cluster->getEneLayer()[getNLayers()-1] ;
    
          }
          else {
          
            eCorrected = -4 ;
    
          }
       }
    }
  }

 cluster->initialize(eCorrected,
   cluster->getEneLayer(),
   cluster->getPosLayer(),
   cluster->getRmsLayer(),
   cluster->getRmsLong(),
   cluster->getRmsTrans(),
   cluster->getDirection(),
   cluster->getTransvOffset()) ;
  return sc;
}
