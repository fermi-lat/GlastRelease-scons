
#include "LastLayerCorrTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/SmartDataPtr.h"
// to access an XML containing Digi parameters file
#include "xmlBase/IFile.h"

static const ToolFactory<LastLayerCorrTool>  s_factory;
const IToolFactory& LastLayerCorrToolFactory = s_factory;

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
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "Initializing LastLayerCorrTool" <<endreq;

    IGlastDetSvc* detSvc;

    IService*   iService = 0;       
    // get pointer to GlastDetSvc
    sc = service("GlastDetSvc", detSvc);
    
    // if GlastDetSvc isn't available - put error message and return
    if(sc.isFailure())
    {
        log << MSG::ERROR << "GlastDetSvc could not be found" <<endreq;
        return sc;
    }
    
    
    // extracting detector geometry constants from xml file
    
    double value;
    if(!detSvc->getNumericConstByName(std::string("CALnLayer"), &value)) 
    {
        log << MSG::ERROR << " constant " << " CALnLayer "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else setNLayers(int(value));
   
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

    sc        = serviceLocator()->getService("EventDataSvc", iService);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);   

     
    return sc;
}


StatusCode LastLayerCorrTool::doEnergyCorr(double eTotal, Event::CalCluster* cluster)

//Purpose and method:
//
//   This function performs 
//   The main actions are:
//      - calculate particle energy by last layer correction method
//          using Leak() function
// 
// TDS input: CalCluster
// TDS output: CalClusters

{

  MsgStream lm(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  
  if(eTotal<500.) {
    setEnergyCorr(-1.);
    return sc;
  }
  else
    {

      Event::TkrVertexCol *tkrRecData =  SmartDataPtr<Event::TkrVertexCol> (m_dataSvc,
      EventModel::TkrRecon::TkrVertexCol);    
           
      // if reconstructed tracker data doesn't exist - put the debugging message
      if (tkrRecData == 0) {
	setEnergyCorr(-2);
	//log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
         return sc;
      }
      
      
      // if data exist and number of tracks not zero 
      // - get information of track position and direction 
      else
	{

	  int rectkr=0;  
	  int ntracks=0;    
	  // First get reconstructed direction from tracker

	  ntracks = tkrRecData->size();
	  //log << MSG::DEBUG << "number of tracks = " << ntracks << endreq;

	  Vector trackDirection;
	  Point trackVertex;
	double slope=0;
	double thetarec=-1;//reconstructed theta
	double phirec=0;  //reconstructed phi
	double px=0;      //projected position on (0,0,0 plane) 
	double py=0;
        
        if (ntracks > 0) {
            rectkr++;
            trackDirection = tkrRecData->front()->getDirection();
            trackVertex = tkrRecData->front()->getPosition();
            slope = fabs(trackDirection.z());
	    thetarec=trackDirection[2];
	    phirec=atan2(trackDirection[1],trackDirection[0]);
	    
	    //find projected positions on the 0,0,0 plane
	    px = trackVertex[0]+sqrt((-trackVertex[2]/thetarec)*(-trackVertex[2]/thetarec)-trackVertex[2]*trackVertex[2])*cos(phirec); 
	   
	    py = trackVertex[1]+sqrt((-trackVertex[2]/thetarec)*(-trackVertex[2]/thetarec)-trackVertex[2]*trackVertex[2])*sin(phirec); 
	    px = px+(m_k0-m_k1*thetarec+m_k2*thetarec*thetarec)*cos(phirec);

	    py = py+(m_k0-m_k1*thetarec+m_k2*thetarec*thetarec)*sin(phirec);

	    //log << MSG::INFO << "track direction = " << thetarec << endreq;
            
        } else {
	  //log << MSG::DEBUG << "No reconstructed tracks " << endreq;
        }	

	//for now valid up to 26 degrees.. 
	if (-thetarec < 0.898 ) {setEnergyCorr(-3.);return sc;}

	//make cuts for the 0-26 degree limit. Will later make the (theta,phi) dependence available

	if (fabs(fabs(fabs(fabs(px)-375)-187.5)-187.5)>55 && fabs(fabs(fabs(fabs(py)-375)-187.5)-187.5)>55) {

        // Evaluation of energy using correlation method
        // Coefficients fitted using GlastSim.
	  /*
        double p0 = -1.49 + 1.72*getTrackSlope();
        double p1 = 0.28 + 0.434 * getTrackSlope();
        double p2 = -15.16 + 11.55 * getTrackSlope();
        double p3 = 13.88 - 10.18 * getTrackSlope();
        double lnE = log(eTotal/1000.);
        double funcoef = (p0 + p1 * lnE )/(1+exp(-p2*(lnE - p3)));
	  */
	  
	double acoef = m_c0 - m_c1*thetarec + m_c2*thetarec*thetarec;

        double lnE = log(eTotal/1000.);
        double funcoef = (acoef + m_c3 * lnE );
	double E1 = eTotal + funcoef* cluster->getEneLayer()[getNLayers()-1];
	double biascoef = m_b0 + m_b1*log(E1/1000);
	double E2 = E1/biascoef;
	funcoef = (acoef + m_c3 * log(E2/1000) );
	///first iteration
	double E3 = eTotal + funcoef* cluster->getEneLayer()[getNLayers()-1];
	biascoef= m_b0 + m_b1*log(E3/1000);
	double E4 = E3/biascoef;
	funcoef = (acoef + m_c3 * log(E4/1000) );
	///second iteration.. probably overkill
	setEnergyCorr( E4);

        //setEnergyCorr( funcoef * cluster->getEneLayer()[getNLayers()-1]);

	} else {setEnergyCorr(-4);}
      }
    }

    return sc;
}


StatusCode LastLayerCorrTool::finalize()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

StatusCode LastLayerCorrTool::execute()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

