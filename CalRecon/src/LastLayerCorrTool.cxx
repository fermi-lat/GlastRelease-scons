
#include "LastLayerCorrTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

static const ToolFactory<LastLayerCorrTool>  s_factory;
const IToolFactory& LastLayerCorrToolFactory = s_factory;


LastLayerCorrTool::LastLayerCorrTool( const std::string& type, 
                                     const std::string& name, 
                                     const IInterface* parent)
                                     : AlgTool(type,name,parent)
{
	// declare base interface for all consecutive concrete classes
	declareInterface<IEnergyCorr>(this);
    // Declare the properties that may be set in the job options file
	
}



StatusCode LastLayerCorrTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "Initializing LastLayerCorrTool" <<endreq;

    IGlastDetSvc* detSvc;

        
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
    if(eTotal<200.) {
        setEnergyCorr(0.);
        return sc;
    }
    else
    {
        // Evaluation of energy using correlation method
        // Coefficients fitted using GlastSim.
        double p0 = -1.49 + 1.72*getTrackSlope();
        double p1 = 0.28 + 0.434 * getTrackSlope();
        double p2 = -15.16 + 11.55 * getTrackSlope();
        double p3 = 13.88 - 10.18 * getTrackSlope();
        double lnE = log(eTotal/1000.);
        double funcoef = (p0 + p1 * lnE )/(1+exp(-p2*(lnE - p3)));
        
        setEnergyCorr( funcoef * cluster->getEneLayer()[getNLayers()-1]);
        
        // Evaluation of energy using correlation with last layer
       	// coefficients fitted using tbsim and valid for ~1GeV<E<~50GeV
       	//double slope = 1.111 + 0.557*log(eTotal/1000.);
        //double intercept = 210. + 112.* log(eTotal/1000.)*log(eTotal/1000.); 
       	//double e_leak = slope * elast + intercept;
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

