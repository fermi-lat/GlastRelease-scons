
#include "SingleClusterTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"


static const ToolFactory<SingleClusterTool>  s_factory;
const IToolFactory& SingleClusterToolFactory = s_factory;


SingleClusterTool::SingleClusterTool( const std::string& type, 
								   const std::string& name, 
								   const IInterface* parent)
								   : AlgTool(type,name,parent)
{
	// declare base interface for all consecutive concrete classes
	declareInterface<ICluster>(this);
    // Declare the properties that may be set in the job options file
	
}



StatusCode SingleClusterTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "Initializing SingleClusterTool" <<endreq;
    
    
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
    } else m_CalnLayers = int(value);
    
    if(!detSvc->getNumericConstByName(std::string("CsIWidth"),&m_CsIWidth))
    {
        log << MSG::ERROR << " constant " << " CsIWidth "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    if(!detSvc->getNumericConstByName(std::string("CsIHeight"),&m_CsIHeight))
    {
        log << MSG::ERROR << " constant " << " CsIHeight "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
         
    return sc;
}


StatusCode SingleClusterTool::findClusters(Event::CalXtalRecCol* calXtalRecCol)

//Purpose and method:
//
//   This function performs the calorimeter cluster reconstruction.
//   The main actions are:
//      - calculate energy sum
//                  energy per layer
//                  average position per layer
//                  quadratic spread per layer
//      - fit the particle direction using Fit_Direction() function
//      - calculate particle energy by profile fitting method
//          using Profile() function
//      - calculate particle energy by last layer correction method
//          using Leak() function
//      - store all calculated quantities in CalCluster object
// 
// TDS input: CalXtalRecCol
// TDS output: CalClustersCol


{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //get pointers to the TDS data structures

    const Point p0(0.,0.,0.);  
    Vector zero(0.,0.,0.);
    
    double ene = 0;                 // total energy in calorimeter
    Vector pCluster = p0;            // cluster position
    int nLayers = m_CalnLayers;     // number of layers
    
    // energy per layer
    std::vector<double> eneLayer(nLayers,0.);   
    
    // Vector of average position per layer
    std::vector<Vector> pLayer(nLayers);    
    
    
    // Vector of quadratic spread per layer
    std::vector<Vector> rmsLayer(nLayers);  
    
    
    
    // Compute barycenter and various moments
    
    
    // loop over all hitted crystals in CalXtalRecCol
    for (Event::CalXtalRecCol::const_iterator it = calXtalRecCol->begin();
    it != calXtalRecCol->end(); it++){

        // get pointer to teh reconstructed data for given crystal
		Event::CalXtalRecData* recData = *it;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy(); // crystal energy
        Vector pXtal = recData->getPosition() - p0; // Vector of crystal position
        int layer = (recData->getPackedId()).getLayer(); // layer number
        
        // update energy of corresponding layer
        eneLayer[layer]+=eneXtal;
        
        // update average position of corresponding layer
        Vector ptmp = eneXtal*pXtal;
        pLayer[layer] += ptmp;
        
        // Vector containing squared coordinates, weighted by crystal energy 
        Vector ptmp_sqr(ptmp.x()*pXtal.x(),
                        ptmp.y()*pXtal.y(),
                        ptmp.z()*pXtal.z());
 
        
        
        // update quadratic spread, which is proportional to eneXtal;
        // this means, that position error in one crystal
        // is assumed to be 1/sqrt(eneXtal) 
        rmsLayer[layer] += ptmp_sqr;
        
        
        // update energy sum
        ene  += eneXtal;

        // update cluster position
        pCluster += ptmp;
    }
    
    // Now take the means

    // if energy sum is not zero - normalize cluster position
    if(ene>0.)pCluster *= (1./ene); 

    // if energy is zero - set cluster position to non-physical value
    else pCluster=Vector(-1000., -1000., -1000.);
    int i = 0;
    
    // loop over calorimeter layers
    for( ;i<nLayers;i++){

        // if energy in the layer is not zero - finalize calculations
        if(eneLayer[i]>0)
        {
            // normalize position in the layer
            pLayer[i] *= (1./eneLayer[i]); 
            
            // normalize quadratic spread in the laye
            rmsLayer[i] *= (1./eneLayer[i]);
            
            
            // Vector containing the squared average position in each component
            Vector sqrLayer(pLayer[i].x()*pLayer[i].x(),
                            pLayer[i].y()*pLayer[i].y(),
                            pLayer[i].z()*pLayer[i].z());
            
            
            // the precision of transverse coordinate measurement
            // if there is no fluctuations: 1/sqrt(12) of crystal width
            Vector d; 
            if(i%2 == 1) d = Vector(m_CsIWidth*m_CsIWidth/12.,0.,0.);
            else d = Vector(0.,m_CsIWidth*m_CsIWidth/12.,0.);
            
            // subtracting the  squared average position and adding
            // the square of crystal width, divided by 12
            rmsLayer[i] += d-sqrLayer;
            
        }
        
        // if energy in the layer is zero - reset position and spread Vectors
        else 
        {
            pLayer[i]=p0;
            rmsLayer[i]=p0;
        }
    }
    
    
    // Now sum the different rms to have one transverse and one longitudinal rms
    double rms_trans=0;
    double rms_long=0;
    std::vector<Vector> posrel(nLayers);
    
    
    for(int ilayer=0;ilayer<nLayers;ilayer++)
    {
        
        posrel[ilayer]=pLayer[ilayer]-pCluster;
        
        // Sum alternatively the rms
        if(ilayer%2)
        {
            rms_trans += rmsLayer[ilayer].y();
            rms_long += rmsLayer[ilayer].x();
        }
        else
        {
            rms_trans += rmsLayer[ilayer].x();
            rms_long += rmsLayer[ilayer].y();
        }
    }
 

	    // Fill CalCluster data
    Event::CalCluster* cl = new Event::CalCluster(ene,pCluster+p0);

    cl->initialize(ene,
		   eneLayer,
		   pLayer,
		   rmsLayer,
		   rms_long,
		   rms_trans,
		   zero,
		   0.);

    getClusterCol()->add(cl);

    return sc;
}

StatusCode SingleClusterTool::finalize()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

StatusCode SingleClusterTool::execute()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}



