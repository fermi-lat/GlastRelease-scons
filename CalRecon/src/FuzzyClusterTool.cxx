// $Header$
#include "ICluster.h"
#include "FuzzyClusterTool.h"
#include "CLHEP/Matrix/Vector.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"


static const ToolFactory<FuzzyClusterTool>  s_factory;
const IToolFactory& FuzzyClusterToolFactory = s_factory;


FuzzyClusterTool::FuzzyClusterTool( const std::string& type, 
								   const std::string& name, 
								   const IInterface* parent)
								   : Cluster(type,name,parent)
{
	// declare base interface for all consecutive concrete classes
	declareInterface<ICluster>(this);
        // Declare the properties that may be set in the job options file
	declareProperty ("clusterTool", m_clusterToolName="FuzzyCluster");
	declareProperty ("clusterSetNo", m_clusterSetNo=0);	
}



StatusCode FuzzyClusterTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc
//    - retrieves the clustering tool

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
	log << MSG::INFO << "Initializing FuzzyClusterTool" <<endreq;
    
    
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
    
    // retrieving the clustering tool as private
    sc = toolSvc()->retrieveTool(m_clusterToolName,
    		m_fuzzyClusterTool, this);
    if(sc.isFailure()) {
    	log << MSG::ERROR << "Unable to create " << m_clusterToolName
    		<< endreq;
    	return sc;
    }

    return sc;
}


StatusCode FuzzyClusterTool::findClusters(Event::CalXtalRecCol* calXtalRecCol)

//Purpose and method:
//
//   This function performs the calorimeter cluster reconstruction.
//   The main actions are:
//      - calculate energy sum
//                  energy per layer
//                  average position per layer
//                  quadratic spread per layer
//      - fit the particle direction using Fit_Direction() function
//      - store all calculated quantities in CalCluster objects
// 
// TDS input: CalXtalRecCol
// TDS output: CalClustersCol


{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //get pointers to the TDS data structures

    const Point p0(0.,0.,0.);  
     
    int nLayers = m_CalnLayers;     // number of layers
    
    int nClusters = 1;

    // data set for the fuzzy cluster tool
    std::vector<HepVector> theDataSet; 
    
    // loop over all hitted crystals in CalXtalRecCol
    Event::CalXtalRecCol::const_iterator it;
    for (it = calXtalRecCol->begin(); it != calXtalRecCol->end(); it++){
        // get pointer to the reconstructed data for given crystal
		Event::CalXtalRecData* recData = *it;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy(); // crystal energy
        Vector pXtal = recData->getPosition() - p0; // Vector of crystal position
	HepVector dataPoint(4);
	dataPoint[0] = pXtal.x();
	dataPoint[1] = pXtal.y();
	dataPoint[2] = pXtal.z();
	dataPoint[3] = eneXtal;
	theDataSet.push_back(dataPoint);
    }

int dataSetSize = theDataSet.size();
int maxClusters = m_fuzzyClusterTool->getMaxClusters();

// if there are less than one point per cluster, consider one single cluster
if(dataSetSize >= maxClusters) {

    sc = m_fuzzyClusterTool->setData(theDataSet);
    if(sc.isFailure())
	return sc;

    sc = m_fuzzyClusterTool->perform();
    if(sc.isFailure())
	return sc;

    std::vector<HepVector> clusterCol;
    sc = m_fuzzyClusterTool->getClusters(m_clusterSetNo, clusterCol);
    if(sc.isFailure())
        return sc;

    if(m_fuzzyClusterTool->good())
    	nClusters = clusterCol.size();
    else
	log << MSG::INFO << "FC:  number of Cal hits = " << dataSetSize << 
		" wrong results. Applying one single cluster!" << endreq;
}
else if(dataSetSize != 0)
	log << MSG::INFO << "FC: number of Cal hits = " << dataSetSize << 
	std::endl << " is less than maximum cluster number: " << 
	maxClusters << std::endl << " Applying one single cluster!" << endreq;

    double *ene=new double[nClusters];                 // total energy in calorimeter
    Vector *pCluster=new Vector[nClusters];            // cluster position
    // energy per layer
	std::vector<double> *eneLayer=new std::vector<double>[nClusters];
    // Vector of average position per layer
	std::vector<Vector> *pLayer=new std::vector<Vector>[nClusters];
    
    // Vector of quadratic spread per layer
	std::vector<Vector> *rmsLayer=new std::vector<Vector>[nClusters];
    
    for(int iClusterNo=0; iClusterNo<nClusters; iClusterNo++) {
	    ene[iClusterNo] = 0;
    	    pCluster[iClusterNo] = p0;            // cluster position
	    eneLayer[iClusterNo] = std::vector<double>(nLayers,0.);
	    pLayer[iClusterNo] = std::vector<Vector>(nLayers);
	    rmsLayer[iClusterNo] = std::vector<Vector>(nLayers);
    }

    
    // Compute barycenter and various moments
    
    // loop over all hitted crystals in CalXtalRecCol
    int iDataPoint;
    for (it = calXtalRecCol->begin(), iDataPoint=0; 
		    it != calXtalRecCol->end(); it++, iDataPoint++){

        // get pointer to teh reconstructed data for given crystal
		Event::CalXtalRecData* recData = *it;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy(); // crystal energy
        Vector pXtal = recData->getPosition() - p0; // Vector of crystal position
        int layer = (recData->getPackedId()).getLayer(); // layer number
        
	int clustNo = 0;

	if((dataSetSize >= maxClusters) && (m_fuzzyClusterTool->good())) {

	// get the cluster number of the hit
	    sc = m_fuzzyClusterTool->getClusterNumber(m_clusterSetNo, 
			iDataPoint, clustNo);
	
            if(sc.isFailure())
                 return sc;
	}

        // update energy of corresponding layer
        eneLayer[clustNo][layer]+=eneXtal;
        
        // update average position of corresponding layer
        Vector ptmp = eneXtal*pXtal;
        pLayer[clustNo][layer] += ptmp;
        
        // Vector containing squared coordinates, weighted by crystal energy 
        Vector ptmp_sqr(ptmp.x()*pXtal.x(),
                        ptmp.y()*pXtal.y(),
                        ptmp.z()*pXtal.z());
 
        
        
        // update quadratic spread, which is proportional to eneXtal;
        // this means, that position error in one crystal
        // is assumed to be 1/sqrt(eneXtal) 
        rmsLayer[clustNo][layer] += ptmp_sqr;
        
        
        // update energy sum
        ene[clustNo]  += eneXtal;

        // update cluster position
        pCluster[clustNo] += ptmp;
    }
    
// loop over all clusters
for(int iCluster=0; iCluster<nClusters; iCluster++) {
    // Now take the means

    // if energy sum is not zero - normalize cluster position
      if(ene[iCluster]>0.)pCluster[iCluster] *= (1./ene[iCluster]); 

    	// if energy is zero - set cluster position to non-physical value
      else pCluster[iCluster]=Vector(-1000., -1000., -1000.);

      int i = 0;
    
      // loop over calorimeter layers
      for( ;i<nLayers;i++){

        // if energy in the layer is not zero - finalize calculations
        if(eneLayer[iCluster][i]>0)
        {
            // normalize position in the layer
            pLayer[iCluster][i] *= (1./eneLayer[iCluster][i]); 
            
            // normalize quadratic spread in the laye
            rmsLayer[iCluster][i] *= (1./eneLayer[iCluster][i]);
            
            
            // Vector containing the squared average position in each component
            Vector sqrLayer(pLayer[iCluster][i].x()*pLayer[iCluster][i].x(),
                            pLayer[iCluster][i].y()*pLayer[iCluster][i].y(),
                            pLayer[iCluster][i].z()*pLayer[iCluster][i].z());
            
            
            // the precision of transverse coordinate measurement
            // if there is no fluctuations: 1/sqrt(12) of crystal width
            Vector d; 
            if(i%2 == 1) d = Vector(m_CsIWidth*m_CsIWidth/12.,0.,0.);
            else d = Vector(0.,m_CsIWidth*m_CsIWidth/12.,0.);
            
            // subtracting the  squared average position and adding
            // the square of crystal width, divided by 12
            rmsLayer[iCluster][i] += d-sqrLayer;
            
        }
        
        // if energy in the layer is zero - reset position and spread Vectors
        else 
        {
            pLayer[iCluster][i]=p0;
            rmsLayer[iCluster][i]=p0;
        }
      }
    
    // Now sum the different rms to have one transverse and one longitudinal rms
    double rms_trans = 0;
    double rms_long = 0;
    std::vector<Vector> posrel(nLayers);
    
    
      for(int ilayer=0;ilayer<nLayers;ilayer++)
      {
        
        posrel[ilayer]=pLayer[iCluster][ilayer]-pCluster[iCluster];
        
        // Sum alternatively the rms
        if(ilayer%2)
        {
            rms_trans += rmsLayer[iCluster][ilayer].y();
            rms_long += rmsLayer[iCluster][ilayer].x();
        }
        else
        {
            rms_trans += rmsLayer[iCluster][ilayer].x();
            rms_long += rmsLayer[iCluster][ilayer].y();
        }
      }
 
    // Compute direction using the positions and rms per layer
    Vector caldir = Fit_Direction(pLayer[iCluster],rmsLayer[iCluster],nLayers);
 
	    // Fill CalCluster data
    Event::CalCluster* cl = new Event::CalCluster(ene[iCluster],pCluster[iCluster]+p0);

    cl->initialize(ene[iCluster],
		   eneLayer[iCluster],
		   pLayer[iCluster],
		   rmsLayer[iCluster],
		   rms_long,
		   rms_trans,
		   caldir,
		   0.);

    getClusterCol()->add(cl);

} 
	delete[] ene;
	delete[] pCluster;
	delete[] eneLayer;
	delete[] pLayer;
	delete[] rmsLayer;
	return sc;
}

StatusCode FuzzyClusterTool::finalize()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

StatusCode FuzzyClusterTool::execute()
{
	StatusCode sc = StatusCode::SUCCESS;

	return sc;
}

Vector FuzzyClusterTool::Fit_Direction(std::vector<Vector> pos,
                                     std::vector<Vector> sigma2,
                                     int nlayers)
//
// Purpose and Method:
//       find the particle direction from average positions in each
//       layer and their quadratic errors. The fit is made independently
//       in XZ and YZ planes for odd and even layers, respectively.
//       Then the 3-vector of particle direction is calculated from these 
//       2 projections.
//       Only position information based on the transversal crystal position
//       is used. The position along the crystal, calculated from signal
//       asymmetry is not used.
//
// Inputs:
//         pos      - average position for each calorimeter layer
//         sigma2   - quadratic error of position measurement for each layer
//                    in each direction  
//         nlayers  - number of calorimeter layers
//
// Returned value:    3-Vector of reconstructred particle direction
//
                                     
{
    
    MsgStream log(msgSvc(), name());
    
    // sigma2.z() is useless here no matter its value.
    double cov_xz = 0;  // covariance x,z
    double cov_yz = 0;  // covariance y,z
    double mx=0;        // mean x
    double my=0;        // mean y
    double mz1=0;       // mean z for x pos
    double mz2=0;       // mean z for y pos
    double norm1=0;     // sum of weights for odd layers
    double norm2=0;     // sum of weights for even layers
    double var_z1=0;    // variance of z for odd layers	
    double var_z2=0;    // variance of z for even layers
    
    // number of layers with non-zero energy deposition
    // in X and Y direction, respectively
    int nlx=0,nly=0;
    
    
    // "non-physical vector of direction, which is returned
    // if fit is imposible due to insufficient number of hitted layers
    Vector nodir(-1000.,-1000.,-1000.);
    
    
    // loop over calorimeter layers
    for(int il=0;il<nlayers;il++)
    {                
        // For the moment forget about longitudinal position
        
        // odd layers used for XZ fit
        if(il%2==1)
        {
            
            // only if the spread in X direction is not zero,
            // which is the case if there is non-zero energy
            // deposition in this layer
            if (sigma2[il].x()>0.)
            {
                nlx++; // counting layers used for the fit in XZ plane 
                
                // calculate weighting coefficient for this layer
                double err = 1/sigma2[il].x(); 
                
                // calculate sums for least square linear fit in XZ plane
                cov_xz += pos[il].x()*pos[il].z()*err;
                var_z1 += pos[il].z()*pos[il].z()*err;
                mx += pos[il].x()*err;
                mz1 += pos[il].z()*err;
                norm1 += err;
            }
        }
        // even layers used for YZ fit
        else
        {
            // only if the spread in Y direction is not zero,
            // which is the case if there is non-zero energy
            // deposition in this layer
            if(sigma2[il].y()>0.)
            {
                
                nly++; // counting layers used for the fit in YZ plane 
                
                // calculate weighting coefficient for this layer
                double err = 1/sigma2[il].y();
                
                
                // calculate sums for least square linear fit in YZ plane
                cov_yz += pos[il].y()*pos[il].z()*err;
                var_z2 += pos[il].z()*pos[il].z()*err;
                my += pos[il].y()*err;
                mz2 += pos[il].z()*err;
                norm2 += err;
            }
        }
    }		
    
    // linear fit requires at least 2 hitted layers in both XZ and YZ planes
    // otherwise non-physical direction is returned
    // which means that direction hasn't been found
    if(nlx <2 || nly < 2 )return nodir;
    
    
    
    
    mx /= norm1;
    mz1 /= norm1;
    cov_xz /= norm1;
    cov_xz -= mx*mz1;
    var_z1 /= norm1;
    var_z1 -= mz1*mz1;
    
    // protection against dividing by 0 in the next statment
    if(var_z1 == 0) return nodir;
    
    // Now we have cov(x,z) and var(z) we can
    // deduce slope in XZ plane
    double tgthx = cov_xz/var_z1;
    
    
    my /= norm2;
    mz2 /= norm2;
    cov_yz /= norm2;
    cov_yz -= my*mz2;
    var_z2 /= norm2;
    var_z2 -= mz2*mz2;
    
    // protection against dividing by 0 in the next statment
    if(var_z2 == 0) return nodir;
    
    // Now we have cov(y,z) and var(z) we can
    // deduce slope in YZ plane
    double tgthy = cov_yz/var_z2;
    
    // combining slope in XZ and YZ planes to get normalized 3-vector
    // of particle direction
    double tgtheta_sqr = tgthx*tgthx+tgthy*tgthy;
    double costheta = 1/sqrt(1+tgtheta_sqr);
    Vector dir(costheta*tgthx,costheta*tgthy,costheta);
    return dir;
}



