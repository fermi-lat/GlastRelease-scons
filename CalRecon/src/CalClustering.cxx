
#include "CalClustering.h"
#include "GaudiKernel/MsgStream.h"

CalClustering::CalClustering
 ( const std::string & type, 
   const std::string & name,
   const IInterface * parent )
 : AlgTool( type, name, parent ), m_calReconSvc(0)
 { declareInterface<ICalClustering>(this) ; }

CalClustering::~CalClustering()
 {} 
    
StatusCode CalClustering::initialize()
 {
  MsgStream log(msgSvc(),name()) ;
  if (service("CalReconSvc",m_calReconSvc,true).isFailure())
   {
    log<<MSG::ERROR<<"Could not find CalReconSvc"<<endreq ;
    return StatusCode::FAILURE ;
   }
  return StatusCode::SUCCESS ;
 }

StatusCode CalClustering::findClusters()
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
  MsgStream log(msgSvc(),name()) ;

  // prepare the initital set of xtals
  XtalDataVec xtals ;
  getXtals(xtals) ;
  
  // find out groups, with at least a zero
  // cluster so that downstream code runs ok
  XtalDataVecVec clusters ;
  if (xtals.size()>0)
   { makeSets(xtals,clusters) ; }
  else
   { clusters.push_back(new XtalDataVec) ; }

  // make the clusters
  setClusters(clusters) ;
  XtalDataVecVec::iterator cluster ;
  for ( cluster = clusters.begin() ;
        cluster != clusters.end() ;
        ++cluster )
    delete (*cluster) ;
  
  return StatusCode::SUCCESS ;
 }

//! Collect CalXtalRecData pointers
void CalClustering::getXtals( XtalDataVec & xtals )
 {
  xtals.clear();
  Event::CalXtalRecCol::const_iterator it ;
  for ( it = m_calReconSvc->getXtalRecs()->begin() ; it != m_calReconSvc->getXtalRecs()->end() ; ++it )
   {
    // get pointer to the reconstructed data for given crystal
	Event::CalXtalRecData * recData = *it ;
    xtals.push_back(recData) ;
   }
 }
 
/// This makes CalClusters out of associated CalXtalRecData pointers
void CalClustering::setClusters
 ( const XtalDataVecVec & clusters )
 {
  m_calReconSvc->getClusters()->delClusters() ;
  int calnLayers = m_calReconSvc->getCalNLayers() ;
  XtalDataVecVec::const_iterator cluster ;
  for ( cluster = clusters.begin() ;
        cluster != clusters.end() ;
        ++cluster )
   {    
    const Point p0(0.,0.,0.);  

    //Initialize variables
    double ene = 0;                                 // Total energy in this cluster
    Vector pCluster(0.,0.,0.);                      // Cluster position
    std::vector<double> eneLayer(calnLayers,0.);     // Energy by layer
    std::vector<Vector> pLayer(calnLayers);       // Position by layer
    std::vector<Vector> rmsLayer(calnLayers);     // rms by layer

    // Compute barycenter and various moments
    
    // loop over all crystals in the current cluster
    XtalDataVec::const_iterator xTalIter;
    for(xTalIter = (**cluster).begin(); xTalIter != (**cluster).end(); xTalIter++)
    {
        // get pointer to the reconstructed data for given crystal
		Event::CalXtalRecData* recData = *xTalIter;
        
        // get reconstructed values
        double eneXtal = recData->getEnergy();                // crystal energy
        Vector pXtal   = recData->getPosition() - p0;         // Vector of crystal position
        int    layer   = (recData->getPackedId()).getLayer(); // layer number

        // update energy of corresponding layer
        eneLayer[layer] += eneXtal;
        
        // update average position of corresponding layer
        Vector ptmp = eneXtal*pXtal;
        pLayer[layer] += ptmp;
        
        // Vector containing squared coordinates, weighted by crystal energy 
        Vector ptmp_sqr(ptmp.x()*pXtal.x(), ptmp.y()*pXtal.y(), ptmp.z()*pXtal.z());
 
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
    if(ene > 0.) pCluster /= ene; 
 	// if energy is zero - set cluster position to non-physical value
    else pCluster = Vector(-1000., -1000., -1000.);
    
    // loop over calorimeter layers
    int i ;
    for( i = 0; i < calnLayers; i++)
    {
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
            double csIWidth = m_calReconSvc->getCalCsIWidth() ;
            if(i%2 == 1) d = Vector(csIWidth*csIWidth/12.,0.,0.);
            else d = Vector(0.,csIWidth*csIWidth/12.,0.);
            
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
    double rms_trans = 0;
    double rms_long = 0;
    std::vector<Vector> posrel(calnLayers);
    
    for(int ilayer = 0; ilayer < calnLayers ; ilayer++)
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
 
    // Compute direction using the positions and rms per layer
    Vector caldir = fitDirection(pLayer, rmsLayer);

    // Fill CalCluster data
    Event::CalCluster* cl = new Event::CalCluster(ene, pCluster + p0);

    cl->initialize(ene, eneLayer, pLayer, rmsLayer, rms_long, rms_trans, caldir, 0.);

    m_calReconSvc->getClusters()->add(cl);

   }
  return ;
 }

Vector CalClustering::fitDirection
 ( std::vector<Vector> pos,
   std::vector<Vector> sigma2 )
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
    for(int il=0;il<m_calReconSvc->getCalNLayers();il++)
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
	
