
#include "CalClustersAlg.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"


/// Glast specific includes
#include "Event/TopLevel/EventModel.h"
#include "GaudiKernel/ObjectVector.h"

//Gamma function and Minuit
#include "TMath.h"
#include "TMinuit.h"

int nbins;  //!< Number of bins used for the fit
std::vector<double> g_elayer;  //!< Energy per layer in GeV
double slope;   //!< slope of the shower direction
double xtalHeight; //!< xtal height in cm
double xtalWidth;  //!< xtal width  in cm


static const AlgFactory<CalClustersAlg>  Factory;
const IAlgFactory& CalClustersAlgFactory = Factory;

using namespace Event;

//! function to compute the true energy deposited in a layer
/*! Uses the incomplete gamma function:
 *       gamma(double,double) implemented in gamma.cxx
*/ 
static double gam_prof(double *par, int i)

// Purpose and method:
//        To find the integral of the shower profile over the layer i.
//
//        Takes into account the shower angle with vertical direction. 
//
//  Inputs:
//           par - pointer to the array of 4 shower profile parameters,
//                 only par[0] and par[2] are used, while 2 other parameters
//                 are evaluated using predefined formulas;
//           i   - layer number
//
//   returned: integral of shower profile over requested layer.
//
{
	double result =0; 

	double length = ((xtalHeight*i+par[2])/1.85)/slope;

	// Evaluation of the parameters of CsI at this energy	
	//	double alpha = par[1];
	//	double lambda = par[3];

	double alpha = 2.65*exp(0.15*log(par[0]));
        double lambda = 2.29*exp(-0.031*log(par[0]));


	double x=length/lambda;
	double dx = xtalHeight / (1.85 *lambda)/slope;
	
	double gamma1 =0;
	double gamma2 = 0;

	// Now we will calculate the gamma incomplete function

	// gamma1 = integration from 0 to x	
	gamma1 = TMath::Gamma(alpha,x);
	x += dx;

	// gamma2 = integration from 0 to x+dx	
	gamma2 = TMath::Gamma(alpha,x);
	
	// the result of integration over Xtal pathlength is E*(gamma2-gamma1)
	result = par[0]*(gamma2 - gamma1);

        return result;
}


//! the fcn function needed by Minuit
/*! Computes the chisquare ie:
\f$ \chi^2= \sum_{i=1}^{8}\frac{(\bar{E}-E_i)^2}{\sigma_i}\f$
*/

static void fcn(int &npar, double *gin, double &f, double *par, int iflag)

// Purpose: calculates the weighted sum of quadratic deviations
//          of energy depositioins in layers from predicted by shower
//          profile function
//          this function is called by Minuit package, which defines
//          the set of parameters. Not all parameters are used by this
//          function
//
// Input: par - array of parameters of fitted profile function
//         
// Output: f - the result of chi2 calculation
{
        int i;
	//calculate 'chisquare'
        double chisq = 0;
        double dlt;
   
        for (i=0;i<nbins; i++) {
                dlt  = (g_elayer[i]-gam_prof(par,i));
                if(g_elayer[i]>0.001) chisq += dlt*dlt/g_elayer[i];
	}
   
        f = chisq;
}

double CalClustersAlg::Leak(double eTotal,double elast)

// Purpose and method:
//              
//              Calculates the energy leakage from calorimeter
//              using the fitted correlation with last layer energy
//               deposition
//  
//    Inputs:
//             eTotal - total energy deposition in the calorimeter
//             elast  - energy deposition in the last layer
//
//    Returned: energy leakage

{
    if(eTotal<200.) return 0.;
    else
    {
        // Evaluation of energy using correlation method
        // Coefficients fitted using GlastSim.
        double p0 = -1.49 + 1.72*slope;
        double p1 = 0.28 + 0.434 * slope;
        double p2 = -15.16 + 11.55 * slope;
        double p3 = 13.88 - 10.18 * slope;
        double lnE = log(eTotal/1000.);
        double funcoef = (p0 + p1 * lnE )/(1+exp(-p2*(lnE - p3)));
        
        double e_leak = funcoef * elast;
        
        // Evaluation of energy using correlation with last layer
       	// coefficients fitted using tbsim and valid for ~1GeV<E<~50GeV
       	//double slope = 1.111 + 0.557*log(eTotal/1000.);
        //double intercept = 210. + 112.* log(eTotal/1000.)*log(eTotal/1000.); 
       	//double e_leak = slope * elast + intercept;
        return e_leak;
    }
}

void CalClustersAlg::Profile(double eTotal, Event::CalCluster* cl)

// Purpose and method:
// 
//               This function fits the parameters of shower profile using
//               the Minuit minimization package and stores the fitted
//               parameters in the CalCluster object
//
//  Inputs:    eTotal - total energy deposition in the calorimeter
//             cl - pointer to the CalCluster object to store the fitted
//                  parameters 
//
//  Output:   parameters fit_energy,ki2,fit_start,fit_alpha,fit_lambda,
//            stored in CalCluster object using initProfile() method

{
    
    if( eTotal<2000. || slope == 0) //algorithm is useless under several GeV
    {
        cl->initProfile(0,0,0,0,0);
    }
    else
    {
        double fit_energy;
        double ki2;
        double fit_alpha;
        double fit_lambda;
        double fit_start;
        double energy_err;
        double alpha_err;
        double lambda_err;
        double start_err;
        
        // Vector of step, initial min and max value
        float vstrt[5];
        float stp[5];
        float bmin[5];
        float bmax[5];
        
        // We are working in GeV
        double eTotal_GeV = eTotal / 1000;
        
        // par[0] is energy
        // par[1] is alpha
        // par[2] is the starting point
        // par[3] is lambda
        
        // Init start values and bounds
        
        // starting value of each parameter
        vstrt[0] = eTotal_GeV;

        
        // parametrisation of alpha
        vstrt[1] = 2.65 * exp(0.15*log(eTotal_GeV));

        // eq to 1X0 in CsI
        vstrt[2] = 1.8f;			
        vstrt[3] = 2.29 * exp(-0.031*log(eTotal_GeV));

        // step value of each parameter
        stp[0] = 0.1f;
        stp[1] = 0.1f;
        stp[2] = 0.02f;
        stp[3] = 0.2f;
        
        // minimum value of each parameter
        bmin[0] = 0.5*eTotal_GeV;
        bmin[1] = 0.5;
        bmin[2] = -5;
        bmin[3] = 0.;
        
        // maximum value of each parameter
        bmax[0] = 5 * eTotal_GeV;
        bmax[1] = 15;
        bmax[2] = 10;
        bmax[3] = 10;
        
        // Those are the arguments for Minuit
        double arglist[10];
        int ierflg = 0;
        
        // Set no output because of tuple output
        arglist[0] = -1;
        minuit->mnexcm("SET PRI", arglist ,1,ierflg);
        
        // idem with warnings ( minuit prints warnings
        // when the Hessian matrix is not positive )
        minuit->mnexcm("SET NOW", arglist ,1,ierflg);
        
        arglist[0] = 1;
        minuit->mnexcm("SET ERR", arglist ,1,ierflg);
        
        // defines the strategy used by minuit
        // 1 is standard
        arglist[0] = 2;
        minuit->mnexcm("SET STR", arglist ,1,ierflg);
        
        
        // Defines parameters
        
        minuit->mnparm(0, "Energy",    vstrt[0], stp[0],
                       bmin[0],bmax[0],ierflg);
        minuit->mnparm(1, "Alpha",    vstrt[1], stp[1],
                       bmin[1],bmax[1],ierflg);
        minuit->mnparm(2, "Starting point",  vstrt[2], stp[2],
                       bmin[2],bmax[2],ierflg);
        minuit->mnparm(3, "Lambda",  vstrt[3], stp[3],
                       bmin[3],bmax[3],ierflg);
        
        // Fix some parameters
        // alpha
        arglist[0] = 2;
        minuit->mnexcm("FIX ", arglist ,1,ierflg);
        
        // lambda
        arglist[0] = 4;
        minuit->mnexcm("FIX ", arglist ,1,ierflg);
        
        
        // Calls Migrad with 500 iterations maximum
        arglist[0]=500;
        arglist[1]=1;
        minuit->mnexcm("MIGRAD",arglist,2,ierflg);
        
        // Retrieve parameter information 
        
        minuit->GetParameter( 0, fit_energy,energy_err ); 
        minuit->GetParameter( 1, fit_alpha,alpha_err ); 
        minuit->GetParameter( 2, fit_start,start_err ); 
        minuit->GetParameter( 3, fit_lambda,lambda_err ); 
        
        // Get chi-square
        double edm,errdef;
        int nvpar,nparx,icstat;
        minuit->mnstat(ki2,edm,errdef,nvpar,nparx,icstat);
        
        // Fills data
        
        cl->initProfile(1000*fit_energy,ki2,fit_start,fit_alpha,fit_lambda);

        // Clear minuit
        arglist[0] = 1;
        minuit->mnexcm("CLEAR", arglist ,1,ierflg);
        
    } 
        
}



Vector CalClustersAlg::Fit_Direction(std::vector<Vector> pos,
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


CalClustersAlg::CalClustersAlg(const std::string& name,
                               ISvcLocator* pSvcLocator):
Algorithm(name, pSvcLocator)
{

    // declaration of parameter needed to distinguish 2 calls
    // of CalClustersAlg:
    // 1st - before TkrRecon and 2nd - after TkrRecon
    
    declareProperty("callNumber",m_callNumber=0);
    
}



StatusCode CalClustersAlg::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc
//    - gets callNumber parameter from jobOptions file
//    - creates minimizer object (Midnight)
//    - sets pointer to the global function to be minimized
//      in Profile() function
//    - clears the global vector( g_elayer) with energies per layer in GeV

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    
    
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
    } else m_CalnLayers = value;
    
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
    
    
    
    
    // get callNumber parameter from jobOptions file    
    setProperties();
    log << MSG::INFO << "CalClustersAlg: callNumber = " 
        << m_callNumber << endreq;
    
    // set global constants with geometry parameters (in cm)
    // used by Profile() function
    
    xtalHeight = m_CsIHeight/10.;  // crystal height in cm
    xtalWidth = m_CsIWidth/10.;    // crystal width in cm
    
    
    // Minuit object
    minuit = new TMinuit(5);
    
    //Sets the function to be minimized
    minuit->SetFCN(fcn);
    
    
    g_elayer.clear();
    
    return sc;
}

StatusCode CalClustersAlg::retrieve()

// Purpose:
//    - to get pointer to existing CalXtalRecCol
//    - to create new CalClusterCol and register it in TDS 
//
//  TDS input:   CalXtalRecCol
//  TDS ourput  CalClusterCol

{
    
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    DataObject* pnode=0;
    

    // get pointer to Calrecon directory in TDS
    sc = eventSvc()->retrieveObject(EventModel::CalRecon::Event, pnode );
    
    // if this directory doesn't yet exist - attempt to create it
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject(EventModel::CalRecon::Event,
                                        new DataObject);

        
        // if can't create - put error message and return
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory"
                << endreq;
            return sc;
        }
    }

    // attempt to get pointer to CalClusterCol, if it exists already
    m_calClusterCol = SmartDataPtr<CalClusterCol> (eventSvc(),
                      EventModel::CalRecon::CalClusterCol);

    
    // if it doesn't exist - create it and register in TDS
    if (!m_calClusterCol )
    {
        m_calClusterCol = 0;
        m_calClusterCol = new CalClusterCol();
        sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,
                                        m_calClusterCol);
    }
    
    // if it exists - delete all clusters
    else
    {
        m_calClusterCol->delClusters();
    }
    
    // get pointer to CalXtalRecCol
    m_calXtalRecCol = SmartDataPtr<CalXtalRecCol>(eventSvc(),
                      EventModel::CalRecon::CalXtalRecCol); 
    
    
    return sc;
}

StatusCode CalClustersAlg::execute()

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
    sc = retrieve();

    const Point p0(0.,0.,0.);
    
    // variable indicating ( if >0) the presence of tracker
    // reconstruction output
    int rectkr=0;  

    int ntracks;
    Vector trackDirection;
    Point trackVertex;
    

    // get pointer to the tracker vertex collection
    SmartDataPtr<TkrVertexCol> tkrRecData(eventSvc(),
                 EventModel::TkrRecon::TkrVertexCol);

    // if reconstructed tracker data doesn't exist - put the debugging message
    if (tkrRecData == 0) {
        log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
        // return sc;
    }

    
    // if data exist and number of tracks not zero 
    // - get information of track position and direction 
    else
    {
        // First get reconstructed direction from tracker
        ntracks = tkrRecData->size();
        log << MSG::INFO << "number of tracks = " << ntracks << endreq;
        
        
        if (ntracks > 0) {
            rectkr++;
            trackDirection = tkrRecData->front()->getDirection();
            trackVertex = tkrRecData->front()->getPosition();
            slope = fabs(trackDirection.z());
            log << MSG::DEBUG << "track direction = " << slope << endreq;
            
        } else {
            log << MSG::INFO << "No reconstructed tracks " << endreq;
        }	
    }

    
    double ene = 0;                 // total energy in calorimeter
    Vector pCluster(p0);            // cluster position
    int nLayers = m_CalnLayers;     // number of layers
    
    // energy per layer
    std::vector<double> eneLayer(nLayers,0.);   
    
    // Vector of average position per layer
    std::vector<Vector> pLayer(nLayers);    
    
    
    // Vector of quadratic spread per layer
    std::vector<Vector> rmsLayer(nLayers);  
    
    
    
    // Compute barycenter and various moments
    
    
    // loop over all hitted crystals in CalXtalRecCol
    for (CalXtalRecCol::const_iterator it = m_calXtalRecCol->begin();
    it != m_calXtalRecCol->end(); it++){

        // get pointer to teh reconstructed data for given crystal
        CalXtalRecData* recData = *it;
        
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
    
    
    // Compute direction using the positions and rms per layer
    Vector caldir = Fit_Direction(posrel,rmsLayer,nLayers);
    
    
    // if no tracker rec then fill slope
    if(!rectkr) slope = caldir.z();
    
    
    // Fill CalCluster data
    CalCluster* cl = new CalCluster(ene,pCluster+p0);
    m_calClusterCol->add(cl);
    
    // Leakage correction
    double eleak = Leak(ene,eneLayer[nLayers-1])+ene;
    
    // iteration
    eleak = Leak(eleak,eneLayer[nLayers-1])+ene;	
    
    
    
    // defines global variable to be used for fcn
    g_elayer.resize(nLayers);
    for ( i =0;i<nLayers;i++)
    {
        // We are working in GeV
        g_elayer[i] = eneLayer[i]/1000.;
    }
    nbins = nLayers;
    
    // Do profile fitting
    Profile(ene,cl);

    // calculating the transverse offset of average position in the calorimeter
    // with respect to the position predicted from tracker information
    double calTransvOffset = 0.;
    if(ntracks>0){
        Vector calOffset = (p0+pCluster) - trackVertex;
        double calLongOffset = trackDirection*calOffset;
        calTransvOffset =sqrt(calOffset.mag2() - calLongOffset*calLongOffset);
        
    }

    // store the calculated quantities in CalCluster object
    cl->initialize(eleak,eneLayer,pLayer,rmsLayer,rms_long,rms_trans,
        caldir,calTransvOffset);
    
    
    // print the reconstruction results for debugging
    m_calClusterCol->writeOut(log << MSG::DEBUG);
    
    return sc;
}

StatusCode CalClustersAlg::finalize()
{
	StatusCode sc = StatusCode::SUCCESS;
	
	
	// delete Minuit object
        delete minuit;

	return sc;
}




