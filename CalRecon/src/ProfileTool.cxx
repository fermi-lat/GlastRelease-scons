
#include "ProfileTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
// to access an XML containing Profile Bias parameters file
#include "xmlBase/IFile.h"

DECLARE_TOOL_FACTORY(ProfileTool) ;

double ProfileTool::m_static_slope=0;

double ProfileTool::m_xtalHeight=0.; //!< xtal height in cm
double ProfileTool::m_xtalWidth=0.;  //!< xtal width  in cm
int ProfileTool::m_nbins=0;  //!< Number of bins used for the fit
std::vector<double> ProfileTool::m_g_elayer;  //!< Energy per layer in GeV
     
//! function to compute the true energy deposited in a layer
/*! Uses the incomplete gamma function:
*       gamma(double,double) implemented in gamma.cxx
*/ 

double ProfileTool::gam_prof(double *par, int i)

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
    
    double length = ((m_xtalHeight*i+par[2])/1.85)/m_static_slope ;
    
    // Evaluation of the parameters of CsI at this energy	
    //	double alpha = par[1];
    //	double lambda = par[3];
    
    double alpha = 2.65*exp(0.15*log(par[0]));
    double lambda = 2.29*exp(-0.031*log(par[0]));
    
    
    double x=length/lambda;
    double dx = m_xtalHeight / (1.85 *lambda)/m_static_slope;
    
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

void ProfileTool::fcn(int & , double *, double &f, double *par, int )
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
    
    for (i=0;i<m_nbins; i++) {
        dlt  = (m_g_elayer[i]-gam_prof(par,i));
        if(m_g_elayer[i]>0.001) chisq += dlt*dlt/m_g_elayer[i];
    }
    
    f = chisq;
}

double ProfileTool::bias( double energy )
// Purpose: calculates the mean bias for the profile method at a given energy 
// and angle. The bias is estimated as the output energy of the fit times a 
// polynomial P( log(fit_energy), cos(incidence) ). It is valid for 
// energies > 2000 MeV and cos( incidence )> 0.4. the fit energy is corrected 
// as follows: 
//      debiased energy= fit energy
//      debiased energy= fit energy + bias( debiased energy ) : repeated thrice 
//
// Input: fit energy - energy in GeV
// Output:  - bias
{
  // check for valid angles
  if( m_static_slope < m_BiasCTLim ) return energy;

  // get value: cos(angle of incidence)
  double cTh= 1./sqrt(1+m_static_slope*m_static_slope);
  double debiasedE= energy;

  // correction repeated thrice
  for( short int times= 0; times<3; ++times ){ 
    double logE= log( debiasedE ); 
    double bias= 0.;

    // evaluate polynomial at logE and cTh
    for( short int biasPar=0; biasPar<m_BiasParNb; biasPar+= m_BiasCTNb ){
      double thetaPar=0.;
      for( short int ctPar= biasPar; ctPar<biasPar+m_BiasCTNb; ++ctPar ){
        thetaPar+= m_BiasVector[ctPar];
        thetaPar*= cTh;
      }
      bias+= thetaPar; 
      bias*= logE;
    }

    // apply correction
    debiasedE= energy+bias*debiasedE;
  }

  return debiasedE;  
}

ProfileTool::ProfileTool( const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : EnergyCorr(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<IEnergyCorr>(this);
    declareProperty ("xmlFile", m_xmlFile="$(CALRECONROOT)/xml/CalProfile.xml");
    // Declare the properties that may be set in the job options file
}



StatusCode ProfileTool::initialize()
// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc
{
    if (EnergyCorr::initialize().isFailure())
     { return StatusCode::FAILURE ; }

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing ProfileTool" <<endreq;
        
    // Minuit object
    m_minuit = new TMinuit(5);
    
    //Sets the function to be minimized
    m_minuit->SetFCN(fcn);
        
    m_g_elayer.clear();
   
    // Read in the parameters from the XML file
    xmlBase::IFile m_ifile(m_xmlFile.c_str());
    if ( m_ifile.contains("profileBias", "CosThetaLimit" ) ){
      m_BiasCTLim= m_ifile.getDouble("profileBias", "CosThetaLimit" );
      log << MSG::INFO << " value for CosThetaLimit " 
          << m_BiasCTLim << endreq;
    } else return StatusCode::FAILURE;

    if ( m_ifile.contains("profileBias", "BiasParNb" ) ){
      m_BiasParNb= m_ifile.getInt("profileBias", "BiasParNb" );
      log << MSG::INFO << " value for BiasParNb" 
          << m_BiasParNb << endreq;
    } else return StatusCode::FAILURE;

    if ( m_ifile.contains("profileBias", "CosThetaParNb" ) ){
      m_BiasCTNb= m_ifile.getInt("profileBias", "CosThetaParNb" );
      log << MSG::INFO << " value for CosThetaParNb " 
          << m_BiasCTNb << endreq;
    } else return StatusCode::FAILURE;

    m_BiasVector= std::vector<double> ( m_BiasParNb );

    char name[10];
    for( short int biasPar=0; biasPar<m_BiasParNb; biasPar+= m_BiasCTNb )
      for( short int ctPar= biasPar; ctPar<biasPar+m_BiasCTNb; ++ctPar ){
        sprintf( name, "bias%d", ctPar );
        if ( m_ifile.contains("profileBias", name) ){
          m_BiasVector[ctPar] = m_ifile.getDouble("profileBias", name);
	        log << MSG::INFO << " value for " << name << " " 
              << m_BiasVector[ctPar] << endreq;
        } else return StatusCode::FAILURE;
      }

    return sc;
}


StatusCode ProfileTool::doEnergyCorr( Event::CalCluster * cluster )
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
    MsgStream lm(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    double eTotal = cluster->getEnergySum() ;
    
    m_static_slope = getKernel()->getSlope(cluster) ;

    m_xtalHeight = getKernel()->getCalCsIHeight()/10.;  // crystal height in cm
    m_xtalWidth = getKernel()->getCalCsIWidth()/10.;    // crystal width in cm

    // defines global variable to be used for fcn
    m_nbins = getKernel()->getCalNLayers() ;
    m_g_elayer.resize(m_nbins);
    for (int i =0;i<m_nbins;i++)
    {
        // We are working in GeV
        m_g_elayer[i] = cluster->getEneLayer(i)/1000.;
    }
    
    if( eTotal<2000. || m_static_slope == 0) //algorithm is useless under several GeV
    {
        cluster->initProfile(0,0,0,0,0);
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
        m_minuit->mnexcm("SET PRI", arglist ,1,ierflg);
        
        // idem with warnings ( minuit prints warnings
        // when the Hessian matrix is not positive )
        m_minuit->mnexcm("SET NOW", arglist ,1,ierflg);
        
        arglist[0] = 1;
        m_minuit->mnexcm("SET ERR", arglist ,1,ierflg);
        
        // defines the strategy used by minuit
        // 1 is standard
        arglist[0] = 2;
        m_minuit->mnexcm("SET STR", arglist ,1,ierflg);
        
        
        // Defines parameters
        
        m_minuit->mnparm(0, "Energy",    vstrt[0], stp[0],
            bmin[0],bmax[0],ierflg);
        m_minuit->mnparm(1, "Alpha",    vstrt[1], stp[1],
            bmin[1],bmax[1],ierflg);
        m_minuit->mnparm(2, "Starting point",  vstrt[2], stp[2],
            bmin[2],bmax[2],ierflg);
        m_minuit->mnparm(3, "Lambda",  vstrt[3], stp[3],
            bmin[3],bmax[3],ierflg);
        
        // Fix some parameters
        // alpha
        arglist[0] = 2;
        m_minuit->mnexcm("FIX ", arglist ,1,ierflg);
        
        // lambda
        arglist[0] = 4;
        m_minuit->mnexcm("FIX ", arglist ,1,ierflg);
        
        
        // Calls Migrad with 500 iterations maximum
        arglist[0]=500;
        arglist[1]=1;
        m_minuit->mnexcm("MIGRAD",arglist,2,ierflg);
        
        // Retrieve parameter information 
        
        m_minuit->GetParameter( 0, fit_energy,energy_err ); 
        m_minuit->GetParameter( 1, fit_alpha,alpha_err ); 
        m_minuit->GetParameter( 2, fit_start,start_err ); 
        m_minuit->GetParameter( 3, fit_lambda,lambda_err ); 
        
        // bias correction
        double fit_energy_opt= bias(fit_energy);
        
        // Get chi-square
        double edm,errdef;
        int nvpar,nparx,icstat;
        m_minuit->mnstat(ki2,edm,errdef,nvpar,nparx,icstat);

        // Fills data
        cluster->initProfile(1000*fit_energy_opt,ki2,
                             fit_start,fit_alpha,fit_lambda);
        
        // Clear minuit
        arglist[0] = 1;
        m_minuit->mnexcm("CLEAR", arglist ,1,ierflg);
    } 
    
    return sc;
}

StatusCode ProfileTool::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;
 
    // delete Minuit object
    delete m_minuit;

    return sc;
}

