
#include "ProfileTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

static const ToolFactory<ProfileTool>  s_factory;
const IToolFactory& ProfileToolFactory = s_factory;

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
    
    double length = ((m_xtalHeight*i+par[2])/1.85)/getStaticSlope();
    
    // Evaluation of the parameters of CsI at this energy	
    //	double alpha = par[1];
    //	double lambda = par[3];
    
    double alpha = 2.65*exp(0.15*log(par[0]));
    double lambda = 2.29*exp(-0.031*log(par[0]));
    
    
    double x=length/lambda;
    double dx = m_xtalHeight / (1.85 *lambda)/getStaticSlope();
    
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


ProfileTool::ProfileTool( const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : AlgTool(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<IEnergyCorr>(this);
    // Declare the properties that may be set in the job options file
    
}



StatusCode ProfileTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing ProfileTool" <<endreq;
    
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
    

    double CsIWidth;
    if(!detSvc->getNumericConstByName(std::string("CsIWidth"),&CsIWidth))
    {
        log << MSG::ERROR << " constant " << " CsIWidth "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    double CsIHeight;
    if(!detSvc->getNumericConstByName(std::string("CsIHeight"),&CsIHeight))
    {
        log << MSG::ERROR << " constant " << " CsIHeight "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 

    m_xtalHeight = CsIHeight/10.;  // crystal height in cm
    m_xtalWidth = CsIWidth/10.;    // crystal width in cm
    
    
    // Minuit object
    m_minuit = new TMinuit(5);
    
    //Sets the function to be minimized
    m_minuit->SetFCN(fcn);
    
    
    m_g_elayer.clear();

    

    return sc;
}


StatusCode ProfileTool::doEnergyCorr(double eTotal, Event::CalCluster* cluster)

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

    // defines global variable to be used for fcn
    m_g_elayer.resize(getNLayers());
    for (int i =0;i<getNLayers();i++)
    {
        // We are working in GeV
        m_g_elayer[i] = eTotal/1000.;
    }
    m_nbins = getNLayers();
    
    if( eTotal<2000. || getTrackSlope() == 0) //algorithm is useless under several GeV
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
        
        // Get chi-square
        double edm,errdef;
        int nvpar,nparx,icstat;
        m_minuit->mnstat(ki2,edm,errdef,nvpar,nparx,icstat);
        
        // Fills data
        
        cluster->initProfile(1000*fit_energy,ki2,fit_start,fit_alpha,fit_lambda);
        
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

StatusCode ProfileTool::execute()
{
    StatusCode sc = StatusCode::SUCCESS;
    
    return sc;
}

     double ProfileTool::m_xtalHeight=0.; //!< xtal height in cm
     double ProfileTool::m_xtalWidth=0.;  //!< xtal width  in cm
     int ProfileTool::m_nbins=0;  //!< Number of bins used for the fit
     std::vector<double> ProfileTool::m_g_elayer;  //!< Energy per layer in GeV
