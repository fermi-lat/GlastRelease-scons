
#include "CalRecon/CalClustersAlg.h"
// #include "Event/dataManager.h"
// #include "Event/messageManager.h"
#include "CalRecon/CsIClusters.h"
#include "CalRecon/CalRecLogs.h"
#include "CalRecon/gamma.h"
#include "CalRecon/Midnight.h"
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

static const AlgFactory<CalClustersAlg>  Factory;
const IAlgFactory& CalClustersFactory = Factory;

int nbins;  //!< Number of bins used for the fit
std::vector<double> g_elayer;  //!< Energy per layer in GeV
double slope;   //!< slope of the shower direction

//! function to compute the true energy deposited in a layer
/*! Uses the incomplete gamma function: gamma(double,double) implemented in gamma.cxx
*/ 
static double gam_prof(double *par, int i)
{
	double result =0; 

	double length = ((2.3*i+par[2])/1.85)/slope;

	// Evaluation of the parameters of CsI at this energy	

	//	double alpha = par[1];

	//	double lambda = par[3];

	double alpha = 2.65*exp(0.15*log(par[0]));
        double lambda = 2.29*exp(-0.031*log(par[0]));


	double x=length/lambda;
	double dx = 2.3 / (1.85 *lambda)/slope;
	
	double gamma1 =0;
	double gamma2 = 0;

	// Now we will calculate the gamma incomplete function

	// gamma1 = integration from 0 to x	
	gamma1 = Gamma(alpha,x);
	x += dx;

	// gamma2 = integration from 0 to x+dx	
	gamma2 = Gamma(alpha,x);
	
	// the result of integration over Xtal pathlength is E*(gamma2 - gamma1)
	result = par[0]*(gamma2 - gamma1);
//	std::cout<< result << "\n";
	return result;
}


//! the fcn function needed by Minuit
/*! Computes the chisquare ie:
\f$ \chi^2= \sum_{i=1}^{8}\frac{(\bar{E}-E_i)^2}{\sigma_i}\f$
*/
static void fcn(int &npar, double *gin, double &f, double *par, int iflag)
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

//! Leakage correction using last layer correlation
/*!The second method uses the correlation between the escaping energy and the energy
*  deposited in the last layer of the calorimeter.  Indeed, the last layer carries
* the most important information concerning the leaking energy:  the total number
* of particles escaping through the back should be nearly proportional to the 
* energy deposited in the last layer.  The measured signal in that layer can therefore
* be modified to account for the leaking energy. 
*
* We used the Monte Carlo simulation of the GLAST beam test configuration to deter
* mine this correlation at several energies, from 2 GeV up to 40 GeV.  For one par
* ticular incident energy, the bidimensionnal distribution of the energy escaping 
* and the energy deposited in the last layer can be fitted by a simple linear function:
* \f[ E_{leak} = \alpha E_{last} + \beta \f]
* The \f$\alpha\f$ and \f$\beta\f$ parameters are proportional to the logarithm of the 
* incident energy and to its square, respectively.  Because the only information we 
* have, initially, on the incident energy is the total measured energy \f$E_m\f$, we have
* to use it as the estimator of \f$E_0\f$.  The reconstructed energy is then:
* \f[ E_{rec} = E_m + \alpha(E_m) E_{last} + \beta(E_m) \f]
* To improve the result, one can iterate using the new estimator to determine the 
* correct values of \f$\alpha\f$ and \f$\beta\f$.  
*
*\par The methods takes 2 arguments:
*\param eTotal Total energy measured in the calorimeter in MeV
*\param elast  Total energy measured in the last layer in MeV
*
*\return Corrected energy in MeV
*
* \warning Has been developped on axis only
* \warning Give biased results when the last layergains are misaligned
*
* \author Regis Terrier
*
* \b Revision: 
* - 10/17/00    RT    comments added
* - 09/15/00    RT    bias correction added, tuned on MC data
* - 08/20/00    RT    first implementation
*/

//################################################
double CalClustersAlg::Leak(double eTotal,double elast)
//################################################
{
	if(eTotal<200.) return 0.;
	else
	{
	        // Evaluation of energy using correlation woth last layer
        	// coefficients fitted using tbsim and valid for ~1GeV<E<~50GeV
        	double slope = 1.111 + 0.557*log(eTotal/1000.);
		double intercept = 210. + 112. * log(eTotal/1000.) *log(eTotal/1000.); 
        	double e_leak = slope * elast + intercept;
		return e_leak;
	}
}

//! Longitudinal profile fitting method
/*! It performs a longitudinal profile fitting using the incomplete gamma function
* ( gamma.cxx). The mean energy density per length unit is taken as:
* \f[
f(z) = \frac{1}{\lambda} \frac{exp(-z/\lambda)(z/\lambda)^{(1-\alpha)}}{\Gamma(\alpha)}
\f]
* Thus integrated on the i th Xtal pathlength it gives 
* \f[
E_i = E_{tot}(\Gamma_{inc}(z_i/\lambda,\alpha)-\Gamma_{inc}(z_{i+1}/\lambda,\alpha))
\f]
* the 2 shower parameters here alpha and lambda describes the maximum position and 
* and the exponential decrease of the profile. Those parameters have been estimated using
* a MC and their dependance over E has been fitted by a power law. They are log-normally 
* distributed with a very broad distribution ( hence some of the shower fluctuations ).
* Therefore they should not be included in the fitting process. 
*
* Here we use 4 parameters: total energy
*                           starting point
*                           alpha
*                           lambda
* They can be fixed or released in Profile()
*
*
* The input is:
* \param eTotal total energy measured in the calorimeter in MeV
* \param cl     the CsICluster in which the results are saved
*
* The output ( the 4 fitting parameters and the chi square) of this method is
* stored in the CsICluster.  
* 
* \warning The algorithm works only on axis thats why slope=1 in execute()
* \warning Gives sensible results only at large enough energies ie ~10GeV
*
* \b Revision:
* - 10/17/00    RT    comments added
* - 05/00       RT    first implementation
*/
//################################################
void CalClustersAlg::Profile(double eTotal, CsICluster* cl)
//################################################
{

  if( eTotal<2000. || slope == 0) //algorithm is useless under several GeV

	{
 	 cl->setFitEnergy(0);    // Conversion to MeV
 	 cl->setProfChisq(0);
  	 cl->setCsiStart(0);
  	 cl->setCsiAlpha(0);
  	 cl->setCsiLambda(0);
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
  vstrt[1] = 2.65 * exp(0.15*log(eTotal_GeV));  // parametrisation of alpha
  vstrt[2] = 1.8f;			// eq to 1X0 in CsI
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
    
  // idem with warnings ( minuit prints warnings when the Hessian matrix is not positive )
  minuit->mnexcm("SET NOW", arglist ,1,ierflg);
  
  arglist[0] = 1;
  minuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  // defines the strategy used by minuit
  // 1 is standard
  arglist[0] = 2;
  minuit->mnexcm("SET STR", arglist ,1,ierflg);
  
  
  // Defines parameters
  
  minuit->mnparm(0, "Energy",    vstrt[0], stp[0], bmin[0],bmax[0],ierflg);
  minuit->mnparm(1, "Alpha",    vstrt[1], stp[1], bmin[1],bmax[1],ierflg);
  minuit->mnparm(2, "Starting point",  vstrt[2], stp[2], bmin[2],bmax[2],ierflg);
  minuit->mnparm(3, "Lambda",  vstrt[3], stp[3], bmin[3],bmax[3],ierflg);
  
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
  cl->setFitEnergy(1000*fit_energy);    // Conversion to MeV
  cl->setProfChisq(ki2);
  cl->setCsiStart(fit_start);
  cl->setCsiAlpha(fit_alpha);
  cl->setCsiLambda(fit_lambda); 

  
  // Clear minuit
  arglist[0] = 1;
  //minuit->mnexcm("CLEAR", arglist ,1,ierflg);

 } 

}


//################################################
CalClustersAlg::CalClustersAlg(const std::string& name, ISvcLocator* pSvcLocator):
//################################################
Algorithm(name, pSvcLocator)
{


}

//################################################
StatusCode CalClustersAlg::initialize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
     sc = service("CalGeometrySvc", m_CalGeo);
	
//	m_CsIClusterList = dataManager::instance()->getData("CsIClusterList",m_CsIClusterList);
//	m_CalRecLogs  = dataManager::instance()->getData("CalRecLogs",m_CalRecLogs);
//	m_CsIClusterList->clear();

        // Minuit object
        minuit = new Midnight(5);

        //Sets the function to be minimized
        minuit->SetFCN(fcn);

	
	g_elayer.clear();
   
	return sc;
}

//################################################
StatusCode CalClustersAlg::retrieve()
//################################################
{
	
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	m_CsIClusterList = 0;
	m_CsIClusterList = new CsIClusterList();

	DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( "/Event/CalRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/CalRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }

    

	m_CalRecLogs = SmartDataPtr<CalRecLogs>(eventSvc(),"/Event/CalRecon/CalRecLogs"); 


//	sc = eventSvc()->retrieveObject("/Event/CalRecon/CalADCLogs",m_CalRawLogs);
	 sc = eventSvc()->registerObject("/Event/CalRecon/CsIClusterList",m_CsIClusterList);
	return sc;
}

/*!Performs the reconstruction. 
 * - Energy per layer is computed and stored in CsICluster in MeV
 * - Barycenter per layer is also computed and stored in CsICluster
 * - Performs high energy corrections: see Profile() and Leak() for details
 */
//################################################
StatusCode CalClustersAlg::execute()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
	sc = retrieve();

	int nLogs = m_CalRecLogs->num();
	double ene = 0;
	const Point p0(0.,0.,0.);
	Vector pCluster(p0);
	int nLayers = m_CalGeo->numLayers() * m_CalGeo->numViews();
	
	std::vector<double> eneLayer(nLayers,0.);
	std::vector<Vector> pLayer(nLayers);
		
	
	
	for (int jlog = 0; jlog < nLogs ; jlog++) {
		CalRecLog* recLog = m_CalRecLogs->Log(jlog);

		double eneLog = recLog->energy();
		Vector pLog = recLog->position() - p0;
		int layer = nLayers-1 - (recLog->layer() * 2 + recLog->view());
		eneLayer[layer]+=eneLog;
		pLayer[layer] += eneLog*pLog;
		ene  += eneLog;
		pCluster += eneLog*pLog;				
	}

	pCluster *= (1./ene);
	int i = 0;
	for( ;i<nLayers;i++){
			if(eneLayer[i]>0)pLayer[i] *= (1./eneLayer[i]);
			else pLayer[i]=p0;
	}


	CsICluster* cl = new CsICluster(ene,pCluster+p0);
	m_CsIClusterList->add(cl);
	cl->setEneLayer(eneLayer);
	cl->setPosLayer(pLayer);

	// r138 correction
	//	double fact = 9.35194909495273952e-01;

	double fact = 1.;
	
	// Leakage correction
	double eleak = Leak(fact*ene,fact*eneLayer[nLayers-1])+fact*ene;
	
	// iteration
	//	eleak = Leak(eleak,fact*eneLayer[nLayers-1])+fact*ene;	


	cl->setEneLeak(eleak);

	// defines global variable to be used for fcn
		g_elayer.resize(nLayers);
	for ( i =0;i<nLayers;i++)
	{
		// We are working in GeV
		g_elayer[i] = fact*eneLayer[i]/1000.;
	}
	nbins = nLayers;
	slope = 1;  // temporary

	// Do profile fitting

	Profile(ene,cl);

	m_CsIClusterList->writeOut();

	return sc;
}

/*! Finalization of algorithm 
 *  - minuit deleted
 *  -  CsIClusterList written
 */
//################################################
StatusCode CalClustersAlg::finalize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
	
	
	// delete Minuit object
  delete minuit;


//	m_CsIClusterList->writeOut();

	return sc;
}




