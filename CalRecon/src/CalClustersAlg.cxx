
#include "CalRecon/CalClustersAlg.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "CalRecon/gamma.h"
#include "CalRecon/Midnight.h"
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
	gamma1 = Gamma(alpha,x);
	x += dx;

	// gamma2 = integration from 0 to x+dx	
	gamma2 = Gamma(alpha,x);
	
	// the result of integration over Xtal pathlength is E*(gamma2-gamma1)
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


double CalClustersAlg::Leak(double eTotal,double elast)
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
                                     
{
    
    MsgStream log(msgSvc(), name());
    
    // sigma2.z() is useless here no matter its value.
    double cov_xz = 0;  // covariance x,z
    double cov_yz = 0;  // covariance y,z
    double mx=0;        // mean x
    double my=0;        // mean y
    double mz1=0;       // mean z for x pos
    double mz2=0;       // mean z for y pos
    double norm1=0;
    double norm2=0;
    double var_z1=0;    // variance of z	
    double var_z2=0;
    int nlx=0,nly=0;
    Vector nodir(-1000.,-1000.,-1000.);
    for(int il=0;il<nlayers;il++)
    {                
        // For the moment forget about longitudinal position
        if(il%2==1)
        {
            if (sigma2[il].x()>0.)
            {
                nlx++;
                double err = 1/sigma2[il].x();
                cov_xz += pos[il].x()*pos[il].z()*err;
                var_z1 += pos[il].z()*pos[il].z()*err;
                mx += pos[il].x()*err;
                mz1 += pos[il].z()*err;
                norm1 += err;
            }
        }
        else
        {
            if(sigma2[il].y()>0.)
            {
                
                nly++;
                double err = 1/sigma2[il].y();
                cov_yz += pos[il].y()*pos[il].z()*err;
                var_z2 += pos[il].z()*pos[il].z()*err;
                my += pos[il].y()*err;
                mz2 += pos[il].z()*err;
                norm2 += err;
            }
        }
    }		
    
    
    if(nlx <2 || nly < 2 )return nodir;
    
    
    
    
    mx /= norm1;
    mz1 /= norm1;
    cov_xz /= norm1;
    cov_xz -= mx*mz1;
    var_z1 /= norm1;
    var_z1 -= mz1*mz1;
    if(var_z1 == 0) return nodir;
    double tgthx = cov_xz/var_z1;
    
    
    my /= norm2;
    mz2 /= norm2;
    cov_yz /= norm2;
    cov_yz -= my*mz2;
    var_z2 /= norm2;
    var_z2 -= mz2*mz2;
    
    // Now we have cov(x,z) and var(z) we can deduce slope
    if(var_z2 == 0) return nodir;
    double tgthy = cov_yz/var_z2;
    
    
    double tgtheta_sqr = tgthx*tgthx+tgthy*tgthy;
    double costheta = 1/sqrt(1+tgtheta_sqr);
    Vector dir(costheta*tgthx,costheta*tgthy,costheta);
    return dir;
}


CalClustersAlg::CalClustersAlg(const std::string& name,
                               ISvcLocator* pSvcLocator):
Algorithm(name, pSvcLocator)
{
    declareProperty("callNumber",m_callNumber=0);
    
}

StatusCode CalClustersAlg::initialize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    sc = service("GlastDetSvc", detSvc);
    
    if(sc.isFailure())
    {
        log << MSG::ERROR << "GlastDetSvc could not be found" <<endreq;
        return sc;
    }
    
    double value;
    if(!detSvc->getNumericConstByName(std::string("CALnLayer"), &value)) {
        log << MSG::ERROR << " constant " << " CALnLayer "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } else m_CalnLayers = value;
    
    if(!detSvc->getNumericConstByName(std::string("CsIWidth"),&m_CsIWidth)){
        log << MSG::ERROR << " constant " << " CsIWidth "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    if(!detSvc->getNumericConstByName(std::string("CsIHeight"),&m_CsIHeight)){
        log << MSG::ERROR << " constant " << " CsIHeight "
            <<" not defined" << endreq;
        return StatusCode::FAILURE;
    } 
    
    
    
    
    
    setProperties();
    log << MSG::INFO << "CalClustersAlg: callNumber = " 
        << m_callNumber << endreq;
    
    xtalHeight = m_CsIHeight/10.;  // crystal height in cm
    xtalWidth = m_CsIWidth/10.;    // crystal width in cm
    
    
    // Minuit object
    minuit = new Midnight(5);
    
    //Sets the function to be minimized
    minuit->SetFCN(fcn);
    
    
    g_elayer.clear();
    
    return sc;
}

StatusCode CalClustersAlg::retrieve()
{
    
    StatusCode sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    
    DataObject* pnode=0;
    
    sc = eventSvc()->retrieveObject(EventModel::CalRecon::Event, pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject(EventModel::CalRecon::Event, new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }
    m_calClusterCol = SmartDataPtr<CalClusterCol> (eventSvc(),EventModel::CalRecon::CalClusterCol);
    if (!m_calClusterCol ){
        m_calClusterCol = 0;
        m_calClusterCol = new CalClusterCol();
        sc = eventSvc()->registerObject(EventModel::CalRecon::CalClusterCol,m_calClusterCol);
    } else {
        m_calClusterCol->delClusters();
    }
    
    m_calXtalRecCol = SmartDataPtr<CalXtalRecCol>(eventSvc(),
        EventModel::CalRecon::CalXtalRecCol); 
    
    
    return sc;
}

//################################################
StatusCode CalClustersAlg::execute()
//################################################
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    sc = retrieve();
    const Point p0(0.,0.,0.);
    
    int rectkr=0;  //is tracker recon OK?
    int ntracks;
    Vector trackDirection;
    Point trackVertex;
    
    SmartDataPtr<TkrVertexCol> tkrRecData(eventSvc(),EventModel::TkrRecon::TkrVertexCol);
    if (tkrRecData == 0) {
        log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
        // return sc;
    }
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
    double ene = 0;
    Vector pCluster(p0);
    int nLayers = m_CalnLayers;
    
    std::vector<double> eneLayer(nLayers,0.);
    std::vector<Vector> pLayer(nLayers);
    std::vector<Vector> rmsLayer(nLayers);
    
    
    
    // Compute barycenter and various moments
    
    
    for (CalXtalRecCol::const_iterator it = m_calXtalRecCol->begin();
    it != m_calXtalRecCol->end(); it++){
        CalXtalRecData* recData = *it;
        
        double eneXtal = recData->getEnergy();
        Vector pXtal = recData->getPosition() - p0;
        int layer = (recData->getPackedId()).getLayer();
        
        eneLayer[layer]+=eneXtal;
        
        Vector ptmp = eneXtal*pXtal;
        pLayer[layer] += ptmp;
        
        Vector ptmp_sqr(ptmp.x()*pXtal.x(),ptmp.y()*pXtal.y(),ptmp.z()*pXtal.z());
        rmsLayer[layer] += ptmp_sqr;  // Position error is assumed to be 1/sqrt(eneXtal)
        
        ene  += eneXtal;
        pCluster += ptmp;
    }
    
    // Now take the means
    if(ene>0.)pCluster *= (1./ene); else pCluster=Vector(-1000., -1000., -1000.);
    int i = 0;
    for( ;i<nLayers;i++){
        if(eneLayer[i]>0)
        {
            pLayer[i] *= (1./eneLayer[i]);
            rmsLayer[i] *= (1./eneLayer[i]);
            
            Vector sqrLayer(pLayer[i].x()*pLayer[i].x(),pLayer[i].y()*pLayer[i].y(),pLayer[i].z()*pLayer[i].z());
            
            
            Vector d; 
            if(i%2 == 1) d = Vector(m_CsIWidth*m_CsIWidth/12.,0.,0.);
            else d = Vector(0.,m_CsIWidth*m_CsIWidth/12.,0.);
            
            rmsLayer[i] += d-sqrLayer;
            
        }
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
    
    //slope=1;  // Temporary whilst the algorith appplies to slope=1 showers only.
    
    // Take square roots of RMS
    //	rms_trans = sqrt(rms_trans);
    //	rms_long = sqrt(rms_long);
    
    // Fill CsICluster data
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
    double calTransvOffset = 0.;
    if(ntracks>0){
        Vector calOffset = (p0+pCluster) - trackVertex;
        double calLongOffset = trackDirection*calOffset;
        calTransvOffset =sqrt(calOffset.mag2() - calLongOffset*calLongOffset);
        
    }
    cl->initialize(eleak,eneLayer,pLayer,rmsLayer,rms_long,rms_trans,
        caldir,calTransvOffset);
    
    m_calClusterCol->writeOut(log << MSG::DEBUG);
    
    return sc;
}

//################################################
StatusCode CalClustersAlg::finalize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
	
	
	// delete Minuit object
        delete minuit;



	return sc;
}




