/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/TkrRecon/TkrTree.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include "CalUtil/ICalClusterHitTool.h"

#include "CalUtil/CalDefs.h"

#include <CalRecon/ICalEnergyCorr.h>
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// to access an XML containing Profile Bias parameters file
#include "xmlBase/IFile.h"

#include "FullShowerDevelopmentDescriptionManager.h"
#include "FullShowerProfileParamsManager.h"

//Gamma function and Minuit
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"
#include "TFile.h"

/**   
* @class CalFullProfileTool
*
* Algorithm for calculating energy by fitting the longitudinal
* shower profile using a full (= longitudinal AND radial) description of the shower development in the calorimeter.
*
*
* $Header$
*/


class CalFullProfileTool  : public AlgTool, virtual public ICalEnergyCorr 
{
public:
  //! destructor
  CalFullProfileTool( const std::string& type, const std::string& name, const IInterface* parent);
  ~CalFullProfileTool() {}; 
  
  StatusCode initialize();
  
  //! Longitudinal profile fitting method
  /*! It performs a longitudinal profile fitting using :
   * the FullShowerDevelopmentDescriptionManager to determine the energy deposition in the layers 
   * and the FullShowerProfileParamsManager to constrain the fit parameters during the fitting procedure.
   *
   * \b Revision:
   * - 07/06/05       Philippe Bruel    first implementation
   */     

  Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrTree* );
  
  int doProfileFit(double *pp, double *vv, double tkr_RLn, MsgStream lm, int optntrye);
  
  double GetRadiationLengthInTracker(Event::TkrTree*);
  
  int DetectSaturation(Event::CalCluster* cluster);
  
  StatusCode finalize();
    
private:  
  /// Detector Service
  IGlastDetSvc* m_detSvc; 
  
  /// G4 Propagator tool
  IPropagator* m_G4PropTool; 
    
  /// TkrGeometrySvc used for access to tracker geometry info
  ITkrGeometrySvc* m_tkrGeom;

  /// Tool to loop over the xtals in a cluster.
  ICalClusterHitTool* m_calClusterHitTool;

  // in order to handle saturation
  static int m_Nsaturated;
  
  double m_eTotal;

  // some results of the fit
  double m_amin; // minimum of minimized function
  double m_par0; // alpha
  double m_par1; // tmax
  double m_par2; // energy
  double m_epar0;
  double m_epar1;
  double m_epar2;
  double m_wideningfactor;
  double m_lastx0;
  double m_totx0cal;
  double m_ierflg;
  
  TMinuit* m_minuit;
  static FullShowerProfileParamsManager *m_fsppm;
  static FullShowerDevelopmentDescriptionManager *m_fsddm;
  static double m_elayer_dat[8];
  static double m_elayer_datsat[8];
  static double m_elayer_nsat[8];
  static double m_elayer_fit[8];
  static double m_eelayer_fit[8];
  // function passed to Minuit to minimize
  static void fcn(int & , double *, double &f, double *par, int );
  static double compute_chi2(double *par);
  static double compute_deposited_energy(double *par, double z0, double z1);
  
  static double m_chisq;
  static double m_params_contribution;
  static double m_params_contribution_factor;
  static double m_totchisq;
  static int m_optpr;
  
  static double m_spy_par[3];
  static double m_spy_totchisq;

private:
  double BIAS0;
  double BIAS1;
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalFullProfileTool) ;

int CalFullProfileTool::m_Nsaturated = 0;
FullShowerDevelopmentDescriptionManager * CalFullProfileTool::m_fsddm = 0;
FullShowerProfileParamsManager * CalFullProfileTool::m_fsppm = 0;
double CalFullProfileTool::m_elayer_dat[8];
double CalFullProfileTool::m_elayer_datsat[8];
double CalFullProfileTool::m_elayer_nsat[8];
double CalFullProfileTool::m_elayer_fit[8];
double CalFullProfileTool::m_eelayer_fit[8];
double CalFullProfileTool::m_chisq = 0;
double CalFullProfileTool::m_params_contribution = 0;
double CalFullProfileTool::m_params_contribution_factor = 0;
double CalFullProfileTool::m_totchisq = 0;
int CalFullProfileTool::m_optpr = 0;
double CalFullProfileTool::m_spy_par[3];
double CalFullProfileTool::m_spy_totchisq = 0;

double CalFullProfileTool::compute_deposited_energy(double *par, double z0, double z1)
{
  double beta = (par[0]-1)/par[1];
  return par[2]*(TMath::Gamma(par[0],beta*z1)-TMath::Gamma(par[0],beta*z0));
}

double CalFullProfileTool::compute_chi2(double *par)
{
  if(par[2]<=0.) return 1e10;

  m_params_contribution = m_fsppm->GetChi2Contribution(par);
  if(m_params_contribution<0) return 1e10;
  
  m_fsddm->FillCurrentFSDD(par[1]);

  int i,j;
  
  // fill TMeF
  for(i=0;i<8;++i) m_elayer_fit[i] = 0;
  double e_i;
  for(i=0;i<m_fsddm->CurrentFSDD->NStep;++i)
    {
      e_i = compute_deposited_energy(par,m_fsddm->CurrentFSDD->X0[i],m_fsddm->CurrentFSDD->X0[i+1]);
      for(j=0;j<8;++j) m_elayer_fit[j] += e_i*m_fsddm->CurrentFSDD->layerfraction[j][i];
    }

  double maxlayerenergy = -99999.;
  for(i=0;i<8;++i) if(m_elayer_fit[i]>maxlayerenergy) maxlayerenergy = m_elayer_fit[i];

  double fit_error = 0.1;
  if(maxlayerenergy>0) fit_error = maxlayerenergy*m_fsppm->relerr;

  // calculate chisquare
  m_chisq = 0;
  double delta;
  for(i=0;i<8;++i)
    {
      delta  = (m_elayer_dat[i]-m_elayer_fit[i])/fit_error;
      if(m_elayer_nsat[i]>0 && m_elayer_fit[i]>m_elayer_dat[i])
        delta = 0;
      m_chisq += delta*delta;
//       if(m_optpr)
//         printf("CalFullProfileTool::compute_chi2 m_elayer_dat[%d] = %f %f -> %f -> %f\n",i,m_elayer_dat[i],m_elayer_fit[i],delta,m_chisq);
    }

  m_params_contribution_factor = 3+(double)m_Nsaturated;
  if(par[1]/m_fsddm->CurrentFSDD->lastx0>0.9)
    {
      m_params_contribution_factor += (par[1]/m_fsddm->CurrentFSDD->lastx0 - 0.9);
      if(m_params_contribution_factor>10) m_params_contribution_factor = 10.;
    }

  m_totchisq = m_chisq + m_params_contribution_factor * m_params_contribution;

//   if(m_optpr)
//     printf("CalFullProfileTool::compute_chi2 m_totchisq = %f (%f) (%f GeV)\n",m_totchisq,m_chisq,par[2]);

  if(m_totchisq<m_spy_totchisq)
    {
      m_spy_totchisq = m_totchisq;
      m_spy_par[0] = par[0];
      m_spy_par[1] = par[1];
      m_spy_par[2] = par[2];
    }

  return m_totchisq;
}

void CalFullProfileTool::fcn(int & , double *, double &f, double *par, int )
{
  f = compute_chi2(par);
}

CalFullProfileTool::CalFullProfileTool( const std::string& type, 
                                      const std::string& name, 
                                      const IInterface* parent)
  : AlgTool(type,name,parent)
{
  // declare base interface for all consecutive concrete classes
  declareInterface<ICalEnergyCorr>(this);
  // Declare the properties that may be set in the job options file
}

StatusCode CalFullProfileTool::initialize()
  // This function does following initialization actions:
  //    - extracts geometry constants from xml file using GlastDetSvc
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::DEBUG << "Initializing CalFullProfileTool" <<endreq;
    
    if ((sc = service("GlastDetSvc", m_detSvc, true)).isFailure())
    { 
      throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
  
    // find TkrGeometrySvc service
    if ((sc = service("TkrGeometrySvc", m_tkrGeom, true)).isFailure())
    {
      throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    IToolSvc* toolSvc = 0;
    if((sc = service("ToolSvc", toolSvc, true)).isFailure()) 
    {
        throw GaudiException("Service [ToolSvc] not found", name(), sc);
    }

    if((sc = toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)).isFailure()) 
    {
        throw GaudiException("Tool [G4PropagationTool] not found", name(), sc);
        log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    // Find CalClusterHitTool
    if ((sc = toolSvc->retrieveTool("CalClusterHitTool", m_calClusterHitTool)).isFailure()) {
      throw GaudiException("Tool [CalClusterHitTool] not found", name(), sc);
      log << MSG::ERROR << "Couldn't find the CalClusterHitTool!" << endreq;
      return StatusCode::FAILURE;
    }
    
    BIAS0 = 2.655362;
    BIAS1 = 0.007235;

    m_fsppm = new FullShowerProfileParamsManager();
    m_fsddm = new FullShowerDevelopmentDescriptionManager(m_detSvc,14,2.,1.,1.85,0.80,0.1);
    
    // Minuit object
    m_minuit = new TMinuit(5);
    
    //Sets the function to be minimized
    m_minuit->SetFCN(fcn);

    //m_saturationadc = 4060;

    return sc;
}

double CalFullProfileTool::GetRadiationLengthInTracker(Event::TkrTree* tree)
{
    double tkr_RLn = 0.;

    if(!tree) return tkr_RLn;
    
    Point  x0 = tree->getAxisParams()->getEventPosition();
    Vector vv = tree->getAxisParams()->getEventAxis();
    double costheta = fabs(vv.z());
    
    // copied from CalValsCorrTool
    // Start at the top most point (ie in the top silicon layer) 
    int    topLayer = tree->getHeadNode()->front()->getTreeStartLayer();
    double zTopLyr  = std::max(m_tkrGeom->getLayerZ(topLayer, 0), m_tkrGeom->getLayerZ(topLayer, 1));

    // Translate x0 to this z position
    double arcLen   = vv.z() > 0. ? (zTopLyr - x0.z()) / vv.z() : 0.;

    // Translate the tree start position to middle of top silicon layer
    x0 = x0 + arcLen * vv;

    // Now set up and call propagator to get the radiation lengths to calorimeter
    arcLen = vv.z() > 0. ? (x0.z() - m_tkrGeom->calZTop()) / vv.z() : 0.;
    m_G4PropTool->setStepStart(x0, -vv);
    m_G4PropTool->step(arcLen);
    tkr_RLn = m_G4PropTool->getRadLength(); 

    // Patch for error in KalFitTrack: 1/2 of first radiator left out
    int topPlane = m_tkrGeom->getPlane(zTopLyr);

    if (m_tkrGeom->isTopPlaneInLayer(topPlane)) 
    {
        tkr_RLn += 0.5*m_tkrGeom->getRadLenConv(topLayer) / vv.z();
    }
    
    // add up the rad lens; this could be a local array if you're bothered by the overhead
    //   but hey, compared to the propagator...
    double tkr_radLen_nom = 0.; 
    int layerCount = topLayer;
    for(; layerCount>=0; --layerCount) {
      tkr_radLen_nom += m_tkrGeom->getRadLenConv(layerCount) 
        + m_tkrGeom->getRadLenRest(layerCount);
    }
    if(costheta!=0)
      {
        tkr_radLen_nom /= costheta;
        if(tkr_RLn > tkr_radLen_nom * 1.5)     {tkr_RLn = tkr_radLen_nom * 1.5;}
        else if(tkr_RLn < tkr_radLen_nom * .5) {tkr_RLn  = tkr_radLen_nom * .5;}
      }
    
    return tkr_RLn;
}

Event::CalCorToolResult* CalFullProfileTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrTree* tree)
  //               This function fits the parameters of shower profile using
  //               the Minuit minimization package and stores the fitted
  //               parameters in the CalCluster object
  //
{
    Event::CalCorToolResult* corResult = 0;  
    MsgStream lm(msgSvc(), name());
    
    if (!cluster)
    {
        lm << MSG::DEBUG << "Ending doEnergyCorr: No Cluster" 
           << endreq;
        return corResult;
    }

    double pp[3];
    double vv[3];
    double tkr_RLn = 0;

    int i;
    
    m_eTotal = cluster->getMomParams().getEnergy()/1000.;

    if( m_eTotal<1.)
      {
        lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : m_eTotal<1GeV -> no energy computation" <<endreq;
        return corResult;
      }
    
    // defines global variable to be used for fcn
    for (i=0;i<8;++i)
      {
        // We are working in GeV
        m_elayer_dat[i] = (*cluster)[i].getEnergy()/1000.;
      }

    // Detect saturation must be called before m_fsddm->Compute !!!
    DetectSaturation(cluster);

    tkr_RLn = 0;
    if(tree!=NULL)
      tkr_RLn = GetRadiationLengthInTracker(tree);

    // use cal position and direction
    pp[0] = (cluster->getPosition()).x();
    pp[1] = (cluster->getPosition()).y();
    pp[2] = (cluster->getPosition()).z();
    if((cluster->getDirection()).z()>0)
      {
        vv[0] = -(cluster->getDirection()).x();
        vv[1] = -(cluster->getDirection()).y();
        vv[2] = -(cluster->getDirection()).z();
      }
    else
      {
        vv[0] = (cluster->getDirection()).x();
        vv[1] = (cluster->getDirection()).y();
        vv[2] = (cluster->getDirection()).z();
      }

    int CALFIT = 0;
    if(tree!=NULL)
      CALFIT = doProfileFit(pp,vv,tkr_RLn,lm,0);
    else
      CALFIT = doProfileFit(pp,vv,tkr_RLn,lm,1);

    double CALFIT_fit_energy = 1000.*m_par2;
    double CALFIT_energy_err = 1000.*m_epar2;
    double CALFIT_fitflag = m_ierflg;
    double CALFIT_alpha = m_par0;
    double CALFIT_tmax = m_par1;
    double CALFIT_tkr_RLn = tkr_RLn;
    double CALFIT_lastx0 = m_lastx0;
    double CALFIT_cal_eff_RLn = m_totx0cal;
    double CALFIT_totchisq = m_totchisq;
    double CALFIT_chisq = m_chisq;
    double CALFIT_parcf = m_params_contribution_factor;
    double CALFIT_parc = m_params_contribution;
    double CALFIT_recp0 = pp[0];
    double CALFIT_recp1 = pp[1];
    double CALFIT_recp2 = pp[2];
    double CALFIT_recv0 = vv[0];
    double CALFIT_recv1 = vv[1];
    double CALFIT_recv2 = vv[2];
    double CALFIT_widening = m_wideningfactor;

    // Fill TKRFIT_ variables with CALFIT results before trying TKRFIT
    double TKRFIT_fit_energy = 1000.*m_par2;
    double TKRFIT_energy_err = 1000.*m_epar2;
    double TKRFIT_fitflag = m_ierflg;
    double TKRFIT_alpha = m_par0;
    double TKRFIT_tmax = m_par1;
    double TKRFIT_tkr_RLn = tkr_RLn;
    double TKRFIT_lastx0 = m_lastx0;
    double TKRFIT_cal_eff_RLn = m_totx0cal;
    double TKRFIT_totchisq = m_totchisq;
    double TKRFIT_chisq = m_chisq;
    double TKRFIT_parcf = m_params_contribution_factor;
    double TKRFIT_parc = m_params_contribution;
    double TKRFIT_recp0 = pp[0];
    double TKRFIT_recp1 = pp[1];
    double TKRFIT_recp2 = pp[2];
    double TKRFIT_recv0 = vv[0];
    double TKRFIT_recv1 = vv[1];
    double TKRFIT_recv2 = vv[2];
    double TKRFIT_widening = m_wideningfactor;

    int TKRFIT = 0;
    if(tree!=NULL)
    {
        pp[0] =  tree->getAxisParams()->getEventPosition().x();
        pp[1] =  tree->getAxisParams()->getEventPosition().y();
        pp[2] =  tree->getAxisParams()->getEventPosition().z();
        vv[0] = -tree->getAxisParams()->getEventAxis().x();
        vv[1] = -tree->getAxisParams()->getEventAxis().y();
        vv[2] = -tree->getAxisParams()->getEventAxis().z();

        TKRFIT = doProfileFit(pp,vv,tkr_RLn,lm,1);

        if(TKRFIT>0)
          {
            TKRFIT_fit_energy = 1000.*m_par2;
            TKRFIT_energy_err = 1000.*m_epar2;
            TKRFIT_fitflag = m_ierflg;
            TKRFIT_alpha = m_par0;
            TKRFIT_tmax = m_par1;
            TKRFIT_tkr_RLn = tkr_RLn;
            TKRFIT_lastx0 = m_lastx0;
            TKRFIT_cal_eff_RLn = m_totx0cal;
            TKRFIT_totchisq = m_totchisq;
            TKRFIT_chisq = m_chisq;
            TKRFIT_parcf = m_params_contribution_factor;
            TKRFIT_parc = m_params_contribution;
            TKRFIT_recp0 = pp[0];
            TKRFIT_recp1 = pp[1];
            TKRFIT_recp2 = pp[2];
            TKRFIT_recv0 = vv[0];
            TKRFIT_recv1 = vv[1];
            TKRFIT_recv2 = vv[2];
            TKRFIT_widening = m_wideningfactor;
          }
      }
      
    if(CALFIT==0 && TKRFIT==0)
      {
        lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : No result with cal direction neither with tkr direction." <<endreq;
        return corResult;
      }

    // Ok, fill in the corrected information and exit
    Event::CalMomParams momParams = cluster->getMomParams();
    corResult = new Event::CalCorToolResult();
    corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
    corResult->setCorrectionName(type());
    // if TKRFIT has failed, TKRFIT_ variables have already been filled with CALFIT results
    momParams.setEnergy(TKRFIT_fit_energy);
    momParams.setEnergyErr(TKRFIT_energy_err);
    corResult->setParams(momParams);
    corResult->setChiSquare(TKRFIT_totchisq);

    corResult->insert(Event::CalCorEneValuePair("fit_energy", TKRFIT_fit_energy)); // energy (in MeV)
    corResult->insert(Event::CalCorEneValuePair("energy_err", TKRFIT_energy_err)); // energy error (in MeV)
    corResult->insert(Event::CalCorEneValuePair("fitflag",    TKRFIT_fitflag)); // flag from Minuit
    corResult->insert(Event::CalCorEneValuePair("alpha",      TKRFIT_alpha)); // alpha (profile fit parameter 0)
    corResult->insert(Event::CalCorEneValuePair("tmax",       TKRFIT_tmax)); // tmax (profile fit parameter 1)
    corResult->insert(Event::CalCorEneValuePair("tkr_RLn",    TKRFIT_tkr_RLn)); // radiation length in tracker (0 for cal only events)
    corResult->insert(Event::CalCorEneValuePair("lastx0",     TKRFIT_lastx0)); // total radiation length seen by the shower
    corResult->insert(Event::CalCorEneValuePair("cal_eff_RLn",TKRFIT_cal_eff_RLn)); // effective radiation length in cal (i.e in CsI)
    corResult->insert(Event::CalCorEneValuePair("totchisq",   TKRFIT_totchisq)); // total chisquare = usual chisquare + weight * parameters constraint
    corResult->insert(Event::CalCorEneValuePair("chisq",      TKRFIT_chisq)); // usual chisquare (sum_layers ( (e-efit)/de )^2
    corResult->insert(Event::CalCorEneValuePair("parcf",      TKRFIT_parcf)); // weighting factor of the parameters constraint
    corResult->insert(Event::CalCorEneValuePair("parc",       TKRFIT_parc)); // parameters constraint contribution to the chisquare
    corResult->insert(Event::CalCorEneValuePair("recp0",      TKRFIT_recp0)); // Point and direction information - kept for the moment for debugging purpose
    corResult->insert(Event::CalCorEneValuePair("recp1",      TKRFIT_recp1));
    corResult->insert(Event::CalCorEneValuePair("recp2",      TKRFIT_recp2));
    corResult->insert(Event::CalCorEneValuePair("recv0",      TKRFIT_recv0));
    corResult->insert(Event::CalCorEneValuePair("recv1",      TKRFIT_recv1));
    corResult->insert(Event::CalCorEneValuePair("recv2",      TKRFIT_recv2));
    corResult->insert(Event::CalCorEneValuePair("widening",   TKRFIT_widening));

    if(CALFIT>0)
      {
        corResult->insert(Event::CalCorEneValuePair("calfit_fit_energy", CALFIT_fit_energy)); // energy (in MeV)
        corResult->insert(Event::CalCorEneValuePair("calfit_energy_err", CALFIT_energy_err)); // energy error (in MeV)
        corResult->insert(Event::CalCorEneValuePair("calfit_fitflag",    CALFIT_fitflag)); // flag from Minuit
        corResult->insert(Event::CalCorEneValuePair("calfit_alpha",      CALFIT_alpha)); // alpha (profile fit parameter 0)
        corResult->insert(Event::CalCorEneValuePair("calfit_tmax",       CALFIT_tmax)); // tmax (profile fit parameter 1)
        corResult->insert(Event::CalCorEneValuePair("calfit_tkr_RLn",    CALFIT_tkr_RLn)); // radiation length in tracker (0 for cal only events)
        corResult->insert(Event::CalCorEneValuePair("calfit_lastx0",     CALFIT_lastx0)); // total radiation length seen by the shower
        corResult->insert(Event::CalCorEneValuePair("calfit_cal_eff_RLn",CALFIT_cal_eff_RLn)); // effective radiation length in cal (i.e in CsI)
        corResult->insert(Event::CalCorEneValuePair("calfit_totchisq",   CALFIT_totchisq)); // total chisquare = usual chisquare + weight * parameters constraint
        corResult->insert(Event::CalCorEneValuePair("calfit_chisq",      CALFIT_chisq)); // usual chisquare (sum_layers ( (e-efit)/de )^2
        corResult->insert(Event::CalCorEneValuePair("calfit_parcf",      CALFIT_parcf)); // weighting factor of the parameters constraint
        corResult->insert(Event::CalCorEneValuePair("calfit_parc",       CALFIT_parc)); // parameters constraint contribution to the chisquare
        corResult->insert(Event::CalCorEneValuePair("calfit_recp0",      CALFIT_recp0)); // Point and direction information - kept for the moment for debugging purpose
        corResult->insert(Event::CalCorEneValuePair("calfit_recp1",      CALFIT_recp1));
        corResult->insert(Event::CalCorEneValuePair("calfit_recp2",      CALFIT_recp2));
        corResult->insert(Event::CalCorEneValuePair("calfit_recv0",      CALFIT_recv0));
        corResult->insert(Event::CalCorEneValuePair("calfit_recv1",      CALFIT_recv1));
        corResult->insert(Event::CalCorEneValuePair("calfit_recv2",      CALFIT_recv2));
        corResult->insert(Event::CalCorEneValuePair("calfit_widening",   CALFIT_widening));
      }
    else
      {
        corResult->insert(Event::CalCorEneValuePair("calfit_fit_energy", 0)); // energy (in MeV)
        corResult->insert(Event::CalCorEneValuePair("calfit_energy_err", 0)); // energy error (in MeV)
        corResult->insert(Event::CalCorEneValuePair("calfit_fitflag",    0)); // flag from Minuit
        corResult->insert(Event::CalCorEneValuePair("calfit_alpha",      0)); // alpha (profile fit parameter 0)
        corResult->insert(Event::CalCorEneValuePair("calfit_tmax",       0)); // tmax (profile fit parameter 1)
        corResult->insert(Event::CalCorEneValuePair("calfit_tkr_RLn",    0)); // radiation length in tracker (0 for cal only events)
        corResult->insert(Event::CalCorEneValuePair("calfit_lastx0",     0)); // total radiation length seen by the shower
        corResult->insert(Event::CalCorEneValuePair("calfit_cal_eff_RLn",0)); // effective radiation length in cal (i.e in CsI)
        corResult->insert(Event::CalCorEneValuePair("calfit_totchisq",   -999)); // total chisquare = usual chisquare + weight * parameters constraint
        corResult->insert(Event::CalCorEneValuePair("calfit_chisq",      0)); // usual chisquare (sum_layers ( (e-efit)/de )^2
        corResult->insert(Event::CalCorEneValuePair("calfit_parcf",      0)); // weighting factor of the parameters constraint
        corResult->insert(Event::CalCorEneValuePair("calfit_parc",       0)); // parameters constraint contribution to the chisquare
        corResult->insert(Event::CalCorEneValuePair("calfit_recp0",      0)); // Point and direction information - kept for the moment for debugging purpose
        corResult->insert(Event::CalCorEneValuePair("calfit_recp1",      0));
        corResult->insert(Event::CalCorEneValuePair("calfit_recp2",      0));
        corResult->insert(Event::CalCorEneValuePair("calfit_recv0",      0));
        corResult->insert(Event::CalCorEneValuePair("calfit_recv1",      0));
        corResult->insert(Event::CalCorEneValuePair("calfit_recv2",      0));
        corResult->insert(Event::CalCorEneValuePair("calfit_widening",   0));
      }

    return corResult;
}

int CalFullProfileTool::doProfileFit(double *pp, double *vv, double tkr_RLn, MsgStream lm, int optntrye)
  //
  //               This function fits the parameters of shower profile using
  //               the Minuit minimization package and stores the fitted
  //               parameters in the CalCluster object
{

  int i;
  
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : pp = " <<  pp[0] << " " << pp[1] << " " << pp[2] << " " <<endreq;
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : vv = " <<  vv[0] << " " << vv[1] << " " << vv[2] << " " <<endreq;
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : tkr_RLn = " << tkr_RLn <<endreq;
  

  double wideningfactor = 1.;
  double mytheta;
  double tkrrln1; // when wideningfactor = 1
  if(fabs(vv[2])>0)
    {
      mytheta = atan(-sqrt(vv[0]*vv[0]+vv[1]*vv[1])/vv[2])*180./TMath::Pi();
      if(mytheta>40)
        {
          tkrrln1 = 2-0.05*(mytheta-40.);
          if(tkr_RLn>tkrrln1)
            {
              wideningfactor = 1.+(tkr_RLn-tkrrln1)/1.5;
              if(wideningfactor>2) wideningfactor = 2.;
            }
        }
    }

  m_fsddm->SetWideningFactor(wideningfactor);
  m_wideningfactor = wideningfactor;

  if(!m_fsddm->Compute(pp,vv,tkr_RLn))
    {
      lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : Problem during m_fsddm->Compute. Returning empty corResult." <<endreq;
      return 0;
    }

  if(m_fsddm->mintotx0cal<0.5) // Make sure that lastx0>0 during the fitting procedure
    {
      lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr :  m_fsddm->mintotx0cal<0.5. No energy computation. Returning empty corResult." <<endreq;
      return 0;
    }
  
  for(i=0;i<=m_fsddm->NDevelopment;++i)
    {
      lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : m_fsddm->FSDDX0[" << i <<"] :"
         << " startx0 = " << m_fsddm->FSDDX0[i]->startx0
         << " totx0cal = " << m_fsddm->FSDDX0[i]->totx0cal
         << " lastx0 = " << m_fsddm->FSDDX0[i]->lastx0
         <<endreq;
    }

  for(i=0;i<8;++i)
    {
      lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : m_fsddm->meantotx0lay[" << i <<"] :"
         << m_fsddm->meantotx0lay[i]
         << " meanposx0lay = " << m_fsddm->meanposx0lay[i]
         <<endreq;
    }

  double arglist[10];
  int ierflg = 0;
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  
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
  
  double vstart[3];
  double vstep[3];
  
  double par0,par1,par2;
  double epar0,epar1,epar2;
  
  int ntrye = 10;
  double enemin = m_eTotal;
  double enemax = m_eTotal*10;
  double enestep = (enemax-enemin)/(double)ntrye;
  
  double mychisq[500];
  double myalpha[500];
  double mytmax[500];
  double myenergy[500];
  
  m_optpr = 0;

  double scan0,scan1,scan2;
  double mymin = 99999999;

  if(optntrye==1)
    {
      for(i=0;i<=ntrye;++i)
        {
          vstart[2] = enemin+enestep*((double)i);
          m_fsppm->Fill(vstart[2]);
          vstart[0] = m_fsppm->alpha;
          vstart[1] = m_fsppm->tmax;
          vstep[0] = 0.05;
          vstep[1] = 0.05;
          vstep[2] = vstart[2]*0.05;
          m_minuit->mnparm(0, "a1", vstart[0], vstep[0], 1.1,20,ierflg);
          m_minuit->mnparm(1, "a2", vstart[1], vstep[1], 0.1,20,ierflg);
          m_minuit->mnparm(2, "a3", vstart[2], vstep[2], m_eTotal*0.95,10000,ierflg);
          m_minuit->FixParameter(2);
          arglist[0] = 500;
          arglist[1] = 1.;
          m_optpr = 0;
          m_spy_totchisq = 99999999;
          m_spy_par[0] = -1;
          m_spy_par[1] = -1;
          m_spy_par[2] = -1;
          m_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
          m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
          m_minuit->GetParameter(0,par0,epar0);
          m_minuit->GetParameter(1,par1,epar1);
          if(ierflg==4)
            {
              vstart[0] = m_spy_par[0];
              vstart[1] = m_spy_par[1];
              vstep[0] = 0.05;
              vstep[1] = 0.05;
              m_minuit->mnparm(0, "a1", vstart[0], vstep[0], 1.1,20,ierflg);
              m_minuit->mnparm(1, "a2", vstart[1], vstep[1], 0.1,20,ierflg);
              m_minuit->mnparm(2, "a3", vstart[2], vstep[2], m_eTotal*0.95,10000,ierflg);
              m_minuit->FixParameter(2);
              arglist[0] = 500;
              arglist[1] = 1.;
              m_optpr = 0;
              m_spy_totchisq = 99999999;
              m_spy_par[0] = -1;
              m_spy_par[1] = -1;
              m_spy_par[2] = -1;
              m_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
              m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
              //          m_minuit->mnprin(3,amin);
              m_minuit->GetParameter(0,par0,epar0);
              m_minuit->GetParameter(1,par1,epar1);
            }
          mychisq[i] = amin;
          myalpha[i] = par0;
          mytmax[i] = par1;
          myenergy[i] = vstart[2];
          if(amin<mymin)
            {
              mymin = amin;
              scan0 = par0;
              scan1 = par1;
              scan2 = vstart[2];
            }
          if(amin>5*mymin) break;
        }
    }
  else
    scan2 = m_eTotal;

  m_fsppm->Fill(scan2);
  vstart[0] = scan0;
  vstart[1] = scan1;
  vstart[2] = scan2;

  vstep[0] = 0.05;
  vstep[1] = 0.05;
  vstep[2] = scan2*0.05;
  
  m_minuit->mnparm(0, "a1", vstart[0], vstep[0], 1.1,20,ierflg);
  m_minuit->mnparm(1, "a2", vstart[1], vstep[1], 0.1,20,ierflg);
  m_minuit->mnparm(2, "a3", vstart[2], vstep[2], m_eTotal*0.95,10000,ierflg);
  m_minuit->Release(2);
  
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : Initial parameters : " << vstart[0] << " " << vstart[1] << " " << vstart[2] <<endreq;

  // Calls Migrad with 500 iterations maximum
  arglist[0] = 500;
  arglist[1] = 1.;
  m_optpr = 0;

  m_spy_totchisq = 99999999;
  m_spy_par[0] = -1;
  m_spy_par[1] = -1;
  m_spy_par[2] = -1;
  m_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  if(ierflg==4)
    {
      //      m_minuit->mnprin(3,amin);
      vstart[0] = m_spy_par[0];
      vstart[1] = m_spy_par[1];
      vstart[2] = m_spy_par[2];
      vstep[0] = 0.05;
      vstep[1] = 0.05;
      vstep[2] = m_spy_par[2]*0.05;
      m_minuit->mnparm(0, "a1", vstart[0], vstep[0], 1.1,20,ierflg);
      m_minuit->mnparm(1, "a2", vstart[1], vstep[1], 0.1,20,ierflg);
      m_minuit->mnparm(2, "a3", vstart[2], vstep[2], m_eTotal*0.95,10000,ierflg);
      arglist[0] = 500;
      arglist[1] = 1.;
      m_optpr = 0;
      m_spy_totchisq = 99999999;
      m_spy_par[0] = -1;
      m_spy_par[1] = -1;
      m_spy_par[2] = -1;
      m_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      m_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    }

  //  m_minuit->mnprin(3,amin);

  m_amin = amin;
  m_ierflg = (double)ierflg;

  m_minuit->GetParameter(0,par0,epar0);
  m_minuit->GetParameter(1,par1,epar1);
  m_minuit->GetParameter(2,par2,epar2);

  m_optpr = 1;
  vstart[0] = par0;
  vstart[1] = par1;
  vstart[2] = par2;
  compute_chi2(vstart);
  m_optpr = 0;
 
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : ierflg = " << ierflg <<endreq;
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : totchisq = " << m_totchisq << " = " << m_chisq << " + " << m_params_contribution_factor << " * " << m_params_contribution <<endreq;
  
  m_lastx0 = m_fsddm->CurrentFSDD->lastx0;
  m_totx0cal = m_fsddm->CurrentFSDD->totx0cal;
    
  // Clear minuit
  arglist[0] = 1;
  m_minuit->mnexcm("CLEAR", arglist ,1,ierflg);
    
  // Quick Bias Correction
  double bias = 1.;
  if(par2>0.1)
    bias = log(BIAS1*log(par2)+BIAS0);
  if(bias<0.97) bias = 0.97;
  if(bias>0.99) bias = 0.99;
  par2 /= bias;
  epar2 /= bias;
  
  lm << MSG::DEBUG << "CalFullProfileTool::doEnergyCorr : results : alpha = " << par0 
     << " tmax = " << par1 
     << " energy = " << par2 
     << " fit flag " << ierflg 
     << " totchisq = " << amin <<endreq;

  m_par0 = par0;
  m_par1 = par1;
  m_par2 = par2;
  m_epar0 = epar0;
  m_epar1 = epar1;
  m_epar2 = epar2;
  
  return 1;
}

int CalFullProfileTool::DetectSaturation(Event::CalCluster* cluster)
{
  // Reset variables.
  int i,j,k;

  for (i=0;i<16;++i) {
    for (j=0;j<8;++j) {
      for (k=0;k<12;++k) {
        m_fsddm->OffSatu[i][j][k] = 0;
      }
    }
  }

  for(j=0;j<8;++j) {
    m_elayer_datsat[j] = 0;
    m_elayer_nsat[j] = 0;
  }

  m_Nsaturated = cluster->getXtalsParams().getNumSaturatedXtals();

  // If there are no saturated xtals there is nothing to do.
  if (m_Nsaturated==0) return 0;

  // Otherwise go ahead and loop over the xtals in the cluster.
  m_calClusterHitTool->fillRecDataVec(cluster);
  std::vector<Event::CalXtalRecData*> xtalList = m_calClusterHitTool->getRecDataVec();
  std::vector<Event::CalXtalRecData*>::const_iterator xtalData;

  for (xtalData = xtalList.begin(); xtalData != xtalList.end(); xtalData++){
    int tower  = (*xtalData)->getTower();
    int layer  = (*xtalData)->getLayer();
    int column = (*xtalData)->getColumn();
    if ((*xtalData)->saturated()){
      m_fsddm->OffSatu[tower][layer][column] = 1;
      m_elayer_datsat[layer] += (*xtalData)->getEnergy()/1000;
      m_elayer_nsat[layer] += 1;
    }
  }

  return 0;
}


  StatusCode CalFullProfileTool::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  
  // delete Minuit object
  delete m_minuit;
  
  delete m_fsddm;
  delete m_fsppm; 
  
  return sc;
}

