/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IParticlePropertySvc.h"
#include "GaudiKernel/ParticleProperty.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/TkrRecon/TkrTree.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include "CalUtil/CalDefs.h"
#include "Event/Digi/CalDigi.h"

#include <CalRecon/ICalEnergyCorr.h>
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

// to access an XML containing Profile Bias parameters file
#include "xmlBase/IFile.h"

#include "NewFullShowerDevelopmentDescriptionManager.h"
#include "NewFullShowerProfileParamsManager.h"

//Gamma function and Minuit
#include "TMath.h"
#include "TMinuit.h"
#include "TNtuple.h"
#include "TFile.h"

/**   
* @class NewCalFullProfileTool
*
* Algorithm for calculating energy by fitting the longitudinal
* shower profile using a full (= longitudinal AND radial) description of the shower development in the calorimeter.
*
*
* $Header$
*/


class NewCalFullProfileTool  : public AlgTool, virtual public ICalEnergyCorr 
{
public:
  //! destructor
  NewCalFullProfileTool( const std::string& type, const std::string& name, const IInterface* parent);
  ~NewCalFullProfileTool() {}; 
  
  StatusCode initialize();
  
  //! Longitudinal profile fitting method
  /*! It performs a longitudinal profile fitting using :
   * the NewFullShowerDevelopmentDescriptionManager to determine the energy deposition in the layers 
   * and the NewFullShowerProfileParamsManager to constrain the fit parameters during the fitting procedure.
   *
   * \b Revision:
   * - 07/06/05       Philippe Bruel    first implementation
   */     

  Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrTree* );

  int doProfileFit(double *pp, double *vv, double tkr_RLn, MsgStream lm, int optntrye, double logEradprof=-999, double wf=-999, double *inputfitpar=NULL);

  double GetRadiationLengthInTracker(Event::TkrTree*);

  int DetectSaturation();
  double GetCrystalDistanceToTrajectory(Point p0, Vector v0, Point p1, Vector v1, double xtalhalflength);
  int SelectCloseCrystals(double *pptraj, double *vvtraj, double maxdistance);
  double GetChi2Dist(double *pptraj, double *vvtraj, double fiterror, double *result=NULL);

  StatusCode finalize();
    
private:
  /// Pointer to the Gaudi data provider service
  IDataProviderSvc* nm_dataSvc;
  
  /// Detector Service
  IGlastDetSvc *    nm_detSvc; 

  /// G4 Propagator tool
  IPropagator *     nm_G4PropTool; 
  
  /// TkrGeometrySvc used for access to tracker geometry info
  ITkrGeometrySvc*  nm_tkrGeom;

  // in order to handle saturation
  float nm_saturationadc;
  bool nm_saturated[16][8][12];
  double nm_xtal_energy[16][8][12];
  double nm_xtal_doca[16][8][12];
  int nm_layer_nsat[8];
  double nm_elayer_datsat[8];

  double ParXYcor[6];

  double nm_eTotal;

  // some results of the fit
  double nm_amin; // minimum of minimized function
  double nm_par0;
  double nm_par1;
  double nm_par2; // log energy
  double nm_epar0;
  double nm_epar1;
  double nm_epar2;
  double nm_xcor;
  double nm_ycor;
  double nm_energy;
  double nm_alpha;
  double nm_beta;
  double nm_tmax;
  double nm_wideningfactor;
  double nm_lastx0;
  double nm_totx0cal;
  double nm_ierflg;

  TMinuit* nm_minuit;
  static NewFullShowerProfileParamsManager *nm_fsppm;
  static NewFullShowerDevelopmentDescriptionManager *nm_fsddm;
  static double nm_elayer_dat[8];
  static double nm_elayer_datsel[8];
  static double nm_elayer_fit[8];
  static double nm_eelayer_fit[8];
  static int nm_nxtal;
  static int nm_nxtal_sat;
  static int nm_nsatlayer;
  static double nm_toterrfit;
  static double nm_extal_dat[FSDD_XTAL_NMAX];
  static double nm_extal_fit[FSDD_XTAL_NMAX];
  static double nm_eextal_fit[FSDD_XTAL_NMAX];
  static double nm_xtal_dist[FSDD_XTAL_NMAX];
  static int nm_xtal_itow[FSDD_XTAL_NMAX];
  static int nm_xtal_ilay[FSDD_XTAL_NMAX];
  static int nm_xtal_icol[FSDD_XTAL_NMAX];
  static int nm_xtal_satu[FSDD_XTAL_NMAX];
  // function passed to Minuit to minimize
  static void fcn(int & , double *, double &f, double *par, int );
  static double compute_chi2(double *par);
  static double compute_deposited_energy(double *par, double z0, double z1);
  
  static double nm_chisq;
  static double nm_params_contribution;
  static double nm_params_contribution_factor;
  static double nm_totchisq;
  static int nm_optpr;
  static int nm_optselxtalinfit;

  static double nm_spy_par[3];
  static double nm_spy_totchisq;
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(NewCalFullProfileTool) ;

NewFullShowerDevelopmentDescriptionManager * NewCalFullProfileTool::nm_fsddm = 0;
NewFullShowerProfileParamsManager * NewCalFullProfileTool::nm_fsppm = 0;
double NewCalFullProfileTool::nm_elayer_dat[8];
double NewCalFullProfileTool::nm_elayer_datsel[8];
double NewCalFullProfileTool::nm_elayer_fit[8];
double NewCalFullProfileTool::nm_eelayer_fit[8];
int NewCalFullProfileTool::nm_nxtal;
int NewCalFullProfileTool::nm_nxtal_sat;
int NewCalFullProfileTool::nm_nsatlayer;
double NewCalFullProfileTool::nm_toterrfit;
double NewCalFullProfileTool::nm_extal_dat[FSDD_XTAL_NMAX];
double NewCalFullProfileTool::nm_extal_fit[FSDD_XTAL_NMAX];
double NewCalFullProfileTool::nm_eextal_fit[FSDD_XTAL_NMAX];
double NewCalFullProfileTool::nm_xtal_dist[FSDD_XTAL_NMAX];
int NewCalFullProfileTool::nm_xtal_itow[FSDD_XTAL_NMAX];
int NewCalFullProfileTool::nm_xtal_ilay[FSDD_XTAL_NMAX];
int NewCalFullProfileTool::nm_xtal_icol[FSDD_XTAL_NMAX];
int NewCalFullProfileTool::nm_xtal_satu[FSDD_XTAL_NMAX];
double NewCalFullProfileTool::nm_chisq = 0;
double NewCalFullProfileTool::nm_params_contribution = 0;
double NewCalFullProfileTool::nm_params_contribution_factor = 0;
double NewCalFullProfileTool::nm_totchisq = 0;
int NewCalFullProfileTool::nm_optpr = 0;
int NewCalFullProfileTool::nm_optselxtalinfit = 0;
double NewCalFullProfileTool::nm_spy_par[3];
double NewCalFullProfileTool::nm_spy_totchisq = 0;

double NewCalFullProfileTool::compute_deposited_energy(double *par, double z0, double z1)
{
  return par[2]*(TMath::Gamma(par[0],par[1]*z1)-TMath::Gamma(par[0],par[1]*z0));
}

double NewCalFullProfileTool::compute_chi2(double *par)
{
  nm_totchisq = 1e10;

  if(!nm_fsppm->Fill(par)) return nm_totchisq;

  double par0[3];
  par0[0] = nm_fsppm->alpha;
  par0[1] = nm_fsppm->beta;
  par0[2] = nm_fsppm->energy;

  nm_params_contribution = nm_fsppm->GetChi2Contribution(par);
  
  nm_fsddm->FillCurrentFSDD(nm_fsppm->tmax);

  int i,j;
  
  // fill TMeF
  for(i=0;i<8;++i) nm_elayer_fit[i] = 0;
  for(i=0;i<nm_nxtal;++i) nm_extal_fit[i] = 0;

  double e_i;
  for(i=0;i<nm_fsddm->CurrentFSDD->NStep;++i)
    {
      e_i = compute_deposited_energy(par0,nm_fsddm->CurrentFSDD->X0[i],nm_fsddm->CurrentFSDD->X0[i+1]);
      for(j=0;j<8;++j) nm_elayer_fit[j] += e_i*nm_fsddm->CurrentFSDD->layerfraction[j][i];
      for(j=0;j<nm_nxtal_sat;++j) nm_extal_fit[j] += e_i*nm_fsddm->CurrentFSDD->xtalfraction[j][i];
      if(nm_optselxtalinfit) for(j=nm_nxtal_sat;j<nm_nxtal;++j) nm_extal_fit[j] += e_i*nm_fsddm->CurrentFSDD->xtalfraction[j][i];
    }

  // calculate maxlayerenergy before subtracting xtal energies !!!
  double maxlayerenergy = -99999.;
  int imax = -1;
  for(i=0;i<8;++i)
    if(nm_elayer_fit[i]>maxlayerenergy) 
      {
        maxlayerenergy = nm_elayer_fit[i];
        imax = i;
      }
  
  double relerrfit = nm_fsppm->relerr;
  nm_toterrfit = maxlayerenergy*relerrfit;

  for(j=0;j<nm_nxtal_sat;++j) nm_elayer_fit[nm_xtal_ilay[j]] -= nm_extal_fit[j];
  if(nm_optselxtalinfit)
    for(j=nm_nxtal_sat;j<nm_nxtal;++j) nm_elayer_fit[nm_xtal_ilay[j]] -= nm_extal_fit[j];

  // calculate chisquare
  nm_chisq = 0;
  double delta;
  for(i=0;i<8;++i)
    {
      if(!nm_optselxtalinfit)
        delta  = (nm_elayer_dat[i]-nm_elayer_fit[i])/nm_toterrfit;
      else
        delta  = (nm_elayer_datsel[i]-nm_elayer_fit[i])/nm_toterrfit;

      nm_chisq += delta*delta;
      nm_eelayer_fit[i] = delta;
    }

  for(i=0;i<nm_nxtal;++i)
    {
      delta  = (nm_extal_dat[i]-nm_extal_fit[i])/nm_toterrfit;
      //
      nm_eextal_fit[i] = delta;
      if(nm_xtal_satu[i] && nm_extal_fit[i]>nm_extal_dat[i]) delta = 0;
      //
      if(i<nm_nxtal_sat || nm_optselxtalinfit)
        nm_chisq += delta*delta;
    }

  nm_params_contribution_factor = 1.0;
  //  nm_params_contribution_factor += (double)nm_nsatlayer;

  if(nm_fsppm->tmax/nm_fsddm->CurrentFSDD->lastx0>0.9)
    {
      nm_params_contribution_factor += (nm_fsppm->tmax/nm_fsddm->CurrentFSDD->lastx0 - 0.9);
      if(nm_params_contribution_factor>10) nm_params_contribution_factor = 10.;
    }
  
  nm_totchisq = nm_chisq + nm_params_contribution_factor * nm_params_contribution;

  if(nm_totchisq<nm_spy_totchisq)
    {
      nm_spy_totchisq = nm_totchisq;
      nm_spy_par[0] = par[0];
      nm_spy_par[1] = par[1];
      nm_spy_par[2] = par[2];
    }

  return nm_totchisq;
}

void NewCalFullProfileTool::fcn(int & , double *, double &f, double *par, int )
{
  f = compute_chi2(par);
}

NewCalFullProfileTool::NewCalFullProfileTool( const std::string& type, 
                                      const std::string& name, 
                                      const IInterface* parent)
  : AlgTool(type,name,parent)
{
  // declare base interface for all consecutive concrete classes
  declareInterface<ICalEnergyCorr>(this);
  // Declare the properties that may be set in the job options file
}

StatusCode NewCalFullProfileTool::initialize()
  // This function does following initialization actions:
  //    - extracts geometry constants from xml file using GlastDetSvc
{
  MsgStream log(msgSvc(), name());
  StatusCode sc = StatusCode::SUCCESS;
  log << MSG::DEBUG << "Initializing NewCalFullProfileTool" <<endreq;
  
  //Locate and store a pointer to the data service which allows access to the TDS
  if ((sc = service("EventDataSvc", nm_dataSvc)).isFailure())
    {
      throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
  
  if ((sc = service("GlastDetSvc", nm_detSvc, true)).isFailure())
    { 
      throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }
  
  // find TkrGeometrySvc service
  if ((sc = service("TkrGeometrySvc", nm_tkrGeom, true)).isFailure())
    {
      throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

  IToolSvc* toolSvc = 0;
  if((sc = service("ToolSvc", toolSvc, true)).isFailure())
    {
      throw GaudiException("Service [ToolSvc] not found", name(), sc);
    }

  if((sc = toolSvc->retrieveTool("G4PropagationTool", nm_G4PropTool)).isFailure())
    {
      throw GaudiException("Tool [G4PropagationTool] not found", name(), sc);
      log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
      return StatusCode::FAILURE;
    }

  nm_fsppm = new NewFullShowerProfileParamsManager(0);
  nm_fsddm = new NewFullShowerDevelopmentDescriptionManager(nm_detSvc,14,2.,1.,1.85,0.80,0.1);
  
  // Minuit object
  nm_minuit = new TMinuit(5);
  
  //Sets the function to be minimized
  nm_minuit->SetFCN(fcn);

  nm_saturationadc = 4060;

  ParXYcor[0] = 0.389137;
  ParXYcor[1] = -0.845238;
  ParXYcor[2] = 44.352357;
  ParXYcor[3] = -119.121806;
  ParXYcor[4] = -223.517007;
  ParXYcor[5] = 791.962123;

  return sc;
}

double NewCalFullProfileTool::GetRadiationLengthInTracker(Event::TkrTree* tree)
{
  double tkr_RLn = 0.;

  if(!tree) return tkr_RLn;
  
  Point  x0 = tree->getAxisParams()->getEventPosition();
  Vector vv = tree->getAxisParams()->getEventAxis();
  double costheta = fabs(vv.z());
  
  // copied from CalValsCorrTool
  // Start at the top most point (ie in the top silicon layer) 
  int    topLayer = tree->getHeadNode()->front()->getTreeStartLayer();
  double zTopLyr  = std::max(nm_tkrGeom->getLayerZ(topLayer, 0), nm_tkrGeom->getLayerZ(topLayer, 1));

  // Translate x0 to this z position
  double arcLen   = vv.z() > 0. ? (zTopLyr - x0.z()) / vv.z() : 0.;

  // Translate the tree start position to middle of top silicon layer
  x0 = x0 + arcLen * vv;

  // Now set up and call propagator to get the radiation lengths to calorimeter
  arcLen = vv.z() > 0. ? (x0.z() - nm_tkrGeom->calZTop()) / vv.z() : 0.;
  nm_G4PropTool->setStepStart(x0, -vv);
  nm_G4PropTool->step(arcLen);
  tkr_RLn = nm_G4PropTool->getRadLength(); 

  // Patch for error in KalFitTrack: 1/2 of first radiator left out
  int topPlane = nm_tkrGeom->getPlane(zTopLyr);

  if (nm_tkrGeom->isTopPlaneInLayer(topPlane)) 
    {
      tkr_RLn += 0.5*nm_tkrGeom->getRadLenConv(topLayer) / vv.z();
    }
    
  // add up the rad lens; this could be a local array if you're bothered by the overhead
  //   but hey, compared to the propagator...
  double tkr_radLen_nom = 0.; 
  int layerCount = topLayer;
  for(; layerCount>=0; --layerCount) {
    tkr_radLen_nom += nm_tkrGeom->getRadLenConv(layerCount) 
      + nm_tkrGeom->getRadLenRest(layerCount);
  }
  if(costheta!=0)
    {
      tkr_radLen_nom /= costheta;
      if(tkr_RLn > tkr_radLen_nom * 1.5) {tkr_RLn = tkr_radLen_nom * 1.5;}
      else if(tkr_RLn < tkr_radLen_nom * .5) {tkr_RLn  = tkr_radLen_nom * .5;}
    }
    
  return tkr_RLn;
}

Event::CalCorToolResult* NewCalFullProfileTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrTree* tree)
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
  
  nm_eTotal = cluster->getMomParams().getEnergy()/1000.;

  if( nm_eTotal<0.4) // to get as low as 3 GeV gamma while rejecting most of mips
    {
      lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : nm_eTotal<1GeV -> no energy computation" <<endreq;
      return corResult;
    }
  
  // defines global variable to be used for fcn
  for (i=0;i<8;++i)
    {
      // We are working in GeV
      nm_elayer_dat[i] = (*cluster)[i].getEnergy()/1000.;
      nm_elayer_datsel[i] = nm_elayer_dat[i];
      nm_elayer_datsat[i] = nm_elayer_dat[i];
    }

  // Detect saturation must be called before nm_fsddm->Compute !!!
  DetectSaturation();

  tkr_RLn = 0;
  if(tree!=NULL)
    tkr_RLn = GetRadiationLengthInTracker(tree);

  double CALFIT_fit_energy = 1000.*nm_eTotal;
  double CALFIT_energy_err = 0;
  double CALFIT_fitflag = 0;
  double CALFIT_par0 = 0;
  double CALFIT_par1 = 0;
  double CALFIT_xcor = 0;
  double CALFIT_ycor = 0;
  double CALFIT_alpha = 0;
  double CALFIT_tmax = 0;
  double CALFIT_tkr_RLn = tkr_RLn;
  double CALFIT_lastx0 = 0;
  double CALFIT_cal_eff_RLn = 0;
  double CALFIT_totchisq = 0;
  double CALFIT_chisq = 0;
  double CALFIT_seltotchisq = 0;
  double CALFIT_chisqdist = 0;
  double CALFIT_chidist = 0;
  double CALFIT_parcf = 0;
  double CALFIT_parc = 0;
  double CALFIT_recp0 = 0;
  double CALFIT_recp1 = 0;
  double CALFIT_recp2 = 0;
  double CALFIT_recv0 = 0;
  double CALFIT_recv1 = 0;
  double CALFIT_recv2 = 0;
  double CALFIT_widening = 0;
  int CALFIT_nxtalsel = 0;

  double TKRFIT_fit_energy = 1000.*nm_eTotal;
  double TKRFIT_energy_err = 0;
  double TKRFIT_fitflag = 0;
  double TKRFIT_par0 = 0;
  double TKRFIT_par1 = 0;
  double TKRFIT_xcor = 0;
  double TKRFIT_ycor = 0;
  double TKRFIT_alpha = 0;
  double TKRFIT_tmax = 0;
  double TKRFIT_tkr_RLn = tkr_RLn;
  double TKRFIT_lastx0 = 0;
  double TKRFIT_cal_eff_RLn = 0;
  double TKRFIT_totchisq = 0;
  double TKRFIT_chisq = 0;
  double TKRFIT_seltotchisq = 0;
  double TKRFIT_chisqdist = 0;
  double TKRFIT_chidist = 0;
  double TKRFIT_parcf = 0;
  double TKRFIT_parc = 0;
  double TKRFIT_recp0 = 0;
  double TKRFIT_recp1 = 0;
  double TKRFIT_recp2 = 0;
  double TKRFIT_recv0 = 0;
  double TKRFIT_recv1 = 0;
  double TKRFIT_recv2 = 0;
  double TKRFIT_widening = 0;

  double chi2dist;
  double fitpar[3];
  double result[2];

  nm_wideningfactor = 1;

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

  Vector vaxis(vv[0],vv[1],vv[2]);
  Point corP = cluster->getCorPosition(vaxis);
  pp[0] = corP.x();
  pp[1] = corP.y();
  pp[2] = corP.z();

  SelectCloseCrystals(pp,vv,100);
  CALFIT_nxtalsel = nm_nxtal;

  nm_optselxtalinfit = 0;
  int CALFIT = doProfileFit(pp,vv,0,lm,1);

  CALFIT_fit_energy = 1000.*nm_energy;
  CALFIT_energy_err = 1000.*nm_energy*log(10.)*nm_epar2;
  CALFIT_fitflag = nm_ierflg;
  CALFIT_par0 = nm_par0;
  CALFIT_par1 = nm_par1;
  CALFIT_xcor = nm_xcor;
  CALFIT_ycor = nm_ycor;
  CALFIT_alpha = nm_alpha;
  CALFIT_tmax = nm_tmax;
  CALFIT_tkr_RLn = 0;
  CALFIT_lastx0 = nm_lastx0;
  CALFIT_cal_eff_RLn = nm_totx0cal;
  CALFIT_totchisq = nm_totchisq;
  CALFIT_chisq = nm_chisq;
  CALFIT_parcf = nm_params_contribution_factor;
  CALFIT_parc = nm_params_contribution;
  CALFIT_recp0 = pp[0];
  CALFIT_recp1 = pp[1];
  CALFIT_recp2 = pp[2];
  CALFIT_recv0 = vv[0];
  CALFIT_recv1 = vv[1];
  CALFIT_recv2 = vv[2];
  CALFIT_widening = nm_wideningfactor;
          
  TKRFIT_fit_energy = 1000.*nm_energy;
  TKRFIT_energy_err = 1000.*nm_energy*log(10.)*nm_epar2;
  TKRFIT_fitflag = nm_ierflg;
  TKRFIT_par0 = nm_par0;
  TKRFIT_par1 = nm_par1;
  TKRFIT_xcor = nm_xcor;
  TKRFIT_ycor = nm_ycor;
  TKRFIT_alpha = nm_alpha;
  TKRFIT_tmax = nm_tmax;
  TKRFIT_tkr_RLn = 0;
  TKRFIT_lastx0 = nm_lastx0;
  TKRFIT_cal_eff_RLn = nm_totx0cal;
  TKRFIT_totchisq = nm_totchisq;
  TKRFIT_chisq = nm_chisq;
  TKRFIT_parcf = nm_params_contribution_factor;
  TKRFIT_parc = nm_params_contribution;
  TKRFIT_recp0 = pp[0];
  TKRFIT_recp1 = pp[1];
  TKRFIT_recp2 = pp[2];
  TKRFIT_recv0 = vv[0];
  TKRFIT_recv1 = vv[1];
  TKRFIT_recv2 = vv[2];
  TKRFIT_widening = nm_wideningfactor;

  fitpar[0] = nm_par0;
  fitpar[1] = nm_par1;
  fitpar[2] = nm_par2;
  nm_optselxtalinfit = 1;
  compute_chi2(fitpar);
  CALFIT_seltotchisq = nm_totchisq;
  TKRFIT_seltotchisq = nm_totchisq;
  nm_optselxtalinfit = 0;
  chi2dist = GetChi2Dist(pp,vv,nm_toterrfit,result);
  CALFIT_chisqdist = chi2dist;
  TKRFIT_chisqdist = chi2dist;
  CALFIT_chidist = result[1];
  TKRFIT_chidist = result[1];

  int TKRFIT = 0;

  double mynorm;
  double ppc[3];
  double ppp[3];
  if(tree!=NULL)
    {
      // switch to neutral direction (head of teh track -> cluster centroid)
      // Not using the neutral axis since Tracy improvement of tree axis precision 2012/06/30
      // so no need of the cluster centroid correction
      pp[0] =  tree->getAxisParams()->getEventPosition().x();
      pp[1] =  tree->getAxisParams()->getEventPosition().y();
      pp[2] =  tree->getAxisParams()->getEventPosition().z();
      vv[0] = -tree->getAxisParams()->getEventAxis().x();
      vv[1] = -tree->getAxisParams()->getEventAxis().y();
      vv[2] = -tree->getAxisParams()->getEventAxis().z();
      if(vv[2]>0)
        {
          vv[0] = -vv[0];
          vv[1] = -vv[1];
          vv[2] = -vv[2];
        }

      SelectCloseCrystals(pp,vv,100);

      nm_optselxtalinfit = 0;
      TKRFIT = doProfileFit(pp,vv,tkr_RLn,lm,1);
      //
      TKRFIT_fit_energy = 1000.*nm_energy;
      TKRFIT_energy_err = 1000.*nm_energy*log(10.)*nm_epar2;
      TKRFIT_fitflag = nm_ierflg;
      TKRFIT_par0 = nm_par0;
      TKRFIT_par1 = nm_par1;
      TKRFIT_xcor = nm_xcor;
      TKRFIT_ycor = nm_ycor;
      TKRFIT_alpha = nm_alpha;
      TKRFIT_tmax = nm_tmax;
      TKRFIT_tkr_RLn = tkr_RLn;
      TKRFIT_lastx0 = nm_lastx0;
      TKRFIT_cal_eff_RLn = nm_totx0cal;
      TKRFIT_totchisq = nm_totchisq;
      TKRFIT_chisq = nm_chisq;
      TKRFIT_parcf = nm_params_contribution_factor;
      TKRFIT_parc = nm_params_contribution;
      TKRFIT_recp0 = pp[0];
      TKRFIT_recp1 = pp[1];
      TKRFIT_recp2 = pp[2];
      TKRFIT_recv0 = vv[0];
      TKRFIT_recv1 = vv[1];
      TKRFIT_recv2 = vv[2];
      TKRFIT_widening = nm_wideningfactor;

      fitpar[0] = nm_par0;
      fitpar[1] = nm_par1;
      fitpar[2] = nm_par2;
      nm_optselxtalinfit = 1;
      compute_chi2(fitpar);
      TKRFIT_seltotchisq = nm_totchisq;
      nm_optselxtalinfit = 0;
      chi2dist = GetChi2Dist(pp,vv,nm_toterrfit,result);
      TKRFIT_chisqdist = chi2dist;
      TKRFIT_chidist = result[1];
    }
    
  if(CALFIT==0 && TKRFIT==0)
    {
      lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : No result with cal direction neither with tkr direction." <<endreq;
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
  corResult->insert(Event::CalCorEneValuePair("par0",       TKRFIT_par0));
  corResult->insert(Event::CalCorEneValuePair("par1",       TKRFIT_par1));
  corResult->insert(Event::CalCorEneValuePair("xcor",       TKRFIT_xcor));
  corResult->insert(Event::CalCorEneValuePair("ycor",       TKRFIT_ycor));
  corResult->insert(Event::CalCorEneValuePair("alpha",      TKRFIT_alpha)); // alpha (profile fit parameter 0)
  corResult->insert(Event::CalCorEneValuePair("tmax",       TKRFIT_tmax)); // tmax (profile fit parameter 1)
  corResult->insert(Event::CalCorEneValuePair("tkr_RLn",    TKRFIT_tkr_RLn)); // radiation length in tracker (0 for cal only events)
  corResult->insert(Event::CalCorEneValuePair("lastx0",     TKRFIT_lastx0)); // total radiation length seen by the shower
  corResult->insert(Event::CalCorEneValuePair("cal_eff_RLn",TKRFIT_cal_eff_RLn)); // effective radiation length in cal (i.e in CsI)
  corResult->insert(Event::CalCorEneValuePair("totchisq",   TKRFIT_totchisq)); // total chisquare = usual chisquare + weight * parameters constraint
  corResult->insert(Event::CalCorEneValuePair("seltotchisq",   TKRFIT_seltotchisq));
  corResult->insert(Event::CalCorEneValuePair("chisqdist",   TKRFIT_chisqdist));
  corResult->insert(Event::CalCorEneValuePair("chidist",   TKRFIT_chidist));
  corResult->insert(Event::CalCorEneValuePair("chisq",      TKRFIT_chisq)); // usual chisquare (sunm_layers ( (e-efit)/de )^2
  corResult->insert(Event::CalCorEneValuePair("parcf",      TKRFIT_parcf)); // weighting factor of the parameters constraint
  corResult->insert(Event::CalCorEneValuePair("parc",       TKRFIT_parc)); // parameters constraint contribution to the chisquare
  corResult->insert(Event::CalCorEneValuePair("recp0",      TKRFIT_recp0)); // Point and direction information - kept for the moment for debugging purpose
  corResult->insert(Event::CalCorEneValuePair("recp1",      TKRFIT_recp1));
  corResult->insert(Event::CalCorEneValuePair("recp2",      TKRFIT_recp2));
  corResult->insert(Event::CalCorEneValuePair("recv0",      TKRFIT_recv0));
  corResult->insert(Event::CalCorEneValuePair("recv1",      TKRFIT_recv1));
  corResult->insert(Event::CalCorEneValuePair("recv2",      TKRFIT_recv2));
  corResult->insert(Event::CalCorEneValuePair("widening",   TKRFIT_widening));
  corResult->insert(Event::CalCorEneValuePair("nxtalsat",   nm_nxtal_sat));
  corResult->insert(Event::CalCorEneValuePair("nxtalsel",   nm_nxtal));
  corResult->insert(Event::CalCorEneValuePair("nlayersat",  nm_nsatlayer));

  if(CALFIT>0)
    {
      corResult->insert(Event::CalCorEneValuePair("calfit_fit_energy", CALFIT_fit_energy)); // energy (in MeV)
      corResult->insert(Event::CalCorEneValuePair("calfit_energy_err", CALFIT_energy_err)); // energy error (in MeV)
      corResult->insert(Event::CalCorEneValuePair("calfit_fitflag",    CALFIT_fitflag)); // flag from Minuit
      corResult->insert(Event::CalCorEneValuePair("calfit_par0",       CALFIT_par0));
      corResult->insert(Event::CalCorEneValuePair("calfit_par1",       CALFIT_par1));
      corResult->insert(Event::CalCorEneValuePair("calfit_xcor",       CALFIT_xcor));
      corResult->insert(Event::CalCorEneValuePair("calfit_ycor",       CALFIT_ycor));
      corResult->insert(Event::CalCorEneValuePair("calfit_alpha",      CALFIT_alpha)); // alpha (profile fit parameter 0)
      corResult->insert(Event::CalCorEneValuePair("calfit_tmax",       CALFIT_tmax)); // tmax (profile fit parameter 1)
      corResult->insert(Event::CalCorEneValuePair("calfit_tkr_RLn",    CALFIT_tkr_RLn)); // radiation length in tracker (0 for cal only events)
      corResult->insert(Event::CalCorEneValuePair("calfit_lastx0",     CALFIT_lastx0)); // total radiation length seen by the shower
      corResult->insert(Event::CalCorEneValuePair("calfit_cal_eff_RLn",CALFIT_cal_eff_RLn)); // effective radiation length in cal (i.e in CsI)
      corResult->insert(Event::CalCorEneValuePair("calfit_totchisq",   CALFIT_totchisq)); // total chisquare = usual chisquare + weight * parameters constraint
      corResult->insert(Event::CalCorEneValuePair("calfit_seltotchisq",   CALFIT_seltotchisq));
      corResult->insert(Event::CalCorEneValuePair("calfit_chisqdist",   CALFIT_chisqdist));
      corResult->insert(Event::CalCorEneValuePair("calfit_chidist",   CALFIT_chidist));
      corResult->insert(Event::CalCorEneValuePair("calfit_chisq",      CALFIT_chisq)); // usual chisquare (sunm_layers ( (e-efit)/de )^2
      corResult->insert(Event::CalCorEneValuePair("calfit_parcf",      CALFIT_parcf)); // weighting factor of the parameters constraint
      corResult->insert(Event::CalCorEneValuePair("calfit_parc",       CALFIT_parc)); // parameters constraint contribution to the chisquare
      corResult->insert(Event::CalCorEneValuePair("calfit_recp0",      CALFIT_recp0)); // Point and direction information - kept for the moment for debugging purpose
      corResult->insert(Event::CalCorEneValuePair("calfit_recp1",      CALFIT_recp1));
      corResult->insert(Event::CalCorEneValuePair("calfit_recp2",      CALFIT_recp2));
      corResult->insert(Event::CalCorEneValuePair("calfit_recv0",      CALFIT_recv0));
      corResult->insert(Event::CalCorEneValuePair("calfit_recv1",      CALFIT_recv1));
      corResult->insert(Event::CalCorEneValuePair("calfit_recv2",      CALFIT_recv2));
      corResult->insert(Event::CalCorEneValuePair("calfit_widening",   CALFIT_widening));
      corResult->insert(Event::CalCorEneValuePair("calfit_nxtalsel",   CALFIT_nxtalsel));
    }
  else
    {
      corResult->insert(Event::CalCorEneValuePair("calfit_fit_energy", 0)); // energy (in MeV)
      corResult->insert(Event::CalCorEneValuePair("calfit_energy_err", 0)); // energy error (in MeV)
      corResult->insert(Event::CalCorEneValuePair("calfit_fitflag",    0)); // flag from Minuit
      corResult->insert(Event::CalCorEneValuePair("calfit_par0",       0));
      corResult->insert(Event::CalCorEneValuePair("calfit_par1",       0));
      corResult->insert(Event::CalCorEneValuePair("calfit_xcor",       0));
      corResult->insert(Event::CalCorEneValuePair("calfit_ycor",       0));
      corResult->insert(Event::CalCorEneValuePair("calfit_alpha",      0)); // alpha (profile fit parameter 0)
      corResult->insert(Event::CalCorEneValuePair("calfit_tkr_RLn",    0)); // radiation length in tracker (0 for cal only events)
      corResult->insert(Event::CalCorEneValuePair("calfit_lastx0",     0)); // total radiation length seen by the shower
      corResult->insert(Event::CalCorEneValuePair("calfit_cal_eff_RLn",0)); // effective radiation length in cal (i.e in CsI)
      corResult->insert(Event::CalCorEneValuePair("calfit_totchisq",   -999)); // total chisquare = usual chisquare + weight * parameters constraint
      corResult->insert(Event::CalCorEneValuePair("calfit_seltotchisq",   -999));
      corResult->insert(Event::CalCorEneValuePair("calfit_chisqdist",   -999));
      corResult->insert(Event::CalCorEneValuePair("calfit_chidist",   -999));
      corResult->insert(Event::CalCorEneValuePair("calfit_chisq",      0)); // usual chisquare (sunm_layers ( (e-efit)/de )^2
      corResult->insert(Event::CalCorEneValuePair("calfit_parcf",      0)); // weighting factor of the parameters constraint
      corResult->insert(Event::CalCorEneValuePair("calfit_parc",       0)); // parameters constraint contribution to the chisquare
      corResult->insert(Event::CalCorEneValuePair("calfit_recp0",      0)); // Point and direction information - kept for the moment for debugging purpose
      corResult->insert(Event::CalCorEneValuePair("calfit_recp1",      0));
      corResult->insert(Event::CalCorEneValuePair("calfit_recp2",      0));
      corResult->insert(Event::CalCorEneValuePair("calfit_recv0",      0));
      corResult->insert(Event::CalCorEneValuePair("calfit_recv1",      0));
      corResult->insert(Event::CalCorEneValuePair("calfit_recv2",      0));
      corResult->insert(Event::CalCorEneValuePair("calfit_widening",   0));
      corResult->insert(Event::CalCorEneValuePair("calfit_nxtalsel",   0));
    }

  return corResult;
}

int NewCalFullProfileTool::doProfileFit(double *pp, double *vv, double tkr_RLn, MsgStream lm, int optntrye, double logEradprof, double wf, double *inputfitpar)
  //
  //               This function fits the parameters of shower profile using
  //               the Minuit minimization package and stores the fitted
  //               parameters in the CalCluster object
{

  int i;
  
  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : pp = " <<  pp[0] << " " << pp[1] << " " << pp[2] << " " <<endreq;
  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : vv = " <<  vv[0] << " " << vv[1] << " " << vv[2] << " " <<endreq;
  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : tkr_RLn = " << tkr_RLn <<endreq;

  double wideningfactor = 1.0;
  //  if(wf>0) wideningfactor = wf;

  if(nm_nxtal_sat==0)
    nm_fsddm->m_fsgm->SetOptHE(0);
  else
    nm_fsddm->m_fsgm->SetOptHE(1);
    
//   if(logEradprof==-999)
//     nm_fsddm->m_fsgm->FillRadialProfile(log10(nm_eTotal)+0.2);
//   else
//     nm_fsddm->m_fsgm->FillRadialProfile(logEradprof);

  nm_fsddm->SetWideningFactor(wideningfactor);
  nm_wideningfactor = wideningfactor;

  if(!nm_fsddm->Compute(pp,vv,tkr_RLn))
    {
      lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : Problem during nm_fsddm->Compute. Returning empty corResult." <<endreq;
      return 0;
    }

  if(nm_fsddm->mintotx0cal<0.5) // Make sure that lastx0>0 during the fitting procedure
    {
      lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr :  nm_fsddm->mintotx0cal<0.5. No energy computation. Returning empty corResult." <<endreq;
      return 0;
    }
  
  for(i=0;i<=nm_fsddm->NDevelopment;++i)
    {
      lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : nm_fsddm->FSDDX0[" << i <<"] :"
         << " startx0 = " << nm_fsddm->FSDDX0[i]->startx0
         << " totx0cal = " << nm_fsddm->FSDDX0[i]->totx0cal
         << " lastx0 = " << nm_fsddm->FSDDX0[i]->lastx0
         <<endreq;
    }

  for(i=0;i<8;++i)
    {
      lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : nm_fsddm->meantotx0lay[" << i <<"] :"
         << nm_fsddm->meantotx0lay[i]
         << " meanposx0lay = " << nm_fsddm->meanposx0lay[i]
         <<endreq;
    }

  double arglist[10];
  int ierflg = 0;
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  
  // Set no output because of tuple output
  arglist[0] = -1;
  nm_minuit->mnexcm("SET PRI", arglist ,1,ierflg);
  
  // idem with warnings ( minuit prints warnings
  // when the Hessian matrix is not positive )
  nm_minuit->mnexcm("SET NOW", arglist ,1,ierflg);
  
  arglist[0] = 1;
  nm_minuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  // defines the strategy used by minuit
  // 1 is standard
  // str0
  arglist[0] = 0;
  nm_minuit->mnexcm("SET STR", arglist ,1,ierflg);        
  
  double vstart[3];
  double vstep[3];
  
  double par0,par1,par2;
  double epar0,epar1,epar2;
  
  vstart[0] = 0;
  vstep[0] = 0.5;
  if(inputfitpar && inputfitpar[0]!=-999.)
    {
      vstart[0] = inputfitpar[0];
      nm_minuit->mnparm(0,"a1",vstart[0],vstep[0],-5.,5.,ierflg);
//       nm_minuit->mnparm(0,"a1",vstart[0],vstep[0],0.,0.,ierflg);
//       nm_minuit->FixParameter(0);
    }
  else
    nm_minuit->mnparm(0,"a1",vstart[0],vstep[0],-5.,5.,ierflg);

  vstart[1] = 0;
  vstep[1] = 0.5;
  if(inputfitpar && inputfitpar[1]!=-999.)
    {
      vstart[1] = inputfitpar[1];
      nm_minuit->mnparm(1,"a2",vstart[1],vstep[1],-5.,5.,ierflg);
//       nm_minuit->mnparm(1,"a2",vstart[1],vstep[1],0.,0.,ierflg);
//       nm_minuit->FixParameter(1);
    }
  else
    nm_minuit->mnparm(1,"a2",vstart[1],vstep[1],-5.,5.,ierflg);

  vstart[2] = log10(nm_eTotal)+0.2;
  vstep[2] = 0.1;
  if(inputfitpar && inputfitpar[2]!=-999.)
    {
      vstart[2] = inputfitpar[2];
      //      nm_minuit->mnparm(2,"a3",vstart[2],vstep[2],log10(nm_eTotal)-0.1,log10(nm_eTotal)+2,ierflg);
      nm_minuit->mnparm(2,"a3",vstart[2],vstep[2],0.,0.,ierflg);
      nm_minuit->FixParameter(2);
    }
  else
    nm_minuit->mnparm(2,"a3",vstart[2],vstep[2],log10(nm_eTotal)-0.1,log10(nm_eTotal)+2,ierflg);
  
  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : Initial parameters : " << vstart[0] << " " << vstart[1] << " " << vstart[2] <<endreq;

  // Calls Migrad with 500 iterations maximum
  arglist[0] = 500;
  arglist[1] = 1.;
  //  nm_optpr = 0;

  nm_spy_totchisq = 99999999;
  nm_spy_par[0] = -1;
  nm_spy_par[1] = -1;
  nm_spy_par[2] = -1;
  nm_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  nm_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  if(ierflg==4)
    {
      //      nm_minuit->mnprin(3,amin);
      vstart[0] = nm_spy_par[0];
      vstart[1] = nm_spy_par[1];
      vstart[2] = nm_spy_par[2];
      vstep[0] = 0.5;
      vstep[1] = 0.5;
      vstep[2] = 0.5;
      nm_minuit->mnparm(0, "a1", vstart[0], vstep[0], 0,0,ierflg);
      nm_minuit->mnparm(1, "a2", vstart[1], vstep[1], 0,0,ierflg);
      nm_minuit->mnparm(2, "a3", vstart[2], vstep[2], log10(nm_eTotal)-0.1,log10(nm_eTotal)+2.,ierflg);
      arglist[0] = 500;
      arglist[1] = 1.;
      nm_optpr = 0;
      nm_spy_totchisq = 99999999;
      nm_spy_par[0] = -1;
      nm_spy_par[1] = -1;
      nm_spy_par[2] = -1;
      nm_minuit->mnexcm("MIGRAD", arglist ,2,ierflg);
      nm_minuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    }

  //  nm_minuit->mnprin(3,amin);

  nm_amin = amin;
  nm_ierflg = (double)ierflg;

  nm_minuit->GetParameter(0,par0,epar0);
  nm_minuit->GetParameter(1,par1,epar1);
  nm_minuit->GetParameter(2,par2,epar2);

  //  nm_optpr = 1;
  vstart[0] = par0;
  vstart[1] = par1;
  vstart[2] = par2;
  compute_chi2(vstart);
  //  nm_optpr = 0;
  
  nm_fsppm->Fill(vstart);

  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : ierflg = " << ierflg <<endreq;
  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : totchisq = " << nm_totchisq << " = " << nm_chisq << " + " << nm_params_contribution_factor << " * " << nm_params_contribution <<endreq;
  
  nm_lastx0 = nm_fsddm->CurrentFSDD->lastx0;
  nm_totx0cal = nm_fsddm->CurrentFSDD->totx0cal;
    
  // Clear minuit
  arglist[0] = 1;
  nm_minuit->mnexcm("CLEAR", arglist ,1,ierflg);
      
  lm << MSG::DEBUG << "NewCalFullProfileTool::doEnergyCorr : results : alpha = " << par0 
     << " tmax = " << par1 
     << " energy = " << par2 
     << " fit flag " << ierflg 
     << " totchisq = " << amin <<endreq;

  nm_par0 = par0;
  nm_par1 = par1;
  nm_par2 = par2;
  nm_epar0 = epar0;
  nm_epar1 = epar1;
  nm_epar2 = epar2;
  nm_xcor = nm_fsppm->xcor;
  nm_ycor = nm_fsppm->ycor;
  nm_energy = nm_fsppm->energy;
  nm_alpha = nm_fsppm->alpha;
  nm_beta = nm_fsppm->beta;
  nm_tmax = nm_fsppm->tmax;

  return 1;
}

int NewCalFullProfileTool::DetectSaturation()
{
  nm_fsddm->NXtal = 0;
  nm_nxtal = 0;
  nm_nxtal_sat = 0;

  int i,j,k;

  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        {
          nm_xtal_energy[i][j][k] = 0;
        }

  // Reset nm_fsddm variables
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        nm_fsddm->OffSatu[i][j][k] = -1;
  
  for(j=0;j<8;++j)
    {
      nm_layer_nsat[j] = 0;
    }

  for(j=0;j<FSDD_XTAL_NMAX;++j)
    {
      nm_extal_dat[j] = 0;
      nm_extal_fit[j] = 0;
      nm_xtal_satu[j] = 0;
      nm_xtal_ilay[j] = -1;
    }

  // Get caldigicol
  Event::CalDigiCol *calDigiCol = SmartDataPtr<Event::CalDigiCol>(nm_dataSvc,EventModel::Digi::CalDigiCol);
  if(calDigiCol==NULL) return 1;

  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        {
          nm_saturated[i][j][k] = false;
        }

  for (Event::CalDigiCol::const_iterator digiIter = calDigiCol->begin(); digiIter != calDigiCol->end(); digiIter++)
    {
      const Event::CalDigi calDigi = **digiIter;
      idents::CalXtalId id = calDigi.getPackedId();
      CalUtil::XtalIdx xtalIdx(calDigi.getPackedId());
      for (Event::CalDigi::CalXtalReadoutCol::const_iterator ro =  calDigi.getReadoutCol().begin();ro != calDigi.getReadoutCol().end();ro++)
        {
          float adcP = ro->getAdc(idents::CalXtalId::POS);
          float adcN = ro->getAdc(idents::CalXtalId::NEG);
          int rangepos = ro->getRange(idents::CalXtalId::POS);
          int rangeneg = ro->getRange(idents::CalXtalId::NEG);
          if( (rangepos==3 && adcP>=nm_saturationadc) || (rangeneg==3 && adcN>=nm_saturationadc))
            {
              nm_saturated[(int)id.getTower()][(int)id.getLayer()][(int)id.getColumn()] = true;
            }
        }
    }

  Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(nm_dataSvc, EventModel::CalRecon::CalXtalRecCol); 
  if(calXtalRecCol==NULL)  return 1;

  nm_nxtal = 0;
  for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
    {
      Event::CalXtalRecData* xTalData = *xTalIter;
      int itow=xTalData->getPackedId().getTower();
      int ilay=xTalData->getPackedId().getLayer();
      int icol=xTalData->getPackedId().getColumn();
      nm_xtal_energy[itow][ilay][icol] = xTalData->getEnergy()/1000;

      if(nm_saturated[itow][ilay][icol])
        {
          //
          if(nm_saturated[itow][ilay][icol]) nm_layer_nsat[ilay] += 1;
          //
          nm_elayer_dat[ilay] -= xTalData->getEnergy()/1000;
          nm_elayer_datsat[ilay] -= xTalData->getEnergy()/1000;
          nm_elayer_datsel[ilay] -= xTalData->getEnergy()/1000;
          //
          nm_extal_dat[nm_nxtal] = xTalData->getEnergy()/1000;
          if(nm_saturated[itow][ilay][icol])
            {
              nm_xtal_satu[nm_nxtal] = 1;
            }
          nm_xtal_dist[nm_nxtal] = 0;
          nm_xtal_itow[nm_nxtal] = itow;
          nm_xtal_ilay[nm_nxtal] = ilay;
          nm_xtal_icol[nm_nxtal] = icol;
          nm_fsddm->OffSatu[itow][ilay][icol] = nm_nxtal;
          ++nm_nxtal;
          nm_fsddm->NXtal = nm_nxtal;
        }
      if(nm_nxtal==FSDD_XTAL_NMAX) break;
    }

  nm_nxtal_sat = nm_nxtal;

  nm_nsatlayer = 0;
  for(j=0;j<8;++j)
    {
      if(nm_layer_nsat[j]>0)
        {
          nm_nsatlayer += 1;
        }
    }

  return 0;
}

double NewCalFullProfileTool::GetCrystalDistanceToTrajectory(Point p0, Vector v0, Point p1, Vector v1, double xtalhalflength)
{
  // p0,v0 = traj
  // p1,v1 = xtal, with p1 = center of xtal and v1 = long axis

  Vector p0p1 = p1-p0;
  double lambda = 1-(v0*v1)*(v0*v1);
  double mydist2 = 0;
  if(lambda==0) // xtal axis and main axis are parallel
    {
      mydist2 = p0p1*p0p1 - (p0p1*v1)*(p0p1*v1);
      return sqrt(mydist2);
    }
  
  double lambda0 = ((p0p1*v0)-(p0p1*v1)*(v0*v1))/lambda;
  Point p00 = p0+lambda0*v0;
  double lambda1 = lambda0*(v0*v1)-(p0p1*v1);
  if(lambda1<-xtalhalflength)
    lambda1 = -xtalhalflength;
  else if(lambda1>xtalhalflength)
    lambda1 = xtalhalflength;
  Point p11 = p1+lambda1*v1;
  p0p1 = p11-p00;
  mydist2 = p0p1*p0p1;

  return sqrt(mydist2);
}

int NewCalFullProfileTool::SelectCloseCrystals(double *pptraj, double *vvtraj, double maxdistance)
{
  int i,i0,i1,j,k;

  double pp1[3];
  double vv1x[3] = {1,0,0};
  double vv1y[3] = {0,1,0};
  double *vv1;
  double mydist;
  double xtow,ytow;

  Vector v0 = Vector(vvtraj[0],vvtraj[1],vvtraj[2]);
  Point p0 = Point(pptraj[0],pptraj[1],pptraj[2]);

  int icount2;

  // reset the nb of selected xtals to the nb of saturated xtals;
  nm_nxtal = nm_nxtal_sat;
  nm_fsddm->NXtal = nm_nxtal;
  for(j=0;j<8;++j) nm_elayer_datsel[j] = nm_elayer_datsat[j];
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        if(nm_fsddm->OffSatu[i][j][k]>=nm_nxtal_sat) nm_fsddm->OffSatu[i][j][k] = -1;

  for(j=0;j<8;++j)
    {
      if(j%2==0)
        vv1 = vv1x;
      else
        vv1 = vv1y;
      pp1[2] = nm_fsddm->m_fsgm->FSGM_calZTop-nm_fsddm->m_fsgm->FSGM_cellVertPitch*(0.5+(double)j);
      //
      for(i0=0;i0<4;++i0)
        for(i1=0;i1<4;++i1)
          {
            i = 4*i1+i0;
            xtow = nm_fsddm->m_fsgm->FSGM_towerPitch*(-2.0+0.5+(double)i0);
            ytow = nm_fsddm->m_fsgm->FSGM_towerPitch*(-2.0+0.5+(double)i1);
            for(k=0;k<12;++k)
              {
                if(j%2==0)
                  {
                    pp1[0] = xtow;
                    pp1[1] = ytow + nm_fsddm->m_fsgm->FSGM_cellHorPitch*(-6.0+0.5+(double)k);
                  }
                else
                  {
                    pp1[0] = xtow + nm_fsddm->m_fsgm->FSGM_cellHorPitch*(-6.0+0.5+(double)k);
                    pp1[1] = ytow;
                  }
                //
                Vector v1 = Vector(vv1[0],vv1[1],vv1[2]);
                Point p1 = Point(pp1[0],pp1[1],pp1[2]);
                mydist = GetCrystalDistanceToTrajectory(p0,v0,p1,v1,nm_fsddm->m_fsgm->FSGM_CsILength/2);
                //
                nm_xtal_doca[i][j][k] = mydist;
                //
                if(nm_saturated[i][j][k]) // already in xtal list
                  {
                    icount2 = nm_fsddm->OffSatu[i][j][k];
                    nm_xtal_dist[icount2] = mydist;
                    continue;
                  }
                //
                if(mydist>maxdistance) continue;
                //
                if(nm_nxtal<FSDD_XTAL_NMAX)
                  {
                    nm_xtal_itow[nm_nxtal] = i;
                    nm_xtal_ilay[nm_nxtal] = j;
                    nm_xtal_icol[nm_nxtal] = k;
                    nm_fsddm->OffSatu[i][j][k] = nm_nxtal;
                    nm_extal_dat[nm_nxtal] = nm_xtal_energy[i][j][k];
                    nm_elayer_datsel[j] -= nm_xtal_energy[i][j][k];
                    nm_extal_fit[nm_nxtal] = 0.0;
                    nm_xtal_satu[nm_nxtal] = 0;
                    ++nm_nxtal;
                    nm_fsddm->NXtal = nm_nxtal;
                  }
                else
                  break;
              }
          }
    }
  
  return 0;
}

double NewCalFullProfileTool::GetChi2Dist(double *pptraj, double *vvtraj, double fiterror, double *result)
{
  if(result)
    {
      result[0] = 0;
      result[1] = 0;
    }
  double chi2dist = 0;
  double chidist = 0;
  int i;
  
  Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(nm_dataSvc, EventModel::CalRecon::CalXtalRecCol); 
  if(calXtalRecCol==NULL)  return chi2dist;

  double xyz[3];
  double xyzp[3];
  double lambda,mydist;

  double xtalfitenergy = 0;
  double xtaldatenergy = 0;

  int icount;
  int xtalused[FSDD_XTAL_NMAX];
  for(i=0;i<nm_nxtal;++i) xtalused[i] = 0;

  for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
    {
      Event::CalXtalRecData* xTalData = *xTalIter;
      int itow=xTalData->getPackedId().getTower();
      int ilay=xTalData->getPackedId().getLayer();
      int icol=xTalData->getPackedId().getColumn();
      nm_xtal_energy[itow][ilay][icol] = xTalData->getEnergy()/1000;
      xyz[0] = xTalData->getPosition().x();
      xyz[1] = xTalData->getPosition().y();
      xyz[2] = xTalData->getPosition().z();
      //
      xtaldatenergy = nm_xtal_energy[itow][ilay][icol];
      xtalfitenergy = 0;
      icount = nm_fsddm->OffSatu[itow][ilay][icol];
      if(icount>-1)
        {
          xtalfitenergy = nm_extal_fit[icount];
          if(icount<nm_nxtal_sat && xtalfitenergy>xtaldatenergy) xtalfitenergy = xtaldatenergy;
          xtalused[icount] = 1;
        }
      lambda = 0;
      for(i=0;i<3;++i) lambda += (xyz[i]-pptraj[i])*vvtraj[i];
      for(i=0;i<3;++i) xyzp[i] = pptraj[i]+lambda*vvtraj[i];
      mydist = 0;
      for(i=0;i<3;++i) mydist += (xyz[i]-xyzp[i])*(xyz[i]-xyzp[i]);
      mydist = sqrt(mydist);
      //
      chi2dist += (xtaldatenergy-xtalfitenergy)*(xtaldatenergy-xtalfitenergy)/fiterror/fiterror*mydist*mydist;
      chidist += fabs((xtaldatenergy-xtalfitenergy))/fiterror*mydist*mydist;
    }

  for(i=0;i<nm_nxtal;++i)
    {
      if(xtalused[i]) continue;
      xtaldatenergy = nm_xtal_energy[nm_xtal_itow[i]][nm_xtal_ilay[i]][nm_xtal_icol[i]];
      xtalfitenergy = nm_extal_fit[i];
      mydist = nm_xtal_doca[nm_xtal_itow[i]][nm_xtal_ilay[i]][nm_xtal_icol[i]];
      chi2dist += (xtaldatenergy-xtalfitenergy)*(xtaldatenergy-xtalfitenergy)/fiterror/fiterror*mydist*mydist;
      chidist += fabs((xtaldatenergy-xtalfitenergy))/fiterror*mydist*mydist;
    }

  if(result)
    {
      result[0] = chi2dist;
      result[1] = chidist;
    }

  return chi2dist;
}

StatusCode NewCalFullProfileTool::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;
  
  // delete Minuit object
  delete nm_minuit;
  
  delete nm_fsddm;
  delete nm_fsppm; 
  
  return sc;
}

