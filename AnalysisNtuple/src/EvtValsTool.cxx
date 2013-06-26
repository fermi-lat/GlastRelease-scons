
/** @file EvtValsTool.cxx
@brief Calculates the "Event" analysis variables from the other ntuple variables
@author Bill Atwood, Leon Rochester

$Header$
*/

#include "ValBase.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "TkrUtil/ITkrGeometrySvc.h"
#include "geometry/Ray.h"

#include "CalUtil/IUBinterpolateTool.h"
#include "IPsfTool.h"

#include "facilities/Util.h"                // for expandEnvVar
#include "CLHEP/Vector/Rotation.h"

#include <algorithm>
/** @class EvtValsTool
@brief Calculates Event Values from the other ntuple variables
*/

// copied from CalValsCorTool
namespace {
  /// A local erf function: alg. from Receipts in C
  double cvct_erf_cal(double x) {
    double t = 1./(1.+.47047*x);
    double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
    return results;
  }
  double cvct_cal_trans(double x) {
    if(x < 0) return (.5*(1. + cvct_erf_cal(-x)));
    else      return (.5*(1. - cvct_erf_cal( x)));
  }
}

class EvtValsTool : public ValBase
{
public:

    EvtValsTool( const std::string& type, 
        const std::string& name, 
        const IInterface* parent);

    virtual ~EvtValsTool() { }

    StatusCode initialize();

    StatusCode calculate();

    void zeroVals();

    void fillHeaderInfo();

  double GetEnergyUB2Correction(int method, int tkr1firstlayer, double tkr1zdir, double energy);

  void fillPsfInfo(double energy, double theta, bool isFront, double cl_level);

  // Copied from CalValsCorTool
  double aveRadLens(Point x0, Vector t0, double radius, int numSamples, double calcentroidposz);

private:

  // some local constants
  double m_towerPitch;
  double m_calZTop;
  double m_calZBot;
  double m_xLo;
  double m_xHi;
  double m_yLo;
  double m_yHi;

  // Internal Variables for aveRadLens calculation (as in CalValsCorTool)
  double             m_radLen_CsI, m_rms_RL_CsI; //Rad. lengths along axis in CsI Xtals (and rms)
  double             m_radLen_Stuff, m_rms_RL_Stuff; // Rad. lengths of other material along axis (and rms)
  double             m_radLen_Cntr, m_rms_RL_Cntr; // Rad. lengths along axis to centroid (and rms)
  double             m_radLen_CntrStuff, m_rms_RL_CntrStuff; // Rad. length of other material (not CsI) to centroid
  double             m_arcLen_CsI;      // Path length along shower axis in CsI Xtals
  double             m_arcLen_Stuff;    // Path length along shower axis for material other then CsI
  double             m_arcLen_Cntr;     // Path length along shower axis to energy centroid

    //Global Items
    unsigned int EvtRun;
    unsigned int EvtEventId;
    unsigned long long EvtEventId64;
    double EvtElapsedTime;
    double EvtLiveTime;

    float EvtEnergyCorr;
    float EvtEnergyCorrUB;
    float EvtEnergyRaw;
    float EvtDeltaEoE;
    float EvtCalEdgeAngle;
    float EvtTkrEdgeAngle;
    float EvtLogEnergy;
    float EvtTkr1EFrac;
    float EvtVtxKin;
    float EvtVtxEAngle;
    float EvtTkrComptonRatio;
    float EvtETkrComptonRatio;
    float EvtPSFModel; 

    float EvtETkr1Chisq;
    float EvtETkr1FirstChisq;
    float EvtETkr1Qual;
    float EvtTkr1PSFMdRat;

    float EvtECalXtalRatio;
    float EvtECalXtalTrunc;
    float EvtECalTrackDoca;
    //float EvtECalTrackSep;
        float EvtECalTransRms;
        float EvtECalLongRms;
        float EvtECalLRmsAsym;
        float EvtECalTrackAngle;
    float EvtEVtxAngle;
    float EvtEVtxDoca;
    float EvtPsf68;
    unsigned int EvtEventFlags;

  float NewEvtEnergyCorr;
  float NewEvtEnergyCorrUB2;
  float EvtEnergyCorrUB2;
  float CalNewCfpEnergyUB2;

  float EvtJointEnergy;
  float EvtJointLogEnergy;
  float EvtJointWeight;
  
  float EvtCalCsIRLn;

    //test
    //char  EvtEvtNum[20];

    //IValsTool* m_pMcTool;
    IValsTool* m_pGltTool;
    IValsTool* m_pTkrHitTool;
    IValsTool* m_pTkrTool;
    IValsTool* m_pVtxTool;
    IValsTool* m_pCalTool;
    IValsTool* m_pAcdTool;

    ITkrGeometrySvc* m_tkrGeom;
    IGlastDetSvc* m_detSvc; 
  IPropagator* m_G4PropTool; 

        // These are copies of the new DC2 vars from IM.   Left here for reference
// EvtECalXtalRatio        continuous        CalXtalRatio/(2.99-1.19*LogEnergyOpt  + .122*LogEnergyOpt^2)/(.749-.355*Tkr1ZDir)
// EvtECalTrackDoca        continuous        CalTrackDoca/(272-140.5*min(EvtLogEnergy, 3.8) + 18.7*min(EvtLogEnergy, 3.8)^2)/(3.08+2.67*Tkr1ZDir)
// EvtECalXtalsTrunc        continuous        CalXtalsTrunc/(max(1.,(-85.4 + 65.2*LogEnergyOpt - 10.4*LogEnergyOpt^2)) + min(12*(LogEnergyOpt-3.3)^2*(max(3.3,LogEnergyOpt)-3.3)/(LogEnergyOpt-3.3),14.))/(1.2+.224*Tkr1ZDir)
// EvtECalTransRms        continuous        CalTransRms/(12.5+ 22.*exp(-(EvtLogEnergy-2.3)^2/.8))/(1.34+.55*Tkr1ZDir)
// EvtECalLongRms        continuous        CalLongRms/(ifelse(EvtLogEnergy < 3.8, (96.1 - 39*EvtLogEnergy + 5.6*EvtLogEnergy^2),  28.*(.847-.0134*EvtLogEnergy+.015*EvtLogEnergy^2)))/(1.06+.0867*Tkr1ZDir)
// EvtECalLRmsAsym        continuous        CalLRmsAsym/(.012+.06*exp(-(EvtLogEnergy-2.0)^2/1.3))/(1.01-.0718*Tkr1ZDir)
// EvtECalTrackAngle        continuous        CalTrackAngle/(2.3-1.11*min(EvtLogEnergy,4.0)+.138*min(EvtLogEnergy,4.0)^2)/(1.43+.612*Tkr1ZDir)

// EvtETkr1Chisq        continuous        Tkr1Chisq/(6.49 -22.1*Tkr1ZDir -22.8*Tkr1ZDir^2)
// EvtETkr1FirstChisq        continuous        Tkr1FirstChisq/(4.34-.708*EvtLogEnergy+.19*EvtLogEnergy^2)/(.751-1.74*Tkr1ZDir-1.83*Tkr1ZDir^2)
// EvtETkr1Qual        continuous        Tkr1Qual/(56.6 +6.78*Tkr1ZDir + 8.23*Tkr1ZDir^2)
// EvtETkrComptonRatio        continuous        EvtTkrComptonRatio / (-9.34+7.22*EvtLogEnergy - .672*EvtLogEnergy^2)/(-.509 - 3.65*Tkr1ZDir - 2.03*Tkr1ZDir^2)

// EvtVtxEAngle        continuous        VtxAngle*sqrt(EvtEnergySumOpt)/(4.24 -1.98*EvtLogEnergy + .269*EvtLogEnergy^2)/(1.95+2.36*Tkr1ZDir+1.3*Tkr1ZDir^2)

  IUBinterpolateTool* m_ubInterpolateTool;
  IPsfTool* m_pPsfTool;
  std::string m_psfVersion;
  std::string m_psfPath;

  double UB2logemin[3];
  double UB2logemax[3];
  int UB2zdirn;
  double UB2zdir[9];
  double UB2zdirbinwidth;
  double UB2val[9];
  double UB2par[3][18][8][6];
};


// Static factory for instantiation of algtool objects
//static ToolFactory<EvtValsTool> s_factory;
//const IToolFactory& EvtValsToolFactory = s_factory;
DECLARE_TOOL_FACTORY(EvtValsTool);

// Standard Constructor
EvtValsTool::EvtValsTool(const std::string& type, 
                         const std::string& name, 
                         const IInterface* parent)
                         : ValBase( type, name, parent )
{    
    // Declare additional interface
    declareInterface<IValsTool>(this);

    declareProperty("psfVersion",m_psfVersion="P7SOURCE_V6MC");
    declareProperty("psfPath",m_psfPath="$(ANALYSISNTUPLEDATAPATH)");
}

/** @page anatup_vars 
    @section evtvalstool EvtValsTool Variables

These are calculated from combinations of the variables from different tools. 

NOTE
- All EvtEXxx variables are previous variables compensated for energy and angle 

<table>
<tr><th> Variable <th> Type <th> Description
<tr><td> EvtRun 
<td>U<td>   Run number, copied from the event header NEW: replaces Run in the merit ntuple  
<tr><td> EvtEventId 
<td>U<td>   Sequence number of event in the run
<tr><td> EvtEventId64
<td>UL<td>   Sequence number of event in the run, 64-bit unsigned version
<tr><td> EvtElapsedTime 
<td>D<td>   Elapsed time in seconds since t0 (for DC1: 18-July-2005, 
            for the future: mission start 1-Jan-2001)
<tr><td> EvtLivetime 
<td>D<td>  (s) Cumulative live time, from start of run, or mission
<tr><td> EvtEnergyCorr 
<td>F<td>   Event energy formed by adding the corrected tracker energy 
            (TkrEnergyCorr) to the layer-by-layer corrected cal. energy CalEnergyCorr. 
<tr><td> EvtEnergyRaw 
<td>F<td>   TkrEnergy + CalEnergyRaw
<tr><td> EvtCalEdgeAngle 
<td>F<td>   Obsolete; replaced by CalTwrGap 
<tr><td> EvtTkrEdgeAngle 
<td>F<td>   Obsolete; replaced by Tkr1TwrGap 
<tr><td> EvtLogEnergy 
<td>F<td>   log10 of EvtEnergyCorr, pegged between log10(20) and log10(50,000). Was EvtLogESum 
<tr><td> EvtTkr1EFrac 
<td>F<td>   Tkr1ConE/EvtEnergyCorr, roughly, fraction of energy carried by best track 
<tr><td> EvtVtxKin 
<td>F<td>   The vertex opening angle compenstated for the energy split between the tracks. 
<tr><td> EvtVtxEAngle 
<td>F<td>   VtxAngle*EvtEnergyCorr.  Should be approx. constant. 
            However an empirical compensation is provided below (see EvtEVtxAngle) 
<tr><td> EvtTkrComptonRatio 
<td>F<td>   Ratio of TkrTotalHits to twice the number of layers from the head 
         of the best track to the bottom of the TKR 
<tr><td> EvtETkrComptonRatio 
<td>F<td>   EvtTkrComptonRatio, flattened in energy and cos(theta). Was EvtTkrEComptonRatio 
<tr><td> EvtPSFModel 
<td>F<td>   PSF expected from simple model; depends only on energy. 
<tr><td> EvtETkr1Chisq 
<td>F<td>   Tkr1Chisq, compensated for energy and angle. 
<tr><td> EvtETkr1FirstChisq 
<td>F<td>   Tkr1FirstChisq, compensated for energy and angle 
<tr><td> EvtETkr1Qual 
<td>F<td>   Tkr1Qual, compensated for energy and angle 
<tr><td> EvtTkr1PSFMdRat 
<td>F<td>   Ratio of errors from covariance matrix to EvtPSFModel 
<tr><td> EvtECalTransRms 
<td>F<td>   CalTransRms, compensated for energy and angle 
<tr><td> EvtECalLongRms 
<td>F<td>   CalLongRms, compensated for energy and angle 
<tr><td> EvtECalLRmsAsym 
<td>F<td>   CalLRmsAsym, compensated for energy and angle 
<tr><td> EvtECalXtalRatio 
<td>F<td>   CalXtalRatio, compensated for energy and angle 
<tr><td> EvtECalXtalTrunc 
<td>F<td>   CalXtalsTrunc, compensated for energy and angle 
<tr><td> EvtECalTrackDoca 
<td>F<td>   CalTrackDoca, compensated for energy and angle 
<tr><td> EvtECalTrackSep 
<td>F<td>   REMOVED! CalTrackSep, compensated for energy and angle 
<tr><td> EvtEVtxAngle 
<td>F<td>   EvtVtxEAngle, compensated for energy and angle 
<tr><td> EvtEVtxDoca 
<td>F<td>   VtxDOCA, compensated for energy and angle 
<tr><td> EvtEventFlags
<td>U<td> Gleam Event Flags, zero denotes no error bits.  see enums/EventFlags.h
          for a definition of the error bits

<tr><td> EvtJointEnergy 
<td>F<td>   DOCUMENTATION NEEDED
<tr><td> EvtJointLogEnergy 
<td>F<td>   DOCUMENTATION NEEDED
<tr><td> EvtJointWeight 
<td>F<td>   DOCUMENTATION NEEDED
</table>
*/

StatusCode EvtValsTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;
    StatusCode fail = StatusCode::FAILURE;

    MsgStream log(msgSvc(), name());

    if( ValBase::initialize().isFailure()) return StatusCode::FAILURE;

    // get the services

    // try using one tool from another

    if(service( "TkrGeometrySvc", m_tkrGeom, true ).isFailure()) {
        log << MSG::ERROR << "Could not find TkrGeometrySvc" << endreq;
        return fail;
    }

    m_towerPitch = m_tkrGeom->towerPitch();
    m_calZTop = m_tkrGeom->calZTop();
    m_calZBot = m_tkrGeom->calZBot();
    m_xLo = m_tkrGeom->getLATLimit(0, LOW);
    m_xHi = m_tkrGeom->getLATLimit(0, HIGH);
    m_yLo = m_tkrGeom->getLATLimit(1, LOW);
    m_yHi = m_tkrGeom->getLATLimit(1, HIGH);
    // Ph. Bruel : hardcoded modification in order to take into account the CU geometry (tower 1 without tracker)
    if(m_xLo==0) m_xLo = -m_towerPitch;

    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
      log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
      return StatusCode::FAILURE;
    }

    IToolSvc* pToolSvc = 0; 
    sc = service("ToolSvc", pToolSvc, true);
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Can't find ToolSvc, will quit now" << endreq;
        return StatusCode::FAILURE;
    }

    if(!pToolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
      log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
      return StatusCode::FAILURE;
    }

    //m_pMcTool = 0;
    //sc = pToolSvc->retrieveTool("McValsTool", m_pMcTool);
    //if( sc.isFailure() ) {
    //    log << MSG::INFO << "Unable to find tool: " "McValsTool" << endreq;
    //    log << "Will carry on anyway, EvtDeltaEoE will not be calculated" << endreq;
    //}
        
    /*    
    m_pGltTool = 0;
    sc = pToolSvc->retrieveTool("GltValsTool", m_pGltTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "GltValsTool" << endreq;
        return sc;
    }
    m_pTkrHitTool = 0;
    sc = pToolSvc->retrieveTool("TkrHitValsTool", m_pTkrHitTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "TkrHitValsTool" << endreq;
        return sc;
    }
    */

    m_pTkrTool = 0;
    sc = pToolSvc->retrieveTool("TkrValsTool", m_pTkrTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "TkrValsTool" << endreq;
        return sc;
    }

    m_pVtxTool = 0;
    sc = pToolSvc->retrieveTool("VtxValsTool", m_pVtxTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "VtxValsTool" << endreq;
        return sc;
    }

    m_pCalTool = 0;
    sc = pToolSvc->retrieveTool("CalValsTool", m_pCalTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "CalValsTool" << endreq;
        return sc;
    }
    
    m_ubInterpolateTool = 0;
    sc = pToolSvc->retrieveTool("UBinterpolateTool", m_ubInterpolateTool);
    if(sc.isFailure() ) {
      log << MSG::ERROR << "Couldn't retrieve UBinterpolateTool" << endreq;
      return StatusCode::FAILURE;
    }

    /*
    m_pAcdTool = 0;
    sc = pToolSvc->retrieveTool("AcdValsTool", m_pAcdTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "AcdValsTool" << endreq;
        return sc;
    }
    */

    m_pPsfTool = 0;
    sc = pToolSvc->retrieveTool("PsfValsTool", m_pPsfTool);
    if( sc.isFailure() ) {
        log << MSG::ERROR << "Unable to find tool: " "PsfValsTool" << endreq;
        return sc;
    } else {
      facilities::Util::expandEnvVar(&m_psfPath);
      sc = m_pPsfTool->loadPsf(m_psfVersion, m_psfPath);
      if(sc.isFailure()) return sc;
    }


    // Ph.Bruel: UB2. Introduced when unbias NewEvtEnergyCorr and EvtEnergyCorr (after correction of removal leakage correction)
    // UB2 loge boudaries set to the center of the lowest and highest energy bins of the unbias analyses
    UB2logemin[0] = 1.0; 
    UB2logemax[0] = 4.0;
    UB2logemin[1] = 1.0; 
    UB2logemax[1] = 4.0;
    UB2logemin[2] = 2.75; 
    UB2logemax[2] = 5.00;
    UB2zdirn = 8;
    UB2zdirbinwidth = 0.1;
    for(int i=0;i<=UB2zdirn;++i) UB2zdir[i] = -0.95+UB2zdirbinwidth*i;
    
    // First the UB2 parameters for the first costheta bins (from [-1.0,-0.9] to [-0.6,-0.5]). Then the 3 other costhea bins.
   UB2par[0][2][0][0]=3.640238e-01; UB2par[0][2][0][1]=1.200958e-01; UB2par[0][2][0][2]=4.862053e-01; UB2par[0][2][0][3]=-2.893312e-01; UB2par[0][2][0][4]=5.719196e-02; UB2par[0][2][0][5]=-3.810883e-03;
   UB2par[0][2][1][0]=8.186851e-01; UB2par[0][2][1][1]=-1.366298e+00; UB2par[0][2][1][2]=1.983881e+00; UB2par[0][2][1][3]=-9.620639e-01; UB2par[0][2][1][4]=1.993735e-01; UB2par[0][2][1][5]=-1.525329e-02;
   UB2par[0][2][2][0]=1.869024e+00; UB2par[0][2][2][1]=-4.218828e+00; UB2par[0][2][2][2]=4.612385e+00; UB2par[0][2][2][3]=-2.071700e+00; UB2par[0][2][2][4]=4.213301e-01; UB2par[0][2][2][5]=-3.227414e-02;
   UB2par[0][2][3][0]=2.822325e+00; UB2par[0][2][3][1]=-6.743960e+00; UB2par[0][2][3][2]=6.850708e+00; UB2par[0][2][3][3]=-2.981045e+00; UB2par[0][2][3][4]=5.970251e-01; UB2par[0][2][3][5]=-4.533042e-02;
   UB2par[0][2][4][0]=2.971755e+00; UB2par[0][2][4][1]=-7.116596e+00; UB2par[0][2][4][2]=7.044642e+00; UB2par[0][2][4][3]=-3.006104e+00; UB2par[0][2][4][4]=5.937654e-01; UB2par[0][2][4][5]=-4.461176e-02;
   UB2par[1][2][0][0]=-4.372376e+00; UB2par[1][2][0][1]=1.130019e+01; UB2par[1][2][0][2]=-9.041805e+00; UB2par[1][2][0][3]=3.503086e+00; UB2par[1][2][0][4]=-6.624163e-01; UB2par[1][2][0][5]=4.884319e-02;
   UB2par[1][2][1][0]=-1.214873e+00; UB2par[1][2][1][1]=3.522340e+00; UB2par[1][2][1][2]=-2.089379e+00; UB2par[1][2][1][3]=5.981506e-01; UB2par[1][2][1][4]=-8.471794e-02; UB2par[1][2][1][5]=4.736130e-03;
   UB2par[1][2][2][0]=4.366510e-01; UB2par[1][2][2][1]=-8.230063e-01; UB2par[1][2][2][2]=1.863608e+00; UB2par[1][2][2][3]=-1.054956e+00; UB2par[1][2][2][4]=2.424353e-01; UB2par[1][2][2][5]=-2.006198e-02;
   UB2par[1][2][3][0]=2.768427e+00; UB2par[1][2][3][1]=-6.051556e+00; UB2par[1][2][3][2]=6.119176e+00; UB2par[1][2][3][3]=-2.684138e+00; UB2par[1][2][3][4]=5.416632e-01; UB2par[1][2][3][5]=-4.132379e-02;
   UB2par[1][2][4][0]=2.892101e+00; UB2par[1][2][4][1]=-6.255216e+00; UB2par[1][2][4][2]=6.108704e+00; UB2par[1][2][4][3]=-2.609615e+00; UB2par[1][2][4][4]=5.171437e-01; UB2par[1][2][4][5]=-3.896411e-02;

   UB2par[0][3][0][0]=9.952621e-01; UB2par[0][3][0][1]=-2.301380e+00; UB2par[0][3][0][2]=2.975552e+00; UB2par[0][3][0][3]=-1.388348e+00; UB2par[0][3][0][4]=2.797526e-01; UB2par[0][3][0][5]=-2.079522e-02;
   UB2par[0][3][1][0]=1.777296e+00; UB2par[0][3][1][1]=-4.347971e+00; UB2par[0][3][1][2]=4.816375e+00; UB2par[0][3][1][3]=-2.161685e+00; UB2par[0][3][1][4]=4.366672e-01; UB2par[0][3][1][5]=-3.313543e-02;
   UB2par[0][3][2][0]=2.555572e+00; UB2par[0][3][2][1]=-6.153082e+00; UB2par[0][3][2][2]=6.267545e+00; UB2par[0][3][2][3]=-2.704138e+00; UB2par[0][3][2][4]=5.340852e-01; UB2par[0][3][2][5]=-3.989160e-02;
   UB2par[0][3][3][0]=2.473942e+00; UB2par[0][3][3][1]=-5.768634e+00; UB2par[0][3][3][2]=5.718165e+00; UB2par[0][3][3][3]=-2.411214e+00; UB2par[0][3][3][4]=4.691106e-01; UB2par[0][3][3][5]=-3.474188e-02;
   UB2par[0][3][4][0]=2.371645e+00; UB2par[0][3][4][1]=-5.305159e+00; UB2par[0][3][4][2]=5.130992e+00; UB2par[0][3][4][3]=-2.123750e+00; UB2par[0][3][4][4]=4.086719e-01; UB2par[0][3][4][5]=-3.007832e-02;
   UB2par[1][3][0][0]=-3.176406e+00; UB2par[1][3][0][1]=6.915739e+00; UB2par[1][3][0][2]=-4.313811e+00; UB2par[1][3][0][3]=1.308579e+00; UB2par[1][3][0][4]=-1.982154e-01; UB2par[1][3][0][5]=1.205257e-02;
   UB2par[1][3][1][0]=-2.335003e+00; UB2par[1][3][1][1]=4.616989e+00; UB2par[1][3][1][2]=-2.258576e+00; UB2par[1][3][1][3]=4.702864e-01; UB2par[1][3][1][4]=-3.518550e-02; UB2par[1][3][1][5]=-1.837783e-04;
   UB2par[1][3][2][0]=-7.006836e-02; UB2par[1][3][2][1]=-4.658618e-01; UB2par[1][3][2][2]=1.909510e+00; UB2par[1][3][2][3]=-1.142952e+00; UB2par[1][3][2][4]=2.646935e-01; UB2par[1][3][2][5]=-2.173145e-02;
   UB2par[1][3][3][0]=4.177027e-01; UB2par[1][3][3][1]=-1.482214e+00; UB2par[1][3][3][2]=2.609446e+00; UB2par[1][3][3][3]=-1.365846e+00; UB2par[1][3][3][4]=3.001401e-01; UB2par[1][3][3][5]=-2.405324e-02;
   UB2par[1][3][4][0]=4.173118e-01; UB2par[1][3][4][1]=-1.560273e+00; UB2par[1][3][4][2]=2.642797e+00; UB2par[1][3][4][3]=-1.363096e+00; UB2par[1][3][4][4]=2.977245e-01; UB2par[1][3][4][5]=-2.378679e-02;

   UB2par[0][4][0][0]=1.285439e+00; UB2par[0][4][0][1]=-3.239525e+00; UB2par[0][4][0][2]=3.836452e+00; UB2par[0][4][0][3]=-1.740280e+00; UB2par[0][4][0][4]=3.480947e-01; UB2par[0][4][0][5]=-2.593508e-02;
   UB2par[0][4][1][0]=1.781593e+00; UB2par[0][4][1][1]=-4.345123e+00; UB2par[0][4][1][2]=4.675546e+00; UB2par[0][4][1][3]=-2.041562e+00; UB2par[0][4][1][4]=4.020156e-01; UB2par[0][4][1][5]=-2.978173e-02;
   UB2par[0][4][2][0]=2.207106e+00; UB2par[0][4][2][1]=-5.227136e+00; UB2par[0][4][2][2]=5.328840e+00; UB2par[0][4][2][3]=-2.282102e+00; UB2par[0][4][2][4]=4.480221e-01; UB2par[0][4][2][5]=-3.336511e-02;
   UB2par[0][4][3][0]=1.611553e+00; UB2par[0][4][3][1]=-3.556026e+00; UB2par[0][4][3][2]=3.642205e+00; UB2par[0][4][3][3]=-1.521852e+00; UB2par[0][4][3][4]=2.912436e-01; UB2par[0][4][3][5]=-2.122522e-02;
   UB2par[0][4][4][0]=1.338573e+00; UB2par[0][4][4][1]=-2.926073e+00; UB2par[0][4][4][2]=3.119397e+00; UB2par[0][4][4][3]=-1.341403e+00; UB2par[0][4][4][4]=2.657832e-01; UB2par[0][4][4][5]=-2.011407e-02;
   UB2par[1][4][0][0]=-2.497529e+00; UB2par[1][4][0][1]=4.777982e+00; UB2par[1][4][0][2]=-2.325765e+00; UB2par[1][4][0][3]=4.912960e-01; UB2par[1][4][0][4]=-4.141543e-02; UB2par[1][4][0][5]=5.748206e-04;
   UB2par[1][4][1][0]=-9.091682e-01; UB2par[1][4][1][1]=1.493638e+00; UB2par[1][4][1][2]=1.468828e-01; UB2par[1][4][1][3]=-3.966937e-01; UB2par[1][4][1][4]=1.147223e-01; UB2par[1][4][1][5]=-1.024044e-02;
   UB2par[1][4][2][0]=2.226544e-01; UB2par[1][4][2][1]=-8.403933e-01; UB2par[1][4][2][2]=1.926643e+00; UB2par[1][4][2][3]=-1.050867e+00; UB2par[1][4][2][4]=2.330042e-01; UB2par[1][4][2][5]=-1.865737e-02;
   UB2par[1][4][3][0]=-9.225142e-01; UB2par[1][4][3][1]=1.549144e+00; UB2par[1][4][3][2]=-2.549067e-02; UB2par[1][4][3][3]=-2.919517e-01; UB2par[1][4][3][4]=9.338154e-02; UB2par[1][4][3][5]=-8.845150e-03;
   UB2par[1][4][4][0]=-4.368881e-01; UB2par[1][4][4][1]=4.730827e-01; UB2par[1][4][4][2]=8.531009e-01; UB2par[1][4][4][3]=-6.416044e-01; UB2par[1][4][4][4]=1.617384e-01; UB2par[1][4][4][5]=-1.403291e-02;

   UB2par[0][5][0][0]=2.084655e+00; UB2par[0][5][0][1]=-4.748644e+00; UB2par[0][5][0][2]=4.873376e+00; UB2par[0][5][0][3]=-2.077100e+00; UB2par[0][5][0][4]=4.006154e-01; UB2par[0][5][0][5]=-2.908899e-02;
   UB2par[0][5][1][0]=1.829712e+00; UB2par[0][5][1][1]=-4.093840e+00; UB2par[0][5][1][2]=4.226933e+00; UB2par[0][5][1][3]=-1.796473e+00; UB2par[0][5][1][4]=3.467506e-01; UB2par[0][5][1][5]=-2.529631e-02;
   UB2par[0][5][2][0]=1.669298e+00; UB2par[0][5][2][1]=-3.489547e+00; UB2par[0][5][2][2]=3.513147e+00; UB2par[0][5][2][3]=-1.450781e+00; UB2par[0][5][2][4]=2.734387e-01; UB2par[0][5][2][5]=-1.956534e-02;
   UB2par[0][5][3][0]=1.073240e+00; UB2par[0][5][3][1]=-2.077969e+00; UB2par[0][5][3][2]=2.306393e+00; UB2par[0][5][3][3]=-9.902294e-01; UB2par[0][5][3][4]=1.939968e-01; UB2par[0][5][3][5]=-1.452679e-02;
   UB2par[0][5][4][0]=7.517885e-01; UB2par[0][5][4][1]=-1.153011e+00; UB2par[0][5][4][2]=1.430667e+00; UB2par[0][5][4][3]=-6.375793e-01; UB2par[0][5][4][4]=1.304485e-01; UB2par[0][5][4][5]=-1.025061e-02;
   UB2par[1][5][0][0]=-1.951543e+00; UB2par[1][5][0][1]=3.746882e+00; UB2par[1][5][0][2]=-1.552979e+00; UB2par[1][5][0][3]=2.104791e-01; UB2par[1][5][0][4]=7.818597e-03; UB2par[1][5][0][5]=-2.749669e-03;
   UB2par[1][5][1][0]=-3.064224e+00; UB2par[1][5][1][1]=6.149143e+00; UB2par[1][5][1][2]=-3.622579e+00; UB2par[1][5][1][3]=1.055975e+00; UB2par[1][5][1][4]=-1.544207e-01; UB2par[1][5][1][5]=9.074540e-03;
   UB2par[1][5][2][0]=-1.975366e+00; UB2par[1][5][2][1]=3.768114e+00; UB2par[1][5][2][2]=-1.708877e+00; UB2par[1][5][2][3]=3.181369e-01; UB2par[1][5][2][4]=-1.546524e-02; UB2par[1][5][2][5]=-1.151050e-03;
   UB2par[1][5][3][0]=-1.930385e+00; UB2par[1][5][3][1]=3.717191e+00; UB2par[1][5][3][2]=-1.640775e+00; UB2par[1][5][3][3]=2.582524e-01; UB2par[1][5][3][4]=5.391201e-03; UB2par[1][5][3][5]=-3.488753e-03;
   UB2par[1][5][4][0]=-1.565463e+00; UB2par[1][5][4][1]=2.861011e+00; UB2par[1][5][4][2]=-9.071125e-01; UB2par[1][5][4][3]=-3.924718e-02; UB2par[1][5][4][4]=6.281185e-02; UB2par[1][5][4][5]=-7.685000e-03;

   UB2par[0][6][0][0]=2.711059e+00; UB2par[0][6][0][1]=-5.995806e+00; UB2par[0][6][0][2]=5.805260e+00; UB2par[0][6][0][3]=-2.417107e+00; UB2par[0][6][0][4]=4.618836e-01; UB2par[0][6][0][5]=-3.345081e-02;
   UB2par[0][6][1][0]=1.847619e+00; UB2par[0][6][1][1]=-3.826811e+00; UB2par[0][6][1][2]=3.776117e+00; UB2par[0][6][1][3]=-1.550254e+00; UB2par[0][6][1][4]=2.910212e-01; UB2par[0][6][1][5]=-2.076674e-02;
   UB2par[0][6][2][0]=1.423600e+00; UB2par[0][6][2][1]=-2.791669e+00; UB2par[0][6][2][2]=2.831660e+00; UB2par[0][6][2][3]=-1.157030e+00; UB2par[0][6][2][4]=2.153278e-01; UB2par[0][6][2][5]=-1.524053e-02;
   UB2par[0][6][3][0]=7.549201e-01; UB2par[0][6][3][1]=-1.223091e+00; UB2par[0][6][3][2]=1.474812e+00; UB2par[0][6][3][3]=-6.243736e-01; UB2par[0][6][3][4]=1.189283e-01; UB2par[0][6][3][5]=-8.657824e-03;
   UB2par[0][6][4][0]=1.527656e+00; UB2par[0][6][4][1]=-2.772954e+00; UB2par[0][6][4][2]=2.792692e+00; UB2par[0][6][4][3]=-1.208088e+00; UB2par[0][6][4][4]=2.467089e-01; UB2par[0][6][4][5]=-1.938067e-02;
   UB2par[1][6][0][0]=2.015995e+00; UB2par[1][6][0][1]=-4.787487e+00; UB2par[1][6][0][2]=5.322447e+00; UB2par[1][6][0][3]=-2.431759e+00; UB2par[1][6][0][4]=4.970551e-01; UB2par[1][6][0][5]=-3.787044e-02;
   UB2par[1][6][1][0]=1.468450e+00; UB2par[1][6][1][1]=-2.987529e+00; UB2par[1][6][1][2]=3.360204e+00; UB2par[1][6][1][3]=-1.513169e+00; UB2par[1][6][1][4]=3.045097e-01; UB2par[1][6][1][5]=-2.292482e-02;
   UB2par[1][6][2][0]=7.351157e-01; UB2par[1][6][2][1]=-1.544663e+00; UB2par[1][6][2][2]=2.255947e+00; UB2par[1][6][2][3]=-1.109461e+00; UB2par[1][6][2][4]=2.342937e-01; UB2par[1][6][2][5]=-1.820923e-02;
   UB2par[1][6][3][0]=1.356232e+00; UB2par[1][6][3][1]=-2.904044e+00; UB2par[1][6][3][2]=3.390844e+00; UB2par[1][6][3][3]=-1.566863e+00; UB2par[1][6][3][4]=3.235801e-01; UB2par[1][6][3][5]=-2.492171e-02;
   UB2par[1][6][4][0]=2.968540e+00; UB2par[1][6][4][1]=-7.096999e+00; UB2par[1][6][4][2]=7.312686e+00; UB2par[1][6][4][3]=-3.263409e+00; UB2par[1][6][4][4]=6.684197e-01; UB2par[1][6][4][5]=-5.154209e-02;

   UB2par[0][7][0][0]=2.657876e+00; UB2par[0][7][0][1]=-5.410399e+00; UB2par[0][7][0][2]=5.070797e+00; UB2par[0][7][0][3]=-2.060421e+00; UB2par[0][7][0][4]=3.852099e-01; UB2par[0][7][0][5]=-2.733700e-02;
   UB2par[0][7][1][0]=2.044038e+00; UB2par[0][7][1][1]=-3.946625e+00; UB2par[0][7][1][2]=3.758150e+00; UB2par[0][7][1][3]=-1.527071e+00; UB2par[0][7][1][4]=2.865725e-01; UB2par[0][7][1][5]=-2.055097e-02;
   UB2par[0][7][2][0]=1.268726e+00; UB2par[0][7][2][1]=-1.960883e+00; UB2par[0][7][2][2]=1.911090e+00; UB2par[0][7][2][3]=-7.473166e-01; UB2par[0][7][2][4]=1.343111e-01; UB2par[0][7][2][5]=-9.303515e-03;
   UB2par[0][7][3][0]=1.207460e+00; UB2par[0][7][3][1]=-1.813255e+00; UB2par[0][7][3][2]=1.869222e+00; UB2par[0][7][3][3]=-7.934743e-01; UB2par[0][7][3][4]=1.580663e-01; UB2par[0][7][3][5]=-1.216305e-02;
   UB2par[0][7][4][0]=3.553467e-01; UB2par[0][7][4][1]=1.844911e-01; UB2par[0][7][4][2]=1.581532e-01; UB2par[0][7][4][3]=-1.232678e-01; UB2par[0][7][4][4]=3.548751e-02; UB2par[0][7][4][5]=-3.624764e-03;
   UB2par[1][7][0][0]=3.136522e+00; UB2par[1][7][0][1]=-7.243286e+00; UB2par[1][7][0][2]=7.352396e+00; UB2par[1][7][0][3]=-3.230815e+00; UB2par[1][7][0][4]=6.479641e-01; UB2par[1][7][0][5]=-4.887631e-02;
   UB2par[1][7][1][0]=2.418667e+00; UB2par[1][7][1][1]=-5.352323e+00; UB2par[1][7][1][2]=5.472072e+00; UB2par[1][7][1][3]=-2.388867e+00; UB2par[1][7][1][4]=4.760871e-01; UB2par[1][7][1][5]=-3.577361e-02;
   UB2par[1][7][2][0]=1.634905e+00; UB2par[1][7][2][1]=-3.656430e+00; UB2par[1][7][2][2]=4.077435e+00; UB2par[1][7][2][3]=-1.849733e+00; UB2par[1][7][2][4]=3.776857e-01; UB2par[1][7][2][5]=-2.887912e-02;
   UB2par[1][7][3][0]=2.024640e+00; UB2par[1][7][3][1]=-4.729774e+00; UB2par[1][7][3][2]=5.142966e+00; UB2par[1][7][3][3]=-2.339635e+00; UB2par[1][7][3][4]=4.834003e-01; UB2par[1][7][3][5]=-3.747466e-02;
   UB2par[1][7][4][0]=1.818338e+00; UB2par[1][7][4][1]=-4.737713e+00; UB2par[1][7][4][2]=5.445835e+00; UB2par[1][7][4][3]=-2.546778e+00; UB2par[1][7][4][4]=5.347708e-01; UB2par[1][7][4][5]=-4.184603e-02;

   UB2par[0][8][0][0]=2.905820e+00; UB2par[0][8][0][1]=-5.403612e+00; UB2par[0][8][0][2]=4.773328e+00; UB2par[0][8][0][3]=-1.864284e+00; UB2par[0][8][0][4]=3.372965e-01; UB2par[0][8][0][5]=-2.323027e-02;
   UB2par[0][8][1][0]=2.163936e+00; UB2par[0][8][1][1]=-3.568152e+00; UB2par[0][8][1][2]=3.104464e+00; UB2par[0][8][1][3]=-1.179917e+00; UB2par[0][8][1][4]=2.090872e-01; UB2par[0][8][1][5]=-1.425177e-02;
   UB2par[0][8][2][0]=1.417903e+00; UB2par[0][8][2][1]=-1.836760e+00; UB2par[0][8][2][2]=1.643967e+00; UB2par[0][8][2][3]=-6.186578e-01; UB2par[0][8][2][4]=1.093766e-01; UB2par[0][8][2][5]=-7.562962e-03;
   UB2par[0][8][3][0]=9.390700e-01; UB2par[0][8][3][1]=-7.392230e-01; UB2par[0][8][3][2]=7.668078e-01; UB2par[0][8][3][3]=-3.149342e-01; UB2par[0][8][3][4]=6.347251e-02; UB2par[0][8][3][5]=-5.137739e-03;
   UB2par[0][8][4][0]=3.972424e-01; UB2par[0][8][4][1]=6.422562e-01; UB2par[0][8][4][2]=-4.800705e-01; UB2par[0][8][4][3]=1.931318e-01; UB2par[0][8][4][4]=-3.311252e-02; UB2par[0][8][4][5]=1.897884e-03;
   UB2par[1][8][0][0]=3.344293e+00; UB2par[1][8][0][1]=-7.755463e+00; UB2par[1][8][0][2]=7.768755e+00; UB2par[1][8][0][3]=-3.382407e+00; UB2par[1][8][0][4]=6.736050e-01; UB2par[1][8][0][5]=-5.051444e-02;
   UB2par[1][8][1][0]=2.764631e+00; UB2par[1][8][1][1]=-6.289682e+00; UB2par[1][8][1][2]=6.348763e+00; UB2par[1][8][1][3]=-2.765330e+00; UB2par[1][8][1][4]=5.522346e-01; UB2par[1][8][1][5]=-4.165144e-02;
   UB2par[1][8][2][0]=2.153540e+00; UB2par[1][8][2][1]=-5.069325e+00; UB2par[1][8][2][2]=5.392973e+00; UB2par[1][8][2][3]=-2.407226e+00; UB2par[1][8][2][4]=4.882704e-01; UB2par[1][8][2][5]=-3.722349e-02;
   UB2par[1][8][3][0]=2.537559e+00; UB2par[1][8][3][1]=-6.225083e+00; UB2par[1][8][3][2]=6.533697e+00; UB2par[1][8][3][3]=-2.912997e+00; UB2par[1][8][3][4]=5.925981e-01; UB2par[1][8][3][5]=-4.533045e-02;
   UB2par[1][8][4][0]=3.032710e+00; UB2par[1][8][4][1]=-7.641679e+00; UB2par[1][8][4][2]=7.866758e+00; UB2par[1][8][4][3]=-3.474160e+00; UB2par[1][8][4][4]=7.020032e-01; UB2par[1][8][4][5]=-5.335804e-02;

   UB2par[0][9][0][0]=2.346968e+00; UB2par[0][9][0][1]=-3.678269e+00; UB2par[0][9][0][2]=3.071169e+00; UB2par[0][9][0][3]=-1.119256e+00; UB2par[0][9][0][4]=1.864866e-01; UB2par[0][9][0][5]=-1.169278e-02;
   UB2par[0][9][1][0]=2.127543e+00; UB2par[0][9][1][1]=-3.109052e+00; UB2par[0][9][1][2]=2.529113e+00; UB2par[0][9][1][3]=-9.007073e-01; UB2par[0][9][1][4]=1.487469e-01; UB2par[0][9][1][5]=-9.388855e-03;
   UB2par[0][9][2][0]=8.686508e-01; UB2par[0][9][2][1]=-3.471718e-01; UB2par[0][9][2][2]=2.951407e-01; UB2par[0][9][2][3]=-6.583384e-02; UB2par[0][9][2][4]=3.163743e-03; UB2par[0][9][2][5]=2.287783e-04;
   UB2par[0][9][3][0]=2.576244e-01; UB2par[0][9][3][1]=1.130220e+00; UB2par[0][9][3][2]=-9.547378e-01; UB2par[0][9][3][3]=4.053021e-01; UB2par[0][9][3][4]=-7.809268e-02; UB2par[0][9][3][5]=5.491179e-03;
   UB2par[0][9][4][0]=-4.486266e-01; UB2par[0][9][4][1]=2.779397e+00; UB2par[0][9][4][2]=-2.338857e+00; UB2par[0][9][4][3]=9.301942e-01; UB2par[0][9][4][4]=-1.700400e-01; UB2par[0][9][4][5]=1.156551e-02;
   UB2par[1][9][0][0]=3.494064e+00; UB2par[1][9][0][1]=-8.145780e+00; UB2par[1][9][0][2]=8.084078e+00; UB2par[1][9][0][3]=-3.492261e+00; UB2par[1][9][0][4]=6.904550e-01; UB2par[1][9][0][5]=-5.140959e-02;
   UB2par[1][9][1][0]=3.089460e+00; UB2par[1][9][1][1]=-7.176241e+00; UB2par[1][9][1][2]=7.140645e+00; UB2par[1][9][1][3]=-3.084590e+00; UB2par[1][9][1][4]=6.121459e-01; UB2par[1][9][1][5]=-4.589998e-02;
   UB2par[1][9][2][0]=3.177584e+00; UB2par[1][9][2][1]=-7.495381e+00; UB2par[1][9][2][2]=7.472120e+00; UB2par[1][9][2][3]=-3.240005e+00; UB2par[1][9][2][4]=6.469063e-01; UB2par[1][9][2][5]=-4.884229e-02;
   UB2par[1][9][3][0]=2.979414e+00; UB2par[1][9][3][1]=-7.279987e+00; UB2par[1][9][3][2]=7.425361e+00; UB2par[1][9][3][3]=-3.265295e+00; UB2par[1][9][3][4]=6.592579e-01; UB2par[1][9][3][5]=-5.022707e-02;
   UB2par[1][9][4][0]=3.271879e+00; UB2par[1][9][4][1]=-8.184427e+00; UB2par[1][9][4][2]=8.311916e+00; UB2par[1][9][4][3]=-3.644978e+00; UB2par[1][9][4][4]=7.335770e-01; UB2par[1][9][4][5]=-5.565883e-02;

   UB2par[0][10][0][0]=2.177469e+00; UB2par[0][10][0][1]=-2.961661e+00; UB2par[0][10][0][2]=2.331219e+00; UB2par[0][10][0][3]=-7.975313e-01; UB2par[0][10][0][4]=1.229659e-01; UB2par[0][10][0][5]=-6.991903e-03;
   UB2par[0][10][1][0]=6.574309e-01; UB2par[0][10][1][1]=2.107697e-01; UB2par[0][10][1][2]=-2.463292e-01; UB2par[0][10][1][3]=1.943871e-01; UB2par[0][10][1][4]=-5.752245e-02; UB2par[0][10][1][5]=5.559949e-03;
   UB2par[0][10][2][0]=3.836240e-01; UB2par[0][10][2][1]=1.185242e+00; UB2par[0][10][2][2]=-1.221666e+00; UB2par[0][10][2][3]=6.007262e-01; UB2par[0][10][2][4]=-1.328669e-01; UB2par[0][10][2][5]=1.076134e-02;
   UB2par[0][10][3][0]=1.791614e-01; UB2par[0][10][3][1]=1.500440e+00; UB2par[0][10][3][2]=-1.297412e+00; UB2par[0][10][3][3]=5.306637e-01; UB2par[0][10][3][4]=-9.768114e-02; UB2par[0][10][3][5]=6.530267e-03;
   UB2par[0][10][4][0]=-7.909648e-01; UB2par[0][10][4][1]=3.583349e+00; UB2par[0][10][4][2]=-2.928247e+00; UB2par[0][10][4][3]=1.116760e+00; UB2par[0][10][4][4]=-1.967860e-01; UB2par[0][10][4][5]=1.299287e-02;
   UB2par[1][10][0][0]=3.359695e+00; UB2par[1][10][0][1]=-7.765826e+00; UB2par[1][10][0][2]=7.704923e+00; UB2par[1][10][0][3]=-3.322377e+00; UB2par[1][10][0][4]=6.557955e-01; UB2par[1][10][0][5]=-4.878552e-02;
   UB2par[1][10][1][0]=3.145652e+00; UB2par[1][10][1][1]=-7.251625e+00; UB2par[1][10][1][2]=7.151692e+00; UB2par[1][10][1][3]=-3.067256e+00; UB2par[1][10][1][4]=6.049388e-01; UB2par[1][10][1][5]=-4.511736e-02;
   UB2par[1][10][2][0]=2.853610e+00; UB2par[1][10][2][1]=-6.542628e+00; UB2par[1][10][2][2]=6.485972e+00; UB2par[1][10][2][3]=-2.777831e+00; UB2par[1][10][2][4]=5.465230e-01; UB2par[1][10][2][5]=-4.062546e-02;
   UB2par[1][10][3][0]=3.573842e+00; UB2par[1][10][3][1]=-8.622523e+00; UB2par[1][10][3][2]=8.546516e+00; UB2par[1][10][3][3]=-3.709095e+00; UB2par[1][10][3][4]=7.435654e-01; UB2par[1][10][3][5]=-5.641833e-02;
   UB2par[1][10][4][0]=3.733883e+00; UB2par[1][10][4][1]=-9.108403e+00; UB2par[1][10][4][2]=9.028468e+00; UB2par[1][10][4][3]=-3.918286e+00; UB2par[1][10][4][4]=7.847081e-01; UB2par[1][10][4][5]=-5.939244e-02;

   UB2par[0][11][0][0]=1.658758e+00; UB2par[0][11][0][1]=-1.526272e+00; UB2par[0][11][0][2]=1.030357e+00; UB2par[0][11][0][3]=-2.660524e-01; UB2par[0][11][0][4]=2.117715e-02; UB2par[0][11][0][5]=4.641383e-04;
   UB2par[0][11][1][0]=9.335576e-01; UB2par[0][11][1][1]=8.652998e-02; UB2par[0][11][1][2]=-3.560279e-01; UB2par[0][11][1][3]=2.826481e-01; UB2par[0][11][1][4]=-7.895966e-02; UB2par[0][11][1][5]=7.343175e-03;
   UB2par[0][11][2][0]=6.363671e-01; UB2par[0][11][2][1]=8.138053e-01; UB2par[0][11][2][2]=-9.484989e-01; UB2par[0][11][2][3]=4.815849e-01; UB2par[0][11][2][4]=-1.058256e-01; UB2par[0][11][2][5]=8.382566e-03;
   UB2par[0][11][3][0]=1.438236e-01; UB2par[0][11][3][1]=1.903135e+00; UB2par[0][11][3][2]=-1.804015e+00; UB2par[0][11][3][3]=7.789633e-01; UB2par[0][11][3][4]=-1.520503e-01; UB2par[0][11][3][5]=1.099002e-02;
   UB2par[0][11][4][0]=-4.040885e-01; UB2par[0][11][4][1]=3.022838e+00; UB2par[0][11][4][2]=-2.576063e+00; UB2par[0][11][4][3]=9.991099e-01; UB2par[0][11][4][4]=-1.767874e-01; UB2par[0][11][4][5]=1.164926e-02;
   UB2par[1][11][0][0]=3.809638e+00; UB2par[1][11][0][1]=-8.706262e+00; UB2par[1][11][0][2]=8.464219e+00; UB2par[1][11][0][3]=-3.619583e+00; UB2par[1][11][0][4]=7.121904e-01; UB2par[1][11][0][5]=-5.294458e-02;
   UB2par[1][11][1][0]=3.282346e+00; UB2par[1][11][1][1]=-7.434833e+00; UB2par[1][11][1][2]=7.218299e+00; UB2par[1][11][1][3]=-3.063056e+00; UB2par[1][11][1][4]=5.987158e-01; UB2par[1][11][1][5]=-4.427925e-02;
   UB2par[1][11][2][0]=3.456720e+00; UB2par[1][11][2][1]=-7.904231e+00; UB2par[1][11][2][2]=7.631400e+00; UB2par[1][11][2][3]=-3.236607e+00; UB2par[1][11][2][4]=6.348697e-01; UB2par[1][11][2][5]=-4.720990e-02;
   UB2par[1][11][3][0]=3.540873e+00; UB2par[1][11][3][1]=-8.277984e+00; UB2par[1][11][3][2]=8.031433e+00; UB2par[1][11][3][3]=-3.412270e+00; UB2par[1][11][3][4]=6.695029e-01; UB2par[1][11][3][5]=-4.972202e-02;
   UB2par[1][11][4][0]=4.246127e+00; UB2par[1][11][4][1]=-1.028567e+01; UB2par[1][11][4][2]=1.000095e+01; UB2par[1][11][4][3]=-4.286273e+00; UB2par[1][11][4][4]=8.495402e-01; UB2par[1][11][4][5]=-6.370729e-02;

   UB2par[0][12][0][0]=1.562439e+00; UB2par[0][12][0][1]=-1.246427e+00; UB2par[0][12][0][2]=8.059555e-01; UB2par[0][12][0][3]=-1.874989e-01; UB2par[0][12][0][4]=8.375868e-03; UB2par[0][12][0][5]=1.273863e-03;
   UB2par[0][12][1][0]=2.129176e-01; UB2par[0][12][1][1]=1.791816e+00; UB2par[0][12][1][2]=-1.788674e+00; UB2par[0][12][1][3]=8.398845e-01; UB2par[0][12][1][4]=-1.815186e-01; UB2par[0][12][1][5]=1.458473e-02;
   UB2par[0][12][2][0]=-6.737713e-01; UB2par[0][12][2][1]=3.835770e+00; UB2par[0][12][2][2]=-3.469195e+00; UB2par[0][12][2][3]=1.464146e+00; UB2par[0][12][2][4]=-2.878907e-01; UB2par[0][12][2][5]=2.134617e-02;
   UB2par[0][12][3][0]=-5.181144e-01; UB2par[0][12][3][1]=3.403146e+00; UB2par[0][12][3][2]=-2.989920e+00; UB2par[0][12][3][3]=1.207725e+00; UB2par[0][12][3][4]=-2.248841e-01; UB2par[0][12][3][5]=1.571244e-02;
   UB2par[0][12][4][0]=-9.572171e-01; UB2par[0][12][4][1]=4.246036e+00; UB2par[0][12][4][2]=-3.521491e+00; UB2par[0][12][4][3]=1.329973e+00; UB2par[0][12][4][4]=-2.303512e-01; UB2par[0][12][4][5]=1.488615e-02;
   UB2par[1][12][0][0]=4.477527e+00; UB2par[1][12][0][1]=-1.031302e+01; UB2par[1][12][0][2]=9.869055e+00; UB2par[1][12][0][3]=-4.189304e+00; UB2par[1][12][0][4]=8.213700e-01; UB2par[1][12][0][5]=-6.094276e-02;
   UB2par[1][12][1][0]=3.561608e+00; UB2par[1][12][1][1]=-7.928859e+00; UB2par[1][12][1][2]=7.530223e+00; UB2par[1][12][1][3]=-3.153128e+00; UB2par[1][12][1][4]=6.105884e-01; UB2par[1][12][1][5]=-4.482956e-02;
   UB2par[1][12][2][0]=3.243664e+00; UB2par[1][12][2][1]=-7.314095e+00; UB2par[1][12][2][2]=7.053595e+00; UB2par[1][12][2][3]=-2.976904e+00; UB2par[1][12][2][4]=5.803590e-01; UB2par[1][12][2][5]=-4.288280e-02;
   UB2par[1][12][3][0]=3.506284e+00; UB2par[1][12][3][1]=-8.201218e+00; UB2par[1][12][3][2]=7.992399e+00; UB2par[1][12][3][3]=-3.413543e+00; UB2par[1][12][3][4]=6.734527e-01; UB2par[1][12][3][5]=-5.027994e-02;
   UB2par[1][12][4][0]=4.570066e+00; UB2par[1][12][4][1]=-1.087951e+01; UB2par[1][12][4][2]=1.039684e+01; UB2par[1][12][4][3]=-4.407434e+00; UB2par[1][12][4][4]=8.666905e-01; UB2par[1][12][4][5]=-6.461171e-02;

   UB2par[0][13][0][0]=4.651761e-01; UB2par[0][13][0][1]=1.386989e+00; UB2par[0][13][0][2]=-1.482660e+00; UB2par[0][13][0][3]=7.430544e-01; UB2par[0][13][0][4]=-1.716450e-01; UB2par[0][13][0][5]=1.466028e-02;
   UB2par[0][13][1][0]=4.184686e-02; UB2par[0][13][1][1]=2.214082e+00; UB2par[0][13][1][2]=-2.085310e+00; UB2par[0][13][1][3]=9.222258e-01; UB2par[0][13][1][4]=-1.898444e-01; UB2par[0][13][1][5]=1.467492e-02;
   UB2par[0][13][2][0]=-8.343744e-01; UB2par[0][13][2][1]=3.931340e+00; UB2par[0][13][2][2]=-3.334097e+00; UB2par[0][13][2][3]=1.331538e+00; UB2par[0][13][2][4]=-2.492490e-01; UB2par[0][13][2][5]=1.764707e-02;
   UB2par[0][13][3][0]=-6.083975e-01; UB2par[0][13][3][1]=3.385385e+00; UB2par[0][13][3][2]=-2.811087e+00; UB2par[0][13][3][3]=1.086536e+00; UB2par[0][13][3][4]=-1.952061e-01; UB2par[0][13][3][5]=1.322624e-02;
   UB2par[0][13][4][0]=-1.253713e+00; UB2par[0][13][4][1]=4.810240e+00; UB2par[0][13][4][2]=-3.887443e+00; UB2par[0][13][4][3]=1.444026e+00; UB2par[0][13][4][4]=-2.489760e-01; UB2par[0][13][4][5]=1.620136e-02;
   UB2par[1][13][0][0]=4.502975e+00; UB2par[1][13][0][1]=-1.006127e+01; UB2par[1][13][0][2]=9.430209e+00; UB2par[1][13][0][3]=-3.932552e+00; UB2par[1][13][0][4]=7.581105e-01; UB2par[1][13][0][5]=-5.532625e-02;
   UB2par[1][13][1][0]=3.848215e+00; UB2par[1][13][1][1]=-8.514149e+00; UB2par[1][13][1][2]=7.989920e+00; UB2par[1][13][1][3]=-3.324391e+00; UB2par[1][13][1][4]=6.410473e-01; UB2par[1][13][1][5]=-4.691264e-02;
   UB2par[1][13][2][0]=3.499624e+00; UB2par[1][13][2][1]=-7.951410e+00; UB2par[1][13][2][2]=7.657704e+00; UB2par[1][13][2][3]=-3.246164e+00; UB2par[1][13][2][4]=6.371928e-01; UB2par[1][13][2][5]=-4.746401e-02;
   UB2par[1][13][3][0]=3.712818e+00; UB2par[1][13][3][1]=-8.537213e+00; UB2par[1][13][3][2]=8.154837e+00; UB2par[1][13][3][3]=-3.427599e+00; UB2par[1][13][3][4]=6.666713e-01; UB2par[1][13][3][5]=-4.914590e-02;
   UB2par[1][13][4][0]=4.231214e+00; UB2par[1][13][4][1]=-1.000954e+01; UB2par[1][13][4][2]=9.562125e+00; UB2par[1][13][4][3]=-4.024964e+00; UB2par[1][13][4][4]=7.831821e-01; UB2par[1][13][4][5]=-5.767392e-02;

   UB2par[0][14][0][0]=-2.463798e-02; UB2par[0][14][0][1]=2.559053e+00; UB2par[0][14][0][2]=-2.425486e+00; UB2par[0][14][0][3]=1.086200e+00; UB2par[0][14][0][4]=-2.297779e-01; UB2par[0][14][0][5]=1.837265e-02;
   UB2par[0][14][1][0]=-1.212086e-01; UB2par[0][14][1][1]=2.936540e+00; UB2par[0][14][1][2]=-2.881944e+00; UB2par[0][14][1][3]=1.293958e+00; UB2par[0][14][1][4]=-2.690273e-01; UB2par[0][14][1][5]=2.102102e-02;
   UB2par[0][14][2][0]=-5.629532e-01; UB2par[0][14][2][1]=3.455129e+00; UB2par[0][14][2][2]=-2.958781e+00; UB2par[0][14][2][3]=1.181274e+00; UB2par[0][14][2][4]=-2.203677e-01; UB2par[0][14][2][5]=1.554895e-02;
   UB2par[0][14][3][0]=-1.436075e+00; UB2par[0][14][3][1]=5.416960e+00; UB2par[0][14][3][2]=-4.559062e+00; UB2par[0][14][3][3]=1.778548e+00; UB2par[0][14][3][4]=-3.244556e-01; UB2par[0][14][3][5]=2.246640e-02;
   UB2par[0][14][4][0]=-2.479276e+00; UB2par[0][14][4][1]=7.547962e+00; UB2par[0][14][4][2]=-6.148352e+00; UB2par[0][14][4][3]=2.324943e+00; UB2par[0][14][4][4]=-4.131918e-01; UB2par[0][14][4][5]=2.802386e-02;
   UB2par[1][14][0][0]=4.116035e+00; UB2par[1][14][0][1]=-9.116068e+00; UB2par[1][14][0][2]=8.600134e+00; UB2par[1][14][0][3]=-3.595673e+00; UB2par[1][14][0][4]=6.939844e-01; UB2par[1][14][0][5]=-5.069775e-02;
   UB2par[1][14][1][0]=3.732373e+00; UB2par[1][14][1][1]=-8.204526e+00; UB2par[1][14][1][2]=7.700737e+00; UB2par[1][14][1][3]=-3.197893e+00; UB2par[1][14][1][4]=6.146011e-01; UB2par[1][14][1][5]=-4.478623e-02;
   UB2par[1][14][2][0]=3.764316e+00; UB2par[1][14][2][1]=-8.426572e+00; UB2par[1][14][2][2]=7.963087e+00; UB2par[1][14][2][3]=-3.330254e+00; UB2par[1][14][2][4]=6.459151e-01; UB2par[1][14][2][5]=-4.756131e-02;
   UB2par[1][14][3][0]=3.982794e+00; UB2par[1][14][3][1]=-8.869277e+00; UB2par[1][14][3][2]=8.273149e+00; UB2par[1][14][3][3]=-3.431110e+00; UB2par[1][14][3][4]=6.617425e-01; UB2par[1][14][3][5]=-4.851313e-02;
   UB2par[1][14][4][0]=3.848733e+00; UB2par[1][14][4][1]=-8.988308e+00; UB2par[1][14][4][2]=8.606553e+00; UB2par[1][14][4][3]=-3.613069e+00; UB2par[1][14][4][4]=6.995028e-01; UB2par[1][14][4][5]=-5.118224e-02;

   UB2par[0][15][0][0]=-3.433772e-01; UB2par[0][15][0][1]=3.284685e+00; UB2par[0][15][0][2]=-2.997582e+00; UB2par[0][15][0][3]=1.290633e+00; UB2par[0][15][0][4]=-2.636395e-01; UB2par[0][15][0][5]=2.047829e-02;
   UB2par[0][15][1][0]=-1.338058e+00; UB2par[0][15][1][1]=5.360676e+00; UB2par[0][15][1][2]=-4.664573e+00; UB2par[0][15][1][3]=1.908323e+00; UB2par[0][15][1][4]=-3.691079e-01; UB2par[0][15][1][5]=2.721153e-02;
   UB2par[0][15][2][0]=-2.275317e+00; UB2par[0][15][2][1]=7.329818e+00; UB2par[0][15][2][2]=-6.212294e+00; UB2par[0][15][2][3]=2.471882e+00; UB2par[0][15][2][4]=-4.650116e-01; UB2par[0][15][2][5]=3.341979e-02;
   UB2par[0][15][3][0]=-1.690499e+00; UB2par[0][15][3][1]=5.803248e+00; UB2par[0][15][3][2]=-4.753060e+00; UB2par[0][15][3][3]=1.814582e+00; UB2par[0][15][3][4]=-3.252898e-01; UB2par[0][15][3][5]=2.218721e-02;
   UB2par[0][15][4][0]=-2.256235e+00; UB2par[0][15][4][1]=6.869586e+00; UB2par[0][15][4][2]=-5.453160e+00; UB2par[0][15][4][3]=2.007785e+00; UB2par[0][15][4][4]=-3.470750e-01; UB2par[0][15][4][5]=2.288852e-02;
   UB2par[1][15][0][0]=4.208171e+00; UB2par[1][15][0][1]=-9.240764e+00; UB2par[1][15][0][2]=8.643142e+00; UB2par[1][15][0][3]=-3.586463e+00; UB2par[1][15][0][4]=6.867690e-01; UB2par[1][15][0][5]=-4.974183e-02;
   UB2par[1][15][1][0]=3.440967e+00; UB2par[1][15][1][1]=-7.410688e+00; UB2par[1][15][1][2]=6.949274e+00; UB2par[1][15][1][3]=-2.874155e+00; UB2par[1][15][1][4]=5.498822e-01; UB2par[1][15][1][5]=-3.991614e-02;
   UB2par[1][15][2][0]=3.206763e+00; UB2par[1][15][2][1]=-7.005830e+00; UB2par[1][15][2][2]=6.676694e+00; UB2par[1][15][2][3]=-2.791627e+00; UB2par[1][15][2][4]=5.396559e-01; UB2par[1][15][2][5]=-3.956680e-02;
   UB2par[1][15][3][0]=3.485531e+00; UB2par[1][15][3][1]=-7.829809e+00; UB2par[1][15][3][2]=7.450966e+00; UB2par[1][15][3][3]=-3.109693e+00; UB2par[1][15][3][4]=5.995560e-01; UB2par[1][15][3][5]=-4.378617e-02;
   UB2par[1][15][4][0]=4.474056e+00; UB2par[1][15][4][1]=-1.020137e+01; UB2par[1][15][4][2]=9.501161e+00; UB2par[1][15][4][3]=-3.924738e+00; UB2par[1][15][4][4]=7.511434e-01; UB2par[1][15][4][5]=-5.446392e-02;

   UB2par[0][16][0][0]=-7.271022e-01; UB2par[0][16][0][1]=4.251599e+00; UB2par[0][16][0][2]=-3.812971e+00; UB2par[0][16][0][3]=1.603126e+00; UB2par[0][16][0][4]=-3.197050e-01; UB2par[0][16][0][5]=2.430450e-02;
   UB2par[0][16][1][0]=-2.190645e-02; UB2par[0][16][1][1]=2.482472e+00; UB2par[0][16][1][2]=-2.274316e+00; UB2par[0][16][1][3]=9.663633e-01; UB2par[0][16][1][4]=-1.920389e-01; UB2par[0][16][1][5]=1.442458e-02;
   UB2par[0][16][2][0]=-1.537281e+00; UB2par[0][16][2][1]=5.712911e+00; UB2par[0][16][2][2]=-4.849739e+00; UB2par[0][16][2][3]=1.921487e+00; UB2par[0][16][2][4]=-3.579699e-01; UB2par[0][16][2][5]=2.536022e-02;
   UB2par[0][16][3][0]=-1.219828e+00; UB2par[0][16][3][1]=4.874358e+00; UB2par[0][16][3][2]=-4.005396e+00; UB2par[0][16][3][3]=1.519307e+00; UB2par[0][16][3][4]=-2.691336e-01; UB2par[0][16][3][5]=1.808227e-02;
   UB2par[0][16][4][0]=-1.801104e+00; UB2par[0][16][4][1]=5.925238e+00; UB2par[0][16][4][2]=-4.649301e+00; UB2par[0][16][4][3]=1.669341e+00; UB2par[0][16][4][4]=-2.777820e-01; UB2par[0][16][4][5]=1.739444e-02;
   UB2par[1][16][0][0]=4.260523e+00; UB2par[1][16][0][1]=-9.283696e+00; UB2par[1][16][0][2]=8.655766e+00; UB2par[1][16][0][3]=-3.587737e+00; UB2par[1][16][0][4]=6.872517e-01; UB2par[1][16][0][5]=-4.986409e-02;
   UB2par[1][16][1][0]=3.805104e+00; UB2par[1][16][1][1]=-8.302312e+00; UB2par[1][16][1][2]=7.741531e+00; UB2par[1][16][1][3]=-3.194038e+00; UB2par[1][16][1][4]=6.096257e-01; UB2par[1][16][1][5]=-4.411159e-02;
   UB2par[1][16][2][0]=3.353862e+00; UB2par[1][16][2][1]=-7.306169e+00; UB2par[1][16][2][2]=6.903903e+00; UB2par[1][16][2][3]=-2.869741e+00; UB2par[1][16][2][4]=5.520956e-01; UB2par[1][16][2][5]=-4.031738e-02;
   UB2par[1][16][3][0]=3.474128e+00; UB2par[1][16][3][1]=-7.527206e+00; UB2par[1][16][3][2]=7.017126e+00; UB2par[1][16][3][3]=-2.879108e+00; UB2par[1][16][3][4]=5.465570e-01; UB2par[1][16][3][5]=-3.934291e-02;
   UB2par[1][16][4][0]=2.893461e+00; UB2par[1][16][4][1]=-6.649288e+00; UB2par[1][16][4][2]=6.513033e+00; UB2par[1][16][4][3]=-2.735246e+00; UB2par[1][16][4][4]=5.252217e-01; UB2par[1][16][4][5]=-3.797804e-02;

   UB2par[0][17][0][0]=-1.698647e+00; UB2par[0][17][0][1]=5.900868e+00; UB2par[0][17][0][2]=-4.895766e+00; UB2par[0][17][0][3]=1.952275e+00; UB2par[0][17][0][4]=-3.758750e-01; UB2par[0][17][0][5]=2.795584e-02;
   UB2par[0][17][1][0]=-1.254416e+00; UB2par[0][17][1][1]=4.906399e+00; UB2par[0][17][1][2]=-4.078116e+00; UB2par[0][17][1][3]=1.609768e+00; UB2par[0][17][1][4]=-3.029187e-01; UB2par[0][17][1][5]=2.185798e-02;
   UB2par[0][17][2][0]=-1.725802e+00; UB2par[0][17][2][1]=5.562516e+00; UB2par[0][17][2][2]=-4.361190e+00; UB2par[0][17][2][3]=1.623225e+00; UB2par[0][17][2][4]=-2.878108e-01; UB2par[0][17][2][5]=1.960037e-02;
   UB2par[0][17][3][0]=-1.126829e+00; UB2par[0][17][3][1]=3.934458e+00; UB2par[0][17][3][2]=-2.761437e+00; UB2par[0][17][3][3]=8.812602e-01; UB2par[0][17][3][4]=-1.256688e-01; UB2par[0][17][3][5]=6.234254e-03;
   UB2par[0][17][4][0]=-2.512281e-01; UB2par[0][17][4][1]=2.118118e+00; UB2par[0][17][4][2]=-1.265232e+00; UB2par[0][17][4][3]=2.748453e-01; UB2par[0][17][4][4]=-7.572927e-03; UB2par[0][17][4][5]=-2.526716e-03;
   UB2par[1][17][0][0]=4.155101e+00; UB2par[1][17][0][1]=-9.041050e+00; UB2par[1][17][0][2]=8.416989e+00; UB2par[1][17][0][3]=-3.470763e+00; UB2par[1][17][0][4]=6.598355e-01; UB2par[1][17][0][5]=-4.743072e-02;
   UB2par[1][17][1][0]=3.449245e+00; UB2par[1][17][1][1]=-7.378494e+00; UB2par[1][17][1][2]=6.888591e+00; UB2par[1][17][1][3]=-2.832395e+00; UB2par[1][17][1][4]=5.379507e-01; UB2par[1][17][1][5]=-3.872496e-02;
   UB2par[1][17][2][0]=3.398145e+00; UB2par[1][17][2][1]=-7.406104e+00; UB2par[1][17][2][2]=6.996625e+00; UB2par[1][17][2][3]=-2.903070e+00; UB2par[1][17][2][4]=5.562132e-01; UB2par[1][17][2][5]=-4.036110e-02;
   UB2par[1][17][3][0]=3.785680e+00; UB2par[1][17][3][1]=-8.205369e+00; UB2par[1][17][3][2]=7.560917e+00; UB2par[1][17][3][3]=-3.084746e+00; UB2par[1][17][3][4]=5.838110e-01; UB2par[1][17][3][5]=-4.195845e-02;
   UB2par[1][17][4][0]=3.503734e+00; UB2par[1][17][4][1]=-7.561969e+00; UB2par[1][17][4][2]=6.962857e+00; UB2par[1][17][4][3]=-2.802755e+00; UB2par[1][17][4][4]=5.191409e-01; UB2par[1][17][4][5]=-3.630783e-02;

   // For the last 3 bins in costheta:
   // The 2 last bins have the same parameters.
   // The parameters are the same for layers [10,11], [12,13] and [14,15,16,17].
   UB2par[0][2][5][0]=2.432210e+00; UB2par[0][2][5][1]=-5.443225e+00; UB2par[0][2][5][2]=5.201100e+00; UB2par[0][2][5][3]=-2.127771e+00; UB2par[0][2][5][4]=4.056839e-01; UB2par[0][2][5][5]=-2.964222e-02;
   UB2par[0][2][6][0]=1.634909e+00; UB2par[0][2][6][1]=-3.325845e+00; UB2par[0][2][6][2]=3.141647e+00; UB2par[0][2][6][3]=-1.238921e+00; UB2par[0][2][6][4]=2.304975e-01; UB2par[0][2][6][5]=-1.673315e-02;
   UB2par[0][2][7][0]=1.634910e+00; UB2par[0][2][7][1]=-3.325845e+00; UB2par[0][2][7][2]=3.141648e+00; UB2par[0][2][7][3]=-1.238921e+00; UB2par[0][2][7][4]=2.304975e-01; UB2par[0][2][7][5]=-1.673315e-02;
   UB2par[1][2][5][0]=-6.691986e-02; UB2par[1][2][5][1]=-3.694283e-02; UB2par[1][2][5][2]=9.762372e-01; UB2par[1][2][5][3]=-5.758933e-01; UB2par[1][2][5][4]=1.312484e-01; UB2par[1][2][5][5]=-1.075686e-02;
   UB2par[1][2][6][0]=-7.987131e+00; UB2par[1][2][6][1]=1.623524e+01; UB2par[1][2][6][2]=-1.192663e+01; UB2par[1][2][6][3]=4.336764e+00; UB2par[1][2][6][4]=-7.698077e-01; UB2par[1][2][6][5]=5.321102e-02;
   UB2par[1][2][7][0]=-7.988005e+00; UB2par[1][2][7][1]=1.623697e+01; UB2par[1][2][7][2]=-1.192795e+01; UB2par[1][2][7][3]=4.337252e+00; UB2par[1][2][7][4]=-7.698952e-01; UB2par[1][2][7][5]=5.321713e-02;
   UB2par[0][3][5][0]=1.622496e+00; UB2par[0][3][5][1]=-3.226650e+00; UB2par[0][3][5][2]=3.058418e+00; UB2par[0][3][5][3]=-1.204365e+00; UB2par[0][3][5][4]=2.218283e-01; UB2par[0][3][5][5]=-1.579729e-02;
   UB2par[0][3][6][0]=3.901885e-01; UB2par[0][3][6][1]=-1.687626e-01; UB2par[0][3][6][2]=3.333113e-01; UB2par[0][3][6][3]=-1.210109e-01; UB2par[0][3][6][4]=2.404923e-02; UB2par[0][3][6][5]=-2.239301e-03;
   UB2par[0][3][7][0]=3.901886e-01; UB2par[0][3][7][1]=-1.687628e-01; UB2par[0][3][7][2]=3.333115e-01; UB2par[0][3][7][3]=-1.210109e-01; UB2par[0][3][7][4]=2.404924e-02; UB2par[0][3][7][5]=-2.239302e-03;
   UB2par[1][3][5][0]=5.603681e-01; UB2par[1][3][5][1]=-1.802347e+00; UB2par[1][3][5][2]=2.670415e+00; UB2par[1][3][5][3]=-1.308521e+00; UB2par[1][3][5][4]=2.774597e-01; UB2par[1][3][5][5]=-2.171773e-02;
   UB2par[1][3][6][0]=9.217661e-02; UB2par[1][3][6][1]=-8.349794e-01; UB2par[1][3][6][2]=1.815357e+00; UB2par[1][3][6][3]=-9.665695e-01; UB2par[1][3][6][4]=2.154214e-01; UB2par[1][3][6][5]=-1.755899e-02;
   UB2par[1][3][7][0]=9.217667e-02; UB2par[1][3][7][1]=-8.349796e-01; UB2par[1][3][7][2]=1.815357e+00; UB2par[1][3][7][3]=-9.665695e-01; UB2par[1][3][7][4]=2.154214e-01; UB2par[1][3][7][5]=-1.755899e-02;
   UB2par[0][4][5][0]=5.221130e-01; UB2par[0][4][5][1]=-4.807964e-01; UB2par[0][4][5][2]=6.680592e-01; UB2par[0][4][5][3]=-2.657299e-01; UB2par[0][4][5][4]=5.018773e-02; UB2par[0][4][5][5]=-3.879997e-03;
   UB2par[0][4][6][0]=-3.041759e-01; UB2par[0][4][6][1]=1.310374e+00; UB2par[0][4][6][2]=-7.234780e-01; UB2par[0][4][6][3]=2.006525e-01; UB2par[0][4][6][4]=-1.831745e-02; UB2par[0][4][6][5]=-3.687132e-04;
   UB2par[0][4][7][0]=-3.041758e-01; UB2par[0][4][7][1]=1.310374e+00; UB2par[0][4][7][2]=-7.234776e-01; UB2par[0][4][7][3]=2.006523e-01; UB2par[0][4][7][4]=-1.831742e-02; UB2par[0][4][7][5]=-3.687154e-04;
   UB2par[1][4][5][0]=-3.949075e+00; UB2par[1][4][5][1]=7.695501e+00; UB2par[1][4][5][2]=-4.906640e+00; UB2par[1][4][5][3]=1.564977e+00; UB2par[1][4][5][4]=-2.451110e-01; UB2par[1][4][5][5]=1.499430e-02;
   UB2par[1][4][6][0]=-8.655066e+00; UB2par[1][4][6][1]=1.588443e+01; UB2par[1][4][6][2]=-1.047854e+01; UB2par[1][4][6][3]=3.405363e+00; UB2par[1][4][6][4]=-5.411753e-01; UB2par[1][4][6][5]=3.363546e-02;
   UB2par[1][4][7][0]=-8.655066e+00; UB2par[1][4][7][1]=1.588443e+01; UB2par[1][4][7][2]=-1.047854e+01; UB2par[1][4][7][3]=3.405363e+00; UB2par[1][4][7][4]=-5.411753e-01; UB2par[1][4][7][5]=3.363546e-02;
   UB2par[0][5][5][0]=-1.725408e-01; UB2par[0][5][5][1]=1.269309e+00; UB2par[0][5][5][2]=-7.965698e-01; UB2par[0][5][5][3]=2.788663e-01; UB2par[0][5][5][4]=-4.393288e-02; UB2par[0][5][5][5]=2.320609e-03;
   UB2par[0][5][6][0]=-8.372118e-01; UB2par[0][5][6][1]=2.433684e+00; UB2par[0][5][6][2]=-1.457151e+00; UB2par[0][5][6][3]=3.911086e-01; UB2par[0][5][6][4]=-3.674177e-02; UB2par[0][5][6][5]=-1.054216e-04;
   UB2par[0][5][7][0]=-8.372111e-01; UB2par[0][5][7][1]=2.433683e+00; UB2par[0][5][7][2]=-1.457149e+00; UB2par[0][5][7][3]=3.911083e-01; UB2par[0][5][7][4]=-3.674173e-02; UB2par[0][5][7][5]=-1.054226e-04;
   UB2par[1][5][5][0]=-1.250923e+00; UB2par[1][5][5][1]=2.023986e+00; UB2par[1][5][5][2]=-1.629962e-01; UB2par[1][5][5][3]=-3.465045e-01; UB2par[1][5][5][4]=1.230897e-01; UB2par[1][5][5][5]=-1.219968e-02;
   UB2par[1][5][6][0]=2.271053e+00; UB2par[1][5][6][1]=-6.248786e+00; UB2par[1][5][6][2]=6.912723e+00; UB2par[1][5][6][3]=-3.175855e+00; UB2par[1][5][6][4]=6.589580e-01; UB2par[1][5][6][5]=-5.103811e-02;
   UB2par[1][5][7][0]=2.271053e+00; UB2par[1][5][7][1]=-6.248786e+00; UB2par[1][5][7][2]=6.912723e+00; UB2par[1][5][7][3]=-3.175855e+00; UB2par[1][5][7][4]=6.589580e-01; UB2par[1][5][7][5]=-5.103811e-02;
   UB2par[0][6][5][0]=3.264470e-02; UB2par[0][6][5][1]=9.191193e-01; UB2par[0][6][5][2]=-5.452803e-01; UB2par[0][6][5][3]=1.810982e-01; UB2par[0][6][5][4]=-2.489962e-02; UB2par[0][6][5][5]=8.959531e-04;
   UB2par[0][6][6][0]=6.908920e-01; UB2par[0][6][6][1]=-9.596749e-01; UB2par[0][6][6][2]=1.467861e+00; UB2par[0][6][6][3]=-8.192783e-01; UB2par[0][6][6][4]=2.019615e-01; UB2par[0][6][6][5]=-1.809310e-02;
   UB2par[0][6][7][0]=6.908924e-01; UB2par[0][6][7][1]=-9.596757e-01; UB2par[0][6][7][2]=1.467861e+00; UB2par[0][6][7][3]=-8.192786e-01; UB2par[0][6][7][4]=2.019616e-01; UB2par[0][6][7][5]=-1.809310e-02;
   UB2par[1][6][5][0]=3.199551e+00; UB2par[1][6][5][1]=-7.351798e+00; UB2par[1][6][5][2]=7.295437e+00; UB2par[1][6][5][3]=-3.170493e+00; UB2par[1][6][5][4]=6.356926e-01; UB2par[1][6][5][5]=-4.812071e-02;
   UB2par[1][6][6][0]=5.957929e+00; UB2par[1][6][6][1]=-1.433144e+01; UB2par[1][6][6][2]=1.355719e+01; UB2par[1][6][6][3]=-5.747693e+00; UB2par[1][6][6][4]=1.131725e+00; UB2par[1][6][6][5]=-8.434533e-02;
   UB2par[1][6][7][0]=5.957929e+00; UB2par[1][6][7][1]=-1.433144e+01; UB2par[1][6][7][2]=1.355719e+01; UB2par[1][6][7][3]=-5.747693e+00; UB2par[1][6][7][4]=1.131725e+00; UB2par[1][6][7][5]=-8.434533e-02;
   UB2par[0][7][5][0]=-9.764313e-02; UB2par[0][7][5][1]=1.481571e+00; UB2par[0][7][5][2]=-1.074524e+00; UB2par[0][7][5][3]=3.852899e-01; UB2par[0][7][5][4]=-6.037745e-02; UB2par[0][7][5][5]=3.195313e-03;
   UB2par[0][7][6][0]=7.064832e-02; UB2par[0][7][6][1]=5.731117e-01; UB2par[0][7][6][2]=1.785462e-01; UB2par[0][7][6][3]=-3.219503e-01; UB2par[0][7][6][4]=1.113029e-01; UB2par[0][7][6][5]=-1.176575e-02;
   UB2par[0][7][7][0]=7.064838e-02; UB2par[0][7][7][1]=5.731116e-01; UB2par[0][7][7][2]=1.785463e-01; UB2par[0][7][7][3]=-3.219504e-01; UB2par[0][7][7][4]=1.113029e-01; UB2par[0][7][7][5]=-1.176575e-02;
   UB2par[1][7][5][0]=3.215142e+00; UB2par[1][7][5][1]=-8.113223e+00; UB2par[1][7][5][2]=8.372008e+00; UB2par[1][7][5][3]=-3.719578e+00; UB2par[1][7][5][4]=7.559258e-01; UB2par[1][7][5][5]=-5.772487e-02;
   UB2par[1][7][6][0]=5.834812e+00; UB2par[1][7][6][1]=-1.495214e+01; UB2par[1][7][6][2]=1.455298e+01; UB2par[1][7][6][3]=-6.253743e+00; UB2par[1][7][6][4]=1.239546e+00; UB2par[1][7][6][5]=-9.269960e-02;
   UB2par[1][7][7][0]=5.834812e+00; UB2par[1][7][7][1]=-1.495214e+01; UB2par[1][7][7][2]=1.455298e+01; UB2par[1][7][7][3]=-6.253743e+00; UB2par[1][7][7][4]=1.239546e+00; UB2par[1][7][7][5]=-9.269960e-02;
   UB2par[0][8][5][0]=-1.035234e-01; UB2par[0][8][5][1]=1.896869e+00; UB2par[0][8][5][2]=-1.536869e+00; UB2par[0][8][5][3]=5.725712e-01; UB2par[0][8][5][4]=-9.326693e-02; UB2par[0][8][5][5]=5.302340e-03;
   UB2par[0][8][6][0]=1.610861e+00; UB2par[0][8][6][1]=-2.515972e+00; UB2par[0][8][6][2]=2.694462e+00; UB2par[0][8][6][3]=-1.319535e+00; UB2par[0][8][6][4]=3.010353e-01; UB2par[0][8][6][5]=-2.560936e-02;
   UB2par[0][8][7][0]=1.610861e+00; UB2par[0][8][7][1]=-2.515972e+00; UB2par[0][8][7][2]=2.694462e+00; UB2par[0][8][7][3]=-1.319536e+00; UB2par[0][8][7][4]=3.010353e-01; UB2par[0][8][7][5]=-2.560936e-02;
   UB2par[1][8][5][0]=4.146812e+00; UB2par[1][8][5][1]=-1.050384e+01; UB2par[1][8][5][2]=1.047976e+01; UB2par[1][8][5][3]=-4.566534e+00; UB2par[1][8][5][4]=9.153892e-01; UB2par[1][8][5][5]=-6.916728e-02;
   UB2par[1][8][6][0]=6.181927e+00; UB2par[1][8][6][1]=-1.541062e+01; UB2par[1][8][6][2]=1.467646e+01; UB2par[1][8][6][3]=-6.204489e+00; UB2par[1][8][6][4]=1.212444e+00; UB2par[1][8][6][5]=-8.949377e-02;
   UB2par[1][8][7][0]=6.181927e+00; UB2par[1][8][7][1]=-1.541062e+01; UB2par[1][8][7][2]=1.467646e+01; UB2par[1][8][7][3]=-6.204489e+00; UB2par[1][8][7][4]=1.212444e+00; UB2par[1][8][7][5]=-8.949377e-02;
   UB2par[0][9][5][0]=-9.159728e-03; UB2par[0][9][5][1]=2.008523e+00; UB2par[0][9][5][2]=-1.734067e+00; UB2par[0][9][5][3]=6.672596e-01; UB2par[0][9][5][4]=-1.130124e-01; UB2par[0][9][5][5]=6.868374e-03;
   UB2par[0][9][6][0]=5.546412e-01; UB2par[0][9][6][1]=-9.228329e-03; UB2par[0][9][6][2]=5.674406e-01; UB2par[0][9][6][3]=-4.699301e-01; UB2par[0][9][6][4]=1.377968e-01; UB2par[0][9][6][5]=-1.343413e-02;
   UB2par[0][9][7][0]=5.546402e-01; UB2par[0][9][7][1]=-9.225785e-03; UB2par[0][9][7][2]=5.674384e-01; UB2par[0][9][7][3]=-4.699292e-01; UB2par[0][9][7][4]=1.377967e-01; UB2par[0][9][7][5]=-1.343412e-02;
   UB2par[1][9][5][0]=4.703566e+00; UB2par[1][9][5][1]=-1.173234e+01; UB2par[1][9][5][2]=1.144968e+01; UB2par[1][9][5][3]=-4.915648e+00; UB2par[1][9][5][4]=9.736016e-01; UB2par[1][9][5][5]=-7.280733e-02;
   UB2par[1][9][6][0]=6.537609e+00; UB2par[1][9][6][1]=-1.636631e+01; UB2par[1][9][6][2]=1.554318e+01; UB2par[1][9][6][3]=-6.544531e+00; UB2par[1][9][6][4]=1.272391e+00; UB2par[1][9][6][5]=-9.338642e-02;
   UB2par[1][9][7][0]=6.537609e+00; UB2par[1][9][7][1]=-1.636631e+01; UB2par[1][9][7][2]=1.554318e+01; UB2par[1][9][7][3]=-6.544531e+00; UB2par[1][9][7][4]=1.272391e+00; UB2par[1][9][7][5]=-9.338642e-02;
   UB2par[0][10][5][0]=-9.398683e-01; UB2par[0][10][5][1]=4.215773e+00; UB2par[0][10][5][2]=-3.491966e+00; UB2par[0][10][5][3]=1.300033e+00; UB2par[0][10][5][4]=-2.195189e-01; UB2par[0][10][5][5]=1.368453e-02;
   UB2par[0][10][6][0]=8.339200e-01; UB2par[0][10][6][1]=-3.895834e-01; UB2par[0][10][6][2]=9.101594e-01; UB2par[0][10][6][3]=-6.499607e-01; UB2par[0][10][6][4]=1.825436e-01; UB2par[0][10][6][5]=-1.750912e-02;
   UB2par[0][10][7][0]=8.339201e-01; UB2par[0][10][7][1]=-3.895836e-01; UB2par[0][10][7][2]=9.101596e-01; UB2par[0][10][7][3]=-6.499608e-01; UB2par[0][10][7][4]=1.825436e-01; UB2par[0][10][7][5]=-1.750912e-02;
   UB2par[1][10][5][0]=4.653740e+00; UB2par[1][10][5][1]=-1.121746e+01; UB2par[1][10][5][2]=1.072528e+01; UB2par[1][10][5][3]=-4.528485e+00; UB2par[1][10][5][4]=8.841796e-01; UB2par[1][10][5][5]=-6.530611e-02;
   UB2par[1][10][6][0]=5.414247e+00; UB2par[1][10][6][1]=-1.328914e+01; UB2par[1][10][6][2]=1.255841e+01; UB2par[1][10][6][3]=-5.231845e+00; UB2par[1][10][6][4]=1.005434e+00; UB2par[1][10][6][5]=-7.297877e-02;
   UB2par[1][10][7][0]=5.414247e+00; UB2par[1][10][7][1]=-1.328914e+01; UB2par[1][10][7][2]=1.255841e+01; UB2par[1][10][7][3]=-5.231845e+00; UB2par[1][10][7][4]=1.005434e+00; UB2par[1][10][7][5]=-7.297877e-02;
   UB2par[0][11][5][0]=-9.398683e-01; UB2par[0][11][5][1]=4.215773e+00; UB2par[0][11][5][2]=-3.491966e+00; UB2par[0][11][5][3]=1.300033e+00; UB2par[0][11][5][4]=-2.195189e-01; UB2par[0][11][5][5]=1.368453e-02;
   UB2par[0][11][6][0]=8.339200e-01; UB2par[0][11][6][1]=-3.895834e-01; UB2par[0][11][6][2]=9.101594e-01; UB2par[0][11][6][3]=-6.499607e-01; UB2par[0][11][6][4]=1.825436e-01; UB2par[0][11][6][5]=-1.750912e-02;
   UB2par[0][11][7][0]=8.339201e-01; UB2par[0][11][7][1]=-3.895836e-01; UB2par[0][11][7][2]=9.101596e-01; UB2par[0][11][7][3]=-6.499608e-01; UB2par[0][11][7][4]=1.825436e-01; UB2par[0][11][7][5]=-1.750912e-02;
   UB2par[1][11][5][0]=4.653740e+00; UB2par[1][11][5][1]=-1.121746e+01; UB2par[1][11][5][2]=1.072528e+01; UB2par[1][11][5][3]=-4.528485e+00; UB2par[1][11][5][4]=8.841796e-01; UB2par[1][11][5][5]=-6.530611e-02;
   UB2par[1][11][6][0]=5.414247e+00; UB2par[1][11][6][1]=-1.328914e+01; UB2par[1][11][6][2]=1.255841e+01; UB2par[1][11][6][3]=-5.231845e+00; UB2par[1][11][6][4]=1.005434e+00; UB2par[1][11][6][5]=-7.297877e-02;
   UB2par[1][11][7][0]=5.414247e+00; UB2par[1][11][7][1]=-1.328914e+01; UB2par[1][11][7][2]=1.255841e+01; UB2par[1][11][7][3]=-5.231845e+00; UB2par[1][11][7][4]=1.005434e+00; UB2par[1][11][7][5]=-7.297877e-02;
   UB2par[0][12][5][0]=-7.149392e-01; UB2par[0][12][5][1]=4.006948e+00; UB2par[0][12][5][2]=-3.367976e+00; UB2par[0][12][5][3]=1.245384e+00; UB2par[0][12][5][4]=-2.075077e-01; UB2par[0][12][5][5]=1.273396e-02;
   UB2par[0][12][6][0]=2.561168e-01; UB2par[0][12][6][1]=1.249081e+00; UB2par[0][12][6][2]=-6.257025e-01; UB2par[0][12][6][3]=1.422677e-02; UB2par[0][12][6][4]=4.593519e-02; UB2par[0][12][6][5]=-6.688203e-03;
   UB2par[0][12][7][0]=2.561168e-01; UB2par[0][12][7][1]=1.249081e+00; UB2par[0][12][7][2]=-6.257026e-01; UB2par[0][12][7][3]=1.422679e-02; UB2par[0][12][7][4]=4.593518e-02; UB2par[0][12][7][5]=-6.688202e-03;
   UB2par[1][12][5][0]=3.696944e+00; UB2par[1][12][5][1]=-8.704458e+00; UB2par[1][12][5][2]=8.335854e+00; UB2par[1][12][5][3]=-3.476895e+00; UB2par[1][12][5][4]=6.670945e-01; UB2par[1][12][5][5]=-4.832916e-02;
   UB2par[1][12][6][0]=5.135480e+00; UB2par[1][12][6][1]=-1.226542e+01; UB2par[1][12][6][2]=1.136386e+01; UB2par[1][12][6][3]=-4.623381e+00; UB2par[1][12][6][4]=8.647692e-01; UB2par[1][12][6][5]=-6.087925e-02;
   UB2par[1][12][7][0]=5.135480e+00; UB2par[1][12][7][1]=-1.226542e+01; UB2par[1][12][7][2]=1.136386e+01; UB2par[1][12][7][3]=-4.623381e+00; UB2par[1][12][7][4]=8.647692e-01; UB2par[1][12][7][5]=-6.087925e-02;
   UB2par[0][13][5][0]=-7.149392e-01; UB2par[0][13][5][1]=4.006948e+00; UB2par[0][13][5][2]=-3.367976e+00; UB2par[0][13][5][3]=1.245384e+00; UB2par[0][13][5][4]=-2.075077e-01; UB2par[0][13][5][5]=1.273396e-02;
   UB2par[0][13][6][0]=2.561168e-01; UB2par[0][13][6][1]=1.249081e+00; UB2par[0][13][6][2]=-6.257025e-01; UB2par[0][13][6][3]=1.422677e-02; UB2par[0][13][6][4]=4.593519e-02; UB2par[0][13][6][5]=-6.688203e-03;
   UB2par[0][13][7][0]=2.561168e-01; UB2par[0][13][7][1]=1.249081e+00; UB2par[0][13][7][2]=-6.257026e-01; UB2par[0][13][7][3]=1.422679e-02; UB2par[0][13][7][4]=4.593518e-02; UB2par[0][13][7][5]=-6.688202e-03;
   UB2par[1][13][5][0]=3.696944e+00; UB2par[1][13][5][1]=-8.704458e+00; UB2par[1][13][5][2]=8.335854e+00; UB2par[1][13][5][3]=-3.476895e+00; UB2par[1][13][5][4]=6.670945e-01; UB2par[1][13][5][5]=-4.832916e-02;
   UB2par[1][13][6][0]=5.135480e+00; UB2par[1][13][6][1]=-1.226542e+01; UB2par[1][13][6][2]=1.136386e+01; UB2par[1][13][6][3]=-4.623381e+00; UB2par[1][13][6][4]=8.647692e-01; UB2par[1][13][6][5]=-6.087925e-02;
   UB2par[1][13][7][0]=5.135480e+00; UB2par[1][13][7][1]=-1.226542e+01; UB2par[1][13][7][2]=1.136386e+01; UB2par[1][13][7][3]=-4.623381e+00; UB2par[1][13][7][4]=8.647692e-01; UB2par[1][13][7][5]=-6.087925e-02;
   UB2par[0][14][5][0]=-1.097285e+00; UB2par[0][14][5][1]=4.333570e+00; UB2par[0][14][5][2]=-3.225335e+00; UB2par[0][14][5][3]=1.049928e+00; UB2par[0][14][5][4]=-1.497136e-01; UB2par[0][14][5][5]=7.336870e-03;
   UB2par[0][14][6][0]=1.252353e+00; UB2par[0][14][6][1]=-1.542002e+00; UB2par[0][14][6][2]=2.256241e+00; UB2par[0][14][6][3]=-1.353782e+00; UB2par[0][14][6][4]=3.480522e-01; UB2par[0][14][6][5]=-3.173999e-02;
   UB2par[0][14][7][0]=1.252353e+00; UB2par[0][14][7][1]=-1.542002e+00; UB2par[0][14][7][2]=2.256241e+00; UB2par[0][14][7][3]=-1.353782e+00; UB2par[0][14][7][4]=3.480522e-01; UB2par[0][14][7][5]=-3.173999e-02;
   UB2par[1][14][5][0]=3.482902e+00; UB2par[1][14][5][1]=-7.664491e+00; UB2par[1][14][5][2]=7.126493e+00; UB2par[1][14][5][3]=-2.877840e+00; UB2par[1][14][5][4]=5.326174e-01; UB2par[1][14][5][5]=-3.712906e-02;
   UB2par[1][14][6][0]=4.783046e+00; UB2par[1][14][6][1]=-1.120836e+01; UB2par[1][14][6][2]=1.063453e+01; UB2par[1][14][6][3]=-4.440312e+00; UB2par[1][14][6][4]=8.531457e-01; UB2par[1][14][6][5]=-6.176752e-02;
   UB2par[1][14][7][0]=4.783046e+00; UB2par[1][14][7][1]=-1.120836e+01; UB2par[1][14][7][2]=1.063453e+01; UB2par[1][14][7][3]=-4.440312e+00; UB2par[1][14][7][4]=8.531457e-01; UB2par[1][14][7][5]=-6.176752e-02;
   UB2par[0][15][5][0]=-1.097285e+00; UB2par[0][15][5][1]=4.333570e+00; UB2par[0][15][5][2]=-3.225335e+00; UB2par[0][15][5][3]=1.049928e+00; UB2par[0][15][5][4]=-1.497136e-01; UB2par[0][15][5][5]=7.336870e-03;
   UB2par[0][15][6][0]=1.252353e+00; UB2par[0][15][6][1]=-1.542002e+00; UB2par[0][15][6][2]=2.256241e+00; UB2par[0][15][6][3]=-1.353782e+00; UB2par[0][15][6][4]=3.480522e-01; UB2par[0][15][6][5]=-3.173999e-02;
   UB2par[0][15][7][0]=1.252353e+00; UB2par[0][15][7][1]=-1.542002e+00; UB2par[0][15][7][2]=2.256241e+00; UB2par[0][15][7][3]=-1.353782e+00; UB2par[0][15][7][4]=3.480522e-01; UB2par[0][15][7][5]=-3.173999e-02;
   UB2par[1][15][5][0]=3.482902e+00; UB2par[1][15][5][1]=-7.664491e+00; UB2par[1][15][5][2]=7.126493e+00; UB2par[1][15][5][3]=-2.877840e+00; UB2par[1][15][5][4]=5.326174e-01; UB2par[1][15][5][5]=-3.712906e-02;
   UB2par[1][15][6][0]=4.783046e+00; UB2par[1][15][6][1]=-1.120836e+01; UB2par[1][15][6][2]=1.063453e+01; UB2par[1][15][6][3]=-4.440312e+00; UB2par[1][15][6][4]=8.531457e-01; UB2par[1][15][6][5]=-6.176752e-02;
   UB2par[1][15][7][0]=4.783046e+00; UB2par[1][15][7][1]=-1.120836e+01; UB2par[1][15][7][2]=1.063453e+01; UB2par[1][15][7][3]=-4.440312e+00; UB2par[1][15][7][4]=8.531457e-01; UB2par[1][15][7][5]=-6.176752e-02;
   UB2par[0][16][5][0]=-1.097285e+00; UB2par[0][16][5][1]=4.333570e+00; UB2par[0][16][5][2]=-3.225335e+00; UB2par[0][16][5][3]=1.049928e+00; UB2par[0][16][5][4]=-1.497136e-01; UB2par[0][16][5][5]=7.336870e-03;
   UB2par[0][16][6][0]=1.252353e+00; UB2par[0][16][6][1]=-1.542002e+00; UB2par[0][16][6][2]=2.256241e+00; UB2par[0][16][6][3]=-1.353782e+00; UB2par[0][16][6][4]=3.480522e-01; UB2par[0][16][6][5]=-3.173999e-02;
   UB2par[0][16][7][0]=1.252353e+00; UB2par[0][16][7][1]=-1.542002e+00; UB2par[0][16][7][2]=2.256241e+00; UB2par[0][16][7][3]=-1.353782e+00; UB2par[0][16][7][4]=3.480522e-01; UB2par[0][16][7][5]=-3.173999e-02;
   UB2par[1][16][5][0]=3.482902e+00; UB2par[1][16][5][1]=-7.664491e+00; UB2par[1][16][5][2]=7.126493e+00; UB2par[1][16][5][3]=-2.877840e+00; UB2par[1][16][5][4]=5.326174e-01; UB2par[1][16][5][5]=-3.712906e-02;
   UB2par[1][16][6][0]=4.783046e+00; UB2par[1][16][6][1]=-1.120836e+01; UB2par[1][16][6][2]=1.063453e+01; UB2par[1][16][6][3]=-4.440312e+00; UB2par[1][16][6][4]=8.531457e-01; UB2par[1][16][6][5]=-6.176752e-02;
   UB2par[1][16][7][0]=4.783046e+00; UB2par[1][16][7][1]=-1.120836e+01; UB2par[1][16][7][2]=1.063453e+01; UB2par[1][16][7][3]=-4.440312e+00; UB2par[1][16][7][4]=8.531457e-01; UB2par[1][16][7][5]=-6.176752e-02;
   UB2par[0][17][5][0]=-1.097285e+00; UB2par[0][17][5][1]=4.333570e+00; UB2par[0][17][5][2]=-3.225335e+00; UB2par[0][17][5][3]=1.049928e+00; UB2par[0][17][5][4]=-1.497136e-01; UB2par[0][17][5][5]=7.336870e-03;
   UB2par[0][17][6][0]=1.252353e+00; UB2par[0][17][6][1]=-1.542002e+00; UB2par[0][17][6][2]=2.256241e+00; UB2par[0][17][6][3]=-1.353782e+00; UB2par[0][17][6][4]=3.480522e-01; UB2par[0][17][6][5]=-3.173999e-02;
   UB2par[0][17][7][0]=1.252353e+00; UB2par[0][17][7][1]=-1.542002e+00; UB2par[0][17][7][2]=2.256241e+00; UB2par[0][17][7][3]=-1.353782e+00; UB2par[0][17][7][4]=3.480522e-01; UB2par[0][17][7][5]=-3.173999e-02;
   UB2par[1][17][5][0]=3.482902e+00; UB2par[1][17][5][1]=-7.664491e+00; UB2par[1][17][5][2]=7.126493e+00; UB2par[1][17][5][3]=-2.877840e+00; UB2par[1][17][5][4]=5.326174e-01; UB2par[1][17][5][5]=-3.712906e-02;
   UB2par[1][17][6][0]=4.783046e+00; UB2par[1][17][6][1]=-1.120836e+01; UB2par[1][17][6][2]=1.063453e+01; UB2par[1][17][6][3]=-4.440312e+00; UB2par[1][17][6][4]=8.531457e-01; UB2par[1][17][6][5]=-6.176752e-02;
   UB2par[1][17][7][0]=4.783046e+00; UB2par[1][17][7][1]=-1.120836e+01; UB2par[1][17][7][2]=1.063453e+01; UB2par[1][17][7][3]=-4.440312e+00; UB2par[1][17][7][4]=8.531457e-01; UB2par[1][17][7][5]=-6.176752e-02;

   UB2par[2][2][0][0]=9.998756e-01; UB2par[2][2][0][1]=-1.244226e-04;
   UB2par[2][2][1][0]=1.007746e+00; UB2par[2][2][1][1]=2.271169e-03;
   UB2par[2][2][2][0]=1.013753e+00; UB2par[2][2][2][1]=4.400102e-03;
   UB2par[2][2][3][0]=9.974145e-01; UB2par[2][2][3][1]=1.875531e-03;
   UB2par[2][2][4][0]=9.946176e-01; UB2par[2][2][4][1]=3.620695e-03;
   UB2par[2][2][5][0]=9.983828e-01; UB2par[2][2][5][1]=3.487181e-03;
   UB2par[2][2][6][0]=9.785508e-01; UB2par[2][2][6][1]=1.471967e-03;
   UB2par[2][2][7][0]=9.785508e-01; UB2par[2][2][7][1]=1.471967e-03;
   UB2par[2][3][0][0]=1.003201e+00; UB2par[2][3][0][1]=2.306857e-03;
   UB2par[2][3][1][0]=1.009041e+00; UB2par[2][3][1][1]=5.005484e-03;
   UB2par[2][3][2][0]=1.010480e+00; UB2par[2][3][2][1]=6.672522e-03;
   UB2par[2][3][3][0]=9.976688e-01; UB2par[2][3][3][1]=6.195976e-03;
   UB2par[2][3][4][0]=9.957765e-01; UB2par[2][3][4][1]=7.153876e-03;
   UB2par[2][3][5][0]=9.933731e-01; UB2par[2][3][5][1]=8.358109e-03;
   UB2par[2][3][6][0]=9.791636e-01; UB2par[2][3][6][1]=9.645891e-03;
   UB2par[2][3][7][0]=9.791636e-01; UB2par[2][3][7][1]=9.645891e-03;
   UB2par[2][4][0][0]=1.002612e+00; UB2par[2][4][0][1]=3.686907e-03;
   UB2par[2][4][1][0]=1.008200e+00; UB2par[2][4][1][1]=6.155735e-03;
   UB2par[2][4][2][0]=1.006751e+00; UB2par[2][4][2][1]=8.607059e-03;
   UB2par[2][4][3][0]=9.985937e-01; UB2par[2][4][3][1]=9.842822e-03;
   UB2par[2][4][4][0]=9.924803e-01; UB2par[2][4][4][1]=1.060933e-02;
   UB2par[2][4][5][0]=9.836799e-01; UB2par[2][4][5][1]=1.118353e-02;
   UB2par[2][4][6][0]=9.685052e-01; UB2par[2][4][6][1]=1.153618e-02;
   UB2par[2][4][7][0]=9.685052e-01; UB2par[2][4][7][1]=1.153618e-02;
   UB2par[2][5][0][0]=9.993778e-01; UB2par[2][5][0][1]=4.086919e-03;
   UB2par[2][5][1][0]=1.004513e+00; UB2par[2][5][1][1]=7.616465e-03;
   UB2par[2][5][2][0]=1.004512e+00; UB2par[2][5][2][1]=1.064278e-02;
   UB2par[2][5][3][0]=9.936895e-01; UB2par[2][5][3][1]=1.130739e-02;
   UB2par[2][5][4][0]=9.888657e-01; UB2par[2][5][4][1]=1.313947e-02;
   UB2par[2][5][5][0]=9.776391e-01; UB2par[2][5][5][1]=1.437895e-02;
   UB2par[2][5][6][0]=9.607296e-01; UB2par[2][5][6][1]=1.202973e-02;
   UB2par[2][5][7][0]=9.607296e-01; UB2par[2][5][7][1]=1.202973e-02;
   UB2par[2][6][0][0]=9.966637e-01; UB2par[2][6][0][1]=6.326611e-03;
   UB2par[2][6][1][0]=1.004468e+00; UB2par[2][6][1][1]=1.093060e-02;
   UB2par[2][6][2][0]=1.001327e+00; UB2par[2][6][2][1]=1.337524e-02;
   UB2par[2][6][3][0]=9.963099e-01; UB2par[2][6][3][1]=1.540435e-02;
   UB2par[2][6][4][0]=9.844940e-01; UB2par[2][6][4][1]=1.664423e-02;
   UB2par[2][6][5][0]=9.712453e-01; UB2par[2][6][5][1]=1.640343e-02;
   UB2par[2][6][6][0]=9.511229e-01; UB2par[2][6][6][1]=1.501369e-02;
   UB2par[2][6][7][0]=9.511229e-01; UB2par[2][6][7][1]=1.501369e-02;
   UB2par[2][7][0][0]=9.966637e-01; UB2par[2][7][0][1]=6.326611e-03;
   UB2par[2][7][1][0]=1.004468e+00; UB2par[2][7][1][1]=1.093060e-02;
   UB2par[2][7][2][0]=1.001327e+00; UB2par[2][7][2][1]=1.337524e-02;
   UB2par[2][7][3][0]=9.963099e-01; UB2par[2][7][3][1]=1.540435e-02;
   UB2par[2][7][4][0]=9.844940e-01; UB2par[2][7][4][1]=1.664423e-02;
   UB2par[2][7][5][0]=9.712453e-01; UB2par[2][7][5][1]=1.640343e-02;
   UB2par[2][7][6][0]=9.511229e-01; UB2par[2][7][6][1]=1.501369e-02;
   UB2par[2][7][7][0]=9.511229e-01; UB2par[2][7][7][1]=1.501369e-02;
   UB2par[2][8][0][0]=9.991571e-01; UB2par[2][8][0][1]=8.918705e-03;
   UB2par[2][8][1][0]=1.003963e+00; UB2par[2][8][1][1]=1.329114e-02;
   UB2par[2][8][2][0]=1.001717e+00; UB2par[2][8][2][1]=1.519614e-02;
   UB2par[2][8][3][0]=9.907165e-01; UB2par[2][8][3][1]=1.666835e-02;
   UB2par[2][8][4][0]=9.781283e-01; UB2par[2][8][4][1]=1.767932e-02;
   UB2par[2][8][5][0]=9.671014e-01; UB2par[2][8][5][1]=1.727924e-02;
   UB2par[2][8][6][0]=9.443637e-01; UB2par[2][8][6][1]=1.855524e-02;
   UB2par[2][8][7][0]=9.443637e-01; UB2par[2][8][7][1]=1.855525e-02;
   UB2par[2][9][0][0]=9.991571e-01; UB2par[2][9][0][1]=8.918705e-03;
   UB2par[2][9][1][0]=1.003963e+00; UB2par[2][9][1][1]=1.329114e-02;
   UB2par[2][9][2][0]=1.001717e+00; UB2par[2][9][2][1]=1.519614e-02;
   UB2par[2][9][3][0]=9.907165e-01; UB2par[2][9][3][1]=1.666835e-02;
   UB2par[2][9][4][0]=9.781283e-01; UB2par[2][9][4][1]=1.767932e-02;
   UB2par[2][9][5][0]=9.671014e-01; UB2par[2][9][5][1]=1.727924e-02;
   UB2par[2][9][6][0]=9.443637e-01; UB2par[2][9][6][1]=1.855524e-02;
   UB2par[2][9][7][0]=9.443637e-01; UB2par[2][9][7][1]=1.855525e-02;
   UB2par[2][10][0][0]=9.965135e-01; UB2par[2][10][0][1]=1.032083e-02;
   UB2par[2][10][1][0]=1.003783e+00; UB2par[2][10][1][1]=1.563989e-02;
   UB2par[2][10][2][0]=9.967430e-01; UB2par[2][10][2][1]=1.643657e-02;
   UB2par[2][10][3][0]=9.866829e-01; UB2par[2][10][3][1]=1.779283e-02;
   UB2par[2][10][4][0]=9.771008e-01; UB2par[2][10][4][1]=2.070549e-02;
   UB2par[2][10][5][0]=9.588268e-01; UB2par[2][10][5][1]=1.802305e-02;
   UB2par[2][10][6][0]=9.434549e-01; UB2par[2][10][6][1]=1.330256e-02;
   UB2par[2][10][7][0]=9.434549e-01; UB2par[2][10][7][1]=1.330256e-02;
   UB2par[2][11][0][0]=9.965135e-01; UB2par[2][11][0][1]=1.032083e-02;
   UB2par[2][11][1][0]=1.003783e+00; UB2par[2][11][1][1]=1.563989e-02;
   UB2par[2][11][2][0]=9.967430e-01; UB2par[2][11][2][1]=1.643657e-02;
   UB2par[2][11][3][0]=9.866829e-01; UB2par[2][11][3][1]=1.779283e-02;
   UB2par[2][11][4][0]=9.771008e-01; UB2par[2][11][4][1]=2.070549e-02;
   UB2par[2][11][5][0]=9.588268e-01; UB2par[2][11][5][1]=1.802305e-02;
   UB2par[2][11][6][0]=9.434549e-01; UB2par[2][11][6][1]=1.330256e-02;
   UB2par[2][11][7][0]=9.434549e-01; UB2par[2][11][7][1]=1.330256e-02;
   UB2par[2][12][0][0]=9.961386e-01; UB2par[2][12][0][1]=1.147354e-02;
   UB2par[2][12][1][0]=1.005906e+00; UB2par[2][12][1][1]=1.838539e-02;
   UB2par[2][12][2][0]=9.936782e-01; UB2par[2][12][2][1]=1.872676e-02;
   UB2par[2][12][3][0]=9.813275e-01; UB2par[2][12][3][1]=2.014975e-02;
   UB2par[2][12][4][0]=9.765661e-01; UB2par[2][12][4][1]=2.379850e-02;
   UB2par[2][12][5][0]=9.532637e-01; UB2par[2][12][5][1]=2.107804e-02;
   UB2par[2][12][6][0]=9.465004e-01; UB2par[2][12][6][1]=-1.097671e-02;
   UB2par[2][12][7][0]=9.465004e-01; UB2par[2][12][7][1]=-1.097671e-02;
   UB2par[2][13][0][0]=9.961386e-01; UB2par[2][13][0][1]=1.147354e-02;
   UB2par[2][13][1][0]=1.005906e+00; UB2par[2][13][1][1]=1.838539e-02;
   UB2par[2][13][2][0]=9.936782e-01; UB2par[2][13][2][1]=1.872676e-02;
   UB2par[2][13][3][0]=9.813275e-01; UB2par[2][13][3][1]=2.014975e-02;
   UB2par[2][13][4][0]=9.765661e-01; UB2par[2][13][4][1]=2.379850e-02;
   UB2par[2][13][5][0]=9.532637e-01; UB2par[2][13][5][1]=2.107804e-02;
   UB2par[2][13][6][0]=9.465004e-01; UB2par[2][13][6][1]=-1.097671e-02;
   UB2par[2][13][7][0]=9.465004e-01; UB2par[2][13][7][1]=-1.097671e-02;
   UB2par[2][14][0][0]=9.954464e-01; UB2par[2][14][0][1]=1.279427e-02;
   UB2par[2][14][1][0]=1.001686e+00; UB2par[2][14][1][1]=1.989297e-02;
   UB2par[2][14][2][0]=9.961128e-01; UB2par[2][14][2][1]=2.135345e-02;
   UB2par[2][14][3][0]=9.818164e-01; UB2par[2][14][3][1]=2.337861e-02;
   UB2par[2][14][4][0]=9.721932e-01; UB2par[2][14][4][1]=2.652124e-02;
   UB2par[2][14][5][0]=9.507111e-01; UB2par[2][14][5][1]=2.178727e-02;
   UB2par[2][14][6][0]=9.787320e-01; UB2par[2][14][6][1]=6.754815e-03;
   UB2par[2][14][7][0]=9.787320e-01; UB2par[2][14][7][1]=6.754815e-03;
   UB2par[2][15][0][0]=9.954464e-01; UB2par[2][15][0][1]=1.279427e-02;
   UB2par[2][15][1][0]=1.001686e+00; UB2par[2][15][1][1]=1.989297e-02;
   UB2par[2][15][2][0]=9.961128e-01; UB2par[2][15][2][1]=2.135345e-02;
   UB2par[2][15][3][0]=9.818164e-01; UB2par[2][15][3][1]=2.337861e-02;
   UB2par[2][15][4][0]=9.721932e-01; UB2par[2][15][4][1]=2.652124e-02;
   UB2par[2][15][5][0]=9.507111e-01; UB2par[2][15][5][1]=2.178727e-02;
   UB2par[2][15][6][0]=9.787320e-01; UB2par[2][15][6][1]=6.754815e-03;
   UB2par[2][15][7][0]=9.787320e-01; UB2par[2][15][7][1]=6.754815e-03;
   UB2par[2][16][0][0]=9.954889e-01; UB2par[2][16][0][1]=1.469025e-02;
   UB2par[2][16][1][0]=9.962365e-01; UB2par[2][16][1][1]=1.995037e-02;
   UB2par[2][16][2][0]=9.882333e-01; UB2par[2][16][2][1]=2.219731e-02;
   UB2par[2][16][3][0]=9.770941e-01; UB2par[2][16][3][1]=2.564271e-02;
   UB2par[2][16][4][0]=9.608747e-01; UB2par[2][16][4][1]=2.699806e-02;
   UB2par[2][16][5][0]=9.492431e-01; UB2par[2][16][5][1]=2.289600e-02;
   UB2par[2][16][6][0]=9.827603e-01; UB2par[2][16][6][1]=4.128575e-03;
   UB2par[2][16][7][0]=9.827603e-01; UB2par[2][16][7][1]=4.128575e-03;
   UB2par[2][17][0][0]=9.954889e-01; UB2par[2][17][0][1]=1.469025e-02;
   UB2par[2][17][1][0]=9.962365e-01; UB2par[2][17][1][1]=1.995037e-02;
   UB2par[2][17][2][0]=9.882333e-01; UB2par[2][17][2][1]=2.219731e-02;
   UB2par[2][17][3][0]=9.770941e-01; UB2par[2][17][3][1]=2.564271e-02;
   UB2par[2][17][4][0]=9.608747e-01; UB2par[2][17][4][1]=2.699806e-02;
   UB2par[2][17][5][0]=9.492431e-01; UB2par[2][17][5][1]=2.289600e-02;
   UB2par[2][17][6][0]=9.827603e-01; UB2par[2][17][6][1]=4.128575e-03;
   UB2par[2][17][7][0]=9.827603e-01; UB2par[2][17][7][1]=4.128575e-03;

    
    // load up the map

    addItem("EvtRun",           &EvtRun, true);
    addItem("EvtEventId",       &EvtEventId, true);
    addItem("EvtEventId64",     &EvtEventId64, true);
    addItem("EvtElapsedTime",   &EvtElapsedTime, true);
    addItem("EvtLiveTime",      &EvtLiveTime, true);

    addItem("EvtEnergyCorr",    &EvtEnergyCorr,   true);
    addItem("EvtEnergyCorrUB",  &EvtEnergyCorrUB, true);
    addItem("EvtEnergyRaw",     &EvtEnergyRaw);
    //addItem("EvtDeltaEoE",      &EvtDeltaEoE);  // moved to McValsTool
    addItem("EvtCalEdgeAngle",  &EvtCalEdgeAngle);
    addItem("EvtTkrEdgeAngle",  &EvtTkrEdgeAngle);
    addItem("EvtLogEnergy",     &EvtLogEnergy);
    addItem("EvtTkr1EFrac",     &EvtTkr1EFrac, true);
    addItem("EvtVtxKin",        &EvtVtxKin);
    addItem("EvtVtxEAngle",     &EvtVtxEAngle);
    addItem("EvtTkrComptonRatio", &EvtTkrComptonRatio);
        addItem("EvtETkrComptonRatio", &EvtETkrComptonRatio);
    addItem("EvtPSFModel",      &EvtPSFModel);

    addItem("EvtETkr1Chisq",    &EvtETkr1Chisq);
    addItem("EvtETkr1FirstChisq", &EvtETkr1FirstChisq);
    addItem("EvtETkr1Qual",     &EvtETkr1Qual);
    addItem("EvtTkr1PSFMdRat",  &EvtTkr1PSFMdRat);

    addItem("EvtECalXtalRatio", &EvtECalXtalRatio);
    addItem("EvtECalXtalTrunc", &EvtECalXtalTrunc);
    addItem("EvtECalTrackDoca", &EvtECalTrackDoca);
    //addItem("EvtECalTrackSep",  &EvtECalTrackSep);
        addItem("EvtECalTransRms",  &EvtECalTransRms);
    addItem("EvtECalLongRms",   &EvtECalLongRms);
        addItem("EvtECalLRmsAsym",  &EvtECalLRmsAsym);
        addItem("EvtECalTrackAngle",&EvtECalTrackAngle);

    addItem("EvtEVtxAngle",     &EvtEVtxAngle);
    addItem("EvtEVtxDoca",      &EvtEVtxDoca);
    addItem("EvtPsf68",         &EvtPsf68);
    addItem("EvtEventFlags",      &EvtEventFlags);

    addItem("EvtEnergyCorrUB2",&EvtEnergyCorrUB2);
    addItem("NewEvtEnergyCorr",&NewEvtEnergyCorr);
    addItem("NewEvtEnergyCorrUB2",&NewEvtEnergyCorrUB2);

    addItem("CalNewCfpEnergyUB2",&CalNewCfpEnergyUB2);
 
    addItem("EvtJointEnergy",&EvtJointEnergy, true);
    addItem("EvtJointLogEnergy",&EvtJointLogEnergy, true);
    addItem("EvtJointWeight",&EvtJointWeight, true);

    addItem("EvtCalCsIRLn",&EvtCalCsIRLn, true);

    zeroVals();

    //m_ubInterpolateTool->addBiasMap("Param","$(ANALYSISNTUPLEDATAPATH)/BiasMapEvtEnergyCorr.txt");
    m_ubInterpolateTool->addBiasMap("ParamFront","$(CALUTILXMLPATH)/BiasMapEvtEnergyCorr_Front.txt");
    m_ubInterpolateTool->addBiasMap("ParamBack","$(CALUTILXMLPATH)/BiasMapEvtEnergyCorr_Back.txt"); 

    return sc;
}

void EvtValsTool::zeroVals() 
{
    ValBase::zeroVals();
    // restore some interesting event variables, no need to touch them again...
    fillHeaderInfo();
}

StatusCode EvtValsTool::calculate()
{
  StatusCode sc = StatusCode::SUCCESS;
  
  MsgStream log(msgSvc(), name());
  
  // we may want to add TDS stuff to this method, but we haven't needed it yet.
  
  //int firstCheck = m_check;
  int nextCheck = NOCALC;
  
  
  // since we know what's happening here, we can plan a little
  // the idea is to call the first check of each tool with the called check value,
  // and the rest with the no-calc value
  // so be careful when adding calls or moving stuff around!!!!
  
  // So far, we're only using Tkr, Cal, Vtx and Mc quantities. 
  // The corresponding ValsTools have already been called, 
  // so we can use nextCheck throughout. 
  
  float eCalSum = -1.0, eTkr = -1.0;
  if(    m_pCalTool->getVal("CalEnergyRaw", eCalSum, nextCheck).isSuccess()
         && m_pTkrTool->getVal("TkrEnergyCorr", eTkr, nextCheck).isSuccess()) {
    EvtEnergyRaw = eTkr + eCalSum;
  }
  
  float eTkrKalEne, eCalRLn, eTkrBest; //, eCal
  int CAL_Type;
  if (   m_pCalTool->getVal("CalCsIRLn", eCalRLn, nextCheck).isSuccess()
         && m_pTkrTool->getVal("TkrSumKalEne", eTkrKalEne, nextCheck).isSuccess()) 
    {
      eTkrBest = std::max(eTkr+eCalSum, eTkrKalEne);
      if(eCalSum < 350 || eCalRLn < 4) {
        if(eCalSum < 5 || eCalRLn < 4) {
          CAL_Type = 0;
        }
        else {CAL_Type = 1;}
      }
      else {CAL_Type = 2;} 
    }
  
  // EdgeAngle:   distance from tower edge / cos(theta) 
  float tkrEdge, calEdge;
  double tkr1ZDir = -1.;
  if(m_pTkrTool->getVal("Tkr1ZDir",tkr1ZDir, nextCheck).isSuccess()) {
    if(tkr1ZDir == 0.) tkr1ZDir = -1.; 
    if (m_pTkrTool->getVal("TkrTwrEdge", tkrEdge, nextCheck).isSuccess()) {
      EvtTkrEdgeAngle = -tkrEdge/tkr1ZDir;
    }
    if (m_pCalTool->getVal("CalTwrEdge", calEdge, nextCheck).isSuccess()) {
      EvtCalEdgeAngle = -calEdge/tkr1ZDir;
    }
  }
  double tkr1ZDir2 = tkr1ZDir*tkr1ZDir;  
  // Final Energy Selection: (Combo approach) 
  //    Use Kalman energy when not enough in the Cal,
  //    use Last Layer Correction when avaialable (eCalEneLLCorr > 0) or else
  //    use integrated edge and leakage correction from CalValsTool
  float eCalSumCorr = -1.0;
  if(    m_pCalTool->getVal("CalEnergyCorr", eCalSumCorr, nextCheck).isSuccess()
         //&& m_pCalTool->getVal("CalEnergyLLCorr", eCalEneLLCorr, nextCheck).isSuccess()
         )
    {       // NOTE: IT IS STILL LARGELY UNDECIDED WHAT TO DO HERE....!!!!!!!  
            //    if(CAL_Type == 0) EvtEnergyCorr = eTkrBest;
            //        else {
            //        if(eCalEneLLCorr > 0) EvtEnergyCorr = eCalEneLLCorr; //This is wrong! LL has comp.for TkrEne!
            //        else { 
      EvtEnergyCorr = (eTkr + eCalSumCorr);
      if ( EvtEnergyCorr == 0 )
        EvtEnergyCorrUB = 0;
      else {
        float Tkr1FirstLayer, bias;
        if (m_pTkrTool->getVal("Tkr1FirstLayer", Tkr1FirstLayer, nextCheck).isSuccess()) {
          if (Tkr1FirstLayer>5)  bias = m_ubInterpolateTool->interpolate("ParamFront", log10(EvtEnergyCorr), tkr1ZDir);
          else bias = m_ubInterpolateTool->interpolate("ParamBack", log10(EvtEnergyCorr), tkr1ZDir);
        }
        else {bias=0;} 
        
        EvtEnergyCorrUB = bias == 0 ? -1 : EvtEnergyCorr / bias;
        //            std::cout << "EvtValsTool EvtEnergyCorr " << EvtEnergyCorr << " ( " << log10(EvtEnergyCorr) << " ) " << EvtEnergyCorrUB << " ( " << log10(EvtEnergyCorrUB) << " ) " << tkr1ZDir << std::endl;
      }
      //        }
      //        }
    }

  // this variable moved to McValsTools and renamed McDeltaEoE 16-May-2012 LSR
  //float mcEnergy;
  //EvtDeltaEoE = -2.0;
  
  //if(m_pMcTool!=NULL&&m_pMcTool->isLoaded()) {
  //    if(m_pMcTool->getVal("McEnergy", mcEnergy, nextCheck).isSuccess()){
  //        if (mcEnergy>0) { 
  //            EvtDeltaEoE = (EvtEnergyCorr - mcEnergy)/(mcEnergy);
  //        }
  //    } 
  //}
  
  // Model simple for PSF(68%) 
  EvtPSFModel = sqrt(pow((.061/pow((std::max(EvtEnergyCorr*1.,1.)/100),.8)),2) + (.001745*.001745));
  m_pPsfTool->setEstimate(EvtPSFModel*180./M_PI);

  // Ph.Bruel: Increase maximum energy from 1TeV to 10TteV
  // Log(base 10) of measured energy - useful for parameterizing effects
  EvtLogEnergy = log10(std::min(std::max(EvtEnergyCorr,10.f),10000000.f));
  float logE = std::min(std::max(EvtLogEnergy,1.3f), 6.0f);
  float logE2 = logE*logE; 
  
  
  // Ratio of PSF model to event total PSF error. Note PhiErr = sin(theta)*d(phi)
  float tkr1ThetaErr, tkr1PhiErr;
  if (m_pTkrTool->getVal("Tkr1ThetaErr",tkr1ThetaErr, nextCheck).isSuccess() &&
      m_pTkrTool->getVal("Tkr1PhiErr",tkr1PhiErr, nextCheck).isSuccess()) {
    EvtTkr1PSFMdRat = sqrt(std::max(0.0f, tkr1ThetaErr*tkr1ThetaErr + tkr1PhiErr*tkr1PhiErr))/ 
      EvtPSFModel;
  }
  
  // Fraction of energy in Track 1
  float tkr1ConE;
  if (m_pTkrTool->getVal("Tkr1ConEne",tkr1ConE, nextCheck).isSuccess()) {
    if(EvtEnergyCorr>0.0) EvtTkr1EFrac = std::max(0.01f, tkr1ConE/EvtEnergyCorr);
  }
  
  // Vtx kinematic variable:  angle * event energy / Track_1 energy fraction
  float vtxAngle;
  if (m_pVtxTool->getVal("VtxAngle", vtxAngle, nextCheck).isSuccess()) {
    if (EvtTkr1EFrac>0.0) EvtVtxKin = vtxAngle*EvtEnergyCorr/EvtTkr1EFrac;
  }
  
  // Vtx angle x event energy  ~ constant
  EvtVtxEAngle = vtxAngle*EvtEnergyCorr;
  
  // Hit counting around track compared to hits on track.  Compton scatters
  // have a ratio ~< 1 (as do MIPs etc.  - this is similar to the former
  // Surplus Hit Ratio
  float totHits, tkr1First;
  if (m_pTkrTool->getVal("TkrTotalHits", totHits, nextCheck).isSuccess()) {
    if (m_pTkrTool->getVal("Tkr1FirstLayer", tkr1First, nextCheck).isSuccess()){
      EvtTkrComptonRatio = totHits/(2.*(m_tkrGeom->numLayers()-tkr1First));
      EvtETkrComptonRatio = EvtTkrComptonRatio/(-9.34+7.22*logE - .672*logE2)
        /(-.509 - 3.65*tkr1ZDir - 2.03*tkr1ZDir2);
    }
  }
  
  // Energy compensated track 1 chisq
  float tkr1Chisq;
  if (m_pTkrTool->getVal("Tkr1Chisq", tkr1Chisq, nextCheck).isSuccess()) {
    EvtETkr1Chisq = tkr1Chisq/(6.49 -22.1*tkr1ZDir -22.8*tkr1ZDir2);
  }
  float tkr1_1stChisq;
  if (m_pTkrTool->getVal("Tkr1FirstChisq", tkr1_1stChisq, nextCheck).isSuccess()) {
    EvtETkr1FirstChisq = tkr1_1stChisq/(4.34-.708*logE+.19*logE2)
      /(.751-1.74*tkr1ZDir-1.83*tkr1ZDir2);
  }
  float tkr1Qual;
  if (m_pTkrTool->getVal("Tkr1Qual", tkr1Qual, nextCheck).isSuccess()) {
    EvtETkr1Qual = tkr1Qual/(56.6 +6.78*tkr1ZDir + 8.23*tkr1ZDir2);
  }
  
  float calXtalRatio;
  if(m_pCalTool->getVal("CalXtalRatio", calXtalRatio, nextCheck).isSuccess()) {
    EvtECalXtalRatio = calXtalRatio/(2.99-1.19*logE  + .122*logE2)/(.749-.355*tkr1ZDir);
  }
  
  float calXtalsTrunc;
  if(m_pCalTool->getVal("CalXtalsTrunc", calXtalsTrunc, nextCheck).isSuccess()) {
    float term1 = std::max(1., (-85.4 + 65.2*logE - 10.4*logE2));
    float logE33 = logE-3.3; 
    float term2 = (logE33 > 0.) ? std::min(14., 12.*logE33*logE33): 0.;
    EvtECalXtalTrunc = calXtalsTrunc/(term1 + term2)/(.935 - .382*tkr1ZDir - .343*tkr1ZDir2);
  }
  
  float calTrackDoca;
  if(m_pCalTool->getVal("CalTrackDoca", calTrackDoca, nextCheck).isSuccess()) {
    float logETrunc = std::min(EvtLogEnergy, 3.8f); 
    EvtECalTrackDoca = calTrackDoca/(272.-140.5*logETrunc + 18.7*logETrunc*logETrunc)
      /(3.08+2.67*tkr1ZDir);
  }
  
  float calTransRms;
  if(m_pCalTool->getVal("CalTransRms", calTransRms, nextCheck).isSuccess()) {
    float logEGaussSq = (EvtLogEnergy-2.3)*(EvtLogEnergy-2.3);
    EvtECalTransRms = calTransRms/(12.5+ 22.*exp(-logEGaussSq/.8))/(1.34+.55*tkr1ZDir);
  }
  
  float calLongRms;
  if(m_pCalTool->getVal("CalLongRms", calLongRms, nextCheck).isSuccess()) {
    if(logE < 3.8) {
      EvtECalLongRms = calLongRms/(96.1 - 39*logE + 5.6*logE2);
    }
    else {
      EvtECalLongRms = calLongRms/(28.*(.847-.0134*logE+.015*logE2));
    }
    EvtECalLongRms /= (1.06+.0867*tkr1ZDir);
  }
  
  float calLRmsAsym;
  if(m_pCalTool->getVal("CalLRmsAsym", calLRmsAsym, nextCheck).isSuccess()) {
    float logEGaussSq = (logE-2.0)*(logE-2.0); 
    EvtECalLRmsAsym         = calLRmsAsym/(.012+.06*exp(-logEGaussSq/1.3))/(1.01-.0718*tkr1ZDir);
  }
  
  float calTrackAngle;
  if(m_pCalTool->getVal("CalTrackAngle", calTrackAngle, nextCheck).isSuccess()) {
    float logETrunc = std::min(logE, 4.0f); 
    EvtECalTrackAngle =        calTrackAngle/(2.3-1.11*logETrunc +.138*logETrunc*logETrunc)
      /(1.43+.612*tkr1ZDir);
  }
  
  EvtEVtxAngle = vtxAngle*sqrt(EvtEnergyCorr)/(4.24 -1.98*logE + .269*logE2)
    /(1.95+2.36*tkr1ZDir+1.3*tkr1ZDir2);
  
  float vtxDoca;
  if (m_pVtxTool->getVal("VtxDOCA", vtxDoca, nextCheck).isSuccess()) {
    EvtEVtxDoca = vtxDoca/(1.55 - .685*logE+ .0851*logE2) 
      / (2.21 + 3.01*tkr1ZDir + 1.59*tkr1ZDir2);
  }

  float myfTkr1FirstLayer = 0;
  if(!(m_pTkrTool->getVal("Tkr1FirstLayer",myfTkr1FirstLayer,nextCheck).isSuccess())) myfTkr1FirstLayer = 0;
  int myTkr1FirstLayer = (int)myfTkr1FirstLayer;
  double myTkr1XDir = 0;
  if(!(m_pTkrTool->getVal("Tkr1XDir",myTkr1XDir,nextCheck).isSuccess())) myTkr1XDir = 0;
  double myTkr1YDir = 0;
  if(!(m_pTkrTool->getVal("Tkr1YDir",myTkr1YDir,nextCheck).isSuccess())) myTkr1YDir = 0;
  double myTkr1ZDir = 0;
  if(!(m_pTkrTool->getVal("Tkr1ZDir",myTkr1ZDir,nextCheck).isSuccess())) myTkr1ZDir = 0;
  float myTkr1X0 = 0;
  if(!(m_pTkrTool->getVal("Tkr1X0",myTkr1X0,nextCheck).isSuccess())) myTkr1X0 = 0;
  float myTkr1Y0 = 0;
  if(!(m_pTkrTool->getVal("Tkr1Y0",myTkr1Y0,nextCheck).isSuccess())) myTkr1Y0 = 0;
  float myTkr1Z0 = 0;
  if(!(m_pTkrTool->getVal("Tkr1Z0",myTkr1Z0,nextCheck).isSuccess())) myTkr1Z0 = 0;
  float myTkr1StripsEnergyCorr = 0;
  if(!(m_pTkrTool->getVal("Tkr1StripsEnergyCorr",myTkr1StripsEnergyCorr,nextCheck).isSuccess())) myTkr1StripsEnergyCorr = 0;

  float myCalEnergyRaw = 0;
  if(!(m_pCalTool->getVal("CalEnergyRaw",myCalEnergyRaw, nextCheck).isSuccess())) myCalEnergyRaw = 0;
  float myCalEnergyCorr = 0;
  if(!(m_pCalTool->getVal("CalEnergyCorr",myCalEnergyCorr, nextCheck).isSuccess())) myCalEnergyCorr = 0;
  float myCalNewCfpEnergy = 0;
  if(!(m_pCalTool->getVal("CalNewCfpEnergy",myCalNewCfpEnergy, nextCheck).isSuccess())) myCalNewCfpEnergy = 0;
  float myCalNewCfpEnergyUB = 0;
  if(!(m_pCalTool->getVal("CalNewCfpEnergyUB",myCalNewCfpEnergyUB, nextCheck).isSuccess())) myCalNewCfpEnergyUB = 0;
  float myCalNewCfpCalEnergy = 0;
  if(!(m_pCalTool->getVal("CalNewCfpCalEnergy",myCalNewCfpCalEnergy, nextCheck).isSuccess())) myCalNewCfpCalEnergy = 0;
  float myCalNewCfpCalEnergyUB = 0;
  if(!(m_pCalTool->getVal("CalNewCfpCalEnergyUB",myCalNewCfpCalEnergyUB, nextCheck).isSuccess())) myCalNewCfpCalEnergyUB = 0;
  float myCal1MomZCntr = 0;
  if(!(m_pCalTool->getVal("Cal1MomZCntr",myCal1MomZCntr, nextCheck).isSuccess())) myCal1MomZCntr = 0;
  float myCalCsIRLn = 0;
  if(!(m_pCalTool->getVal("CalCsIRLn",myCalCsIRLn, nextCheck).isSuccess())) myCalCsIRLn = 0;

  // Switch to non unbiased energies in case the unbiased energies = 0
  if(myCalNewCfpEnergyUB==0 && myCalNewCfpEnergy>0) myCalNewCfpEnergyUB = myCalNewCfpEnergy;
  if(myCalNewCfpCalEnergyUB==0 && myCalNewCfpCalEnergy>0) myCalNewCfpCalEnergyUB = myCalNewCfpCalEnergy;

  // Switch to CalEnergyRaw in case CalEnergyCorr==0
  if(myCalEnergyCorr<=0) myCalEnergyCorr = myCalEnergyRaw;

  NewEvtEnergyCorr = (myCalEnergyCorr+myTkr1StripsEnergyCorr)/0.85;
  double corfactor = GetEnergyUB2Correction(0,myTkr1FirstLayer,(float)myTkr1ZDir,NewEvtEnergyCorr);
  NewEvtEnergyCorrUB2 = NewEvtEnergyCorr*corfactor;

  double myEvtEnergyCorr = EvtEnergyCorr/0.9;
  corfactor = GetEnergyUB2Correction(1,myTkr1FirstLayer,(float)myTkr1ZDir,myEvtEnergyCorr);
  EvtEnergyCorrUB2 = myEvtEnergyCorr*corfactor;

  // Unbias myCalNewcfpEnergy
  corfactor = GetEnergyUB2Correction(2,myTkr1FirstLayer,(float)myTkr1ZDir,myCalNewCfpEnergy);
  CalNewCfpEnergyUB2 = myCalNewCfpEnergy*corfactor;

  // Building the event energy from NewEvtEnergyCorrUB2 and myCalNewCfpEnergyUB as shown in Ph.Bruel presentation Energy reconstruction (slide 12) at the pass8 workshop in Washington (2012)
  // If no track, switch to myCalNewCfpCalEnergyUB if >0
  // The output of this can be slightly different from what TkrUtil::TkrEnergyTool::getEvtEnergyEstimation returns because Tkr1ZDir, CalEnergyCorr can be different (and because in TkrUtil we don't have CalNewCfpEnergyUB, but it's minor).
  double logeraw,logerawmin,logerawmax,mytkr1zdir;
  double logerawdelta = 0.1;
  logeraw = -10;
  if(myCalEnergyRaw>0) logeraw = log10(myCalEnergyRaw);

  double E0 = NewEvtEnergyCorrUB2;
  double E1 = CalNewCfpEnergyUB2;
  double weight = 0;

  EvtJointEnergy = E0;
  if(myTkr1ZDir>=0 && myCalNewCfpCalEnergyUB>0)
    {
      // if no track, take the new profile energy estimated with the cal axis
      EvtJointEnergy = myCalNewCfpCalEnergyUB;
    }
  else if(myTkr1ZDir<0 && E1>0)
    {
      // if there is a track and the new profile energy is available, take the weighted sum (depending on the raw energy)
      // around the region (in the logEraw,ZDir plane) where the parametric and profile methods have equivalent resolution
      mytkr1zdir = myTkr1ZDir;
      if(mytkr1zdir>-0.2) mytkr1zdir = -0.2;
      logerawmin = 2.85892e+00+2.16216e-02*myTkr1FirstLayer;
      if(mytkr1zdir>-0.65) logerawmin += 1.78271*(mytkr1zdir+0.65);
      logerawmin -= logerawdelta;
      logerawmax = logerawmin+2.*logerawdelta;
      //
      if(logeraw>logerawmax)
        {
          EvtJointEnergy = E1;
          weight = 1;
        }
      else if(logeraw<logerawmin)
        EvtJointEnergy = E0;
      else
        {
          weight = (logeraw-logerawmin)/2./logerawdelta;
          EvtJointEnergy = E0*(1-weight)+E1*weight;
        }
    }

  EvtJointLogEnergy = 0;
  if(EvtJointEnergy>0) EvtJointLogEnergy = log10(EvtJointEnergy);
  EvtJointWeight = weight;

  // Compute CalCsIRLn based on Tkr1
  EvtCalCsIRLn = 0;
  if(myTkr1ZDir<0)
    {
      Point x0(myTkr1X0,myTkr1Y0,myTkr1Z0);
      Vector t0(-myTkr1XDir,-myTkr1YDir,-myTkr1ZDir);
      // Construct Event Axis along which the shower will be evaluated 
      Ray axis(x0, t0); 
      double arc_len = (x0.z()- m_calZTop)/t0.z(); 
      Point cal_top = axis.position(-arc_len);   // Event axis entry point to top of Cal Stack 
      double rm_hard   = 40. + 36.*cvct_cal_trans((myCalEnergyRaw-250)/200.)* arc_len/500.;
      EvtCalCsIRLn = EvtValsTool::aveRadLens(cal_top,-t0,rm_hard/4.,6,myCal1MomZCntr);
    }

  //Compute the expected Psf68, based on Tkr1Theta, Tkr1FirstLayer, and EvtJointEnergy
  EvtPsf68=0.;
  double cl_level=0.68;
  float tkrNumTracks = (m_pTkrTool->getVal("TkrNumTracks",tkrNumTracks, nextCheck).isSuccess())?tkrNumTracks:-1.0;
  if(tkrNumTracks>0) {
    float tkr1Theta = (m_pTkrTool->getVal("Tkr1Theta",tkr1Theta, nextCheck).isSuccess())?tkr1Theta:-1.0;
    bool isFront = (myfTkr1FirstLayer>=m_tkrGeom->numNoConverter()+m_tkrGeom->numSuperGlast())?true:false;
    fillPsfInfo(EvtJointEnergy, tkr1Theta*180./M_PI, isFront, cl_level);
  }
  return sc;
}

double EvtValsTool::GetEnergyUB2Correction(int method, int tkr1firstlayer, double tkr1zdir, double energy)
{
  if(energy<=0) return 1;
  if(method<0||method>=3) return 1;
  if(tkr1firstlayer<2||tkr1firstlayer>=18) return 1;

  double myloge = log10(energy);
  if(myloge<UB2logemin[method]) myloge = UB2logemin[method];
  else if(myloge>UB2logemax[method]) myloge = UB2logemax[method];

  double corfactor = 1;
  double emean = 1;

  int i,j;
  double loge;
  for(i=0;i<UB2zdirn;++i)
    {
      if(method<2)
        {
          loge = 1;
          UB2val[i] = 0;
          for(j=0;j<6;++j)
            {
              UB2val[i] += UB2par[method][tkr1firstlayer][i][j]*loge;
              loge *= myloge;
            }
        }
      else
        UB2val[i] = UB2par[method][tkr1firstlayer][i][0]-UB2par[method][tkr1firstlayer][i][1]*(myloge-5.)*(myloge-5.);
    }
  UB2val[UB2zdirn] = 1;
  
  i = (int)floor((tkr1zdir-UB2zdir[0])/UB2zdirbinwidth);

  if(i<0)
    emean = UB2val[0];
  else if(i>=UB2zdirn)
    emean = UB2val[UB2zdirn];
  else
    emean = (UB2val[i]*(UB2zdir[i+1]-tkr1zdir)+UB2val[i+1]*(tkr1zdir-UB2zdir[i]))/UB2zdirbinwidth;

  if(emean>0)
    {
      corfactor = 1/emean;
      if(corfactor>3) corfactor = 3;
    }

  return corfactor;
}

void EvtValsTool::fillHeaderInfo () 
{
   SmartDataPtr<Event::EventHeader> header(m_pEventSvc, EventModel::EventHeader);

    if(header) {
        EvtRun         = header->run();
        EvtEventId     = header->event();
        EvtEventId64   = header->event();
        EvtElapsedTime = header->time();
        EvtLiveTime    = header->livetime();
        EvtEventFlags  = header->gleamEventFlags();
    }
}

void EvtValsTool::fillPsfInfo(double energy, double theta, bool isFront,  double cl_level)
{
  if(energy<=0) {
    EvtPsf68 = 0; 
    return;}
  MsgStream log(msgSvc(), name());
  EvtPsf68 = m_pPsfTool->computePsf(cl_level, energy, theta, isFront);
  log<<MSG::DEBUG<<"energy:"<<" "<< energy<<"Theta(deg):"<< theta<<" Front? "<< isFront<<" PSF68: "<<EvtPsf68<<endreq;
}

double EvtValsTool::aveRadLens(Point cal_top, Vector t0, double radius, int numSamples, double calcentroidposz)
{ // This method finds the averages and rms for a cylinder of rays passing through 
  // the calorimeter of the radiation lengths in CsI and other material. 
  // The radius of the cylinder is "radius" and the number of rays = numSample (plus the 
  //       the central ray).  
  // Note: this method as to called in sequence.  It depends on various internal variables 
  //       having already been calculated and dumps most of its output to internal vars.
  
  // Initialize the internal transfer variables
  m_radLen_CsI   = 0.;
  m_rms_RL_CsI   = 0.;
  m_radLen_Stuff = 0.;
  m_rms_RL_Stuff = 0.;
  m_radLen_Cntr  = 0.;
  m_rms_RL_Cntr  = 0.; 
  m_radLen_CntrStuff = 0.;
  m_rms_RL_CntrStuff = 0.;
  
  m_arcLen_CsI   = 0.; 
  m_arcLen_Stuff = 0.; 
  m_arcLen_Cntr  = 0.; 
  
  double weights = 0.; 
  
  double xLo = m_xLo;
  double xHi = m_xHi;
  double yLo = m_yLo;
  double yHi = m_yHi;

  // Only do leakage correction for tracks which "hit" the Calorimeter
  if (cal_top.x()<xLo || cal_top.x()>xHi || cal_top.y()<yLo || cal_top.y()>yHi) 
    return 0;
  
  // Make a unit vector perpendicular to Event Axis (t0)
  // Need to protect against t0 = (0. , 0., -1.) case
  double costheta = t0.z();
  double sintheta = sqrt(1.-costheta*costheta);
  double cosphi  = 1.; 
  if(fabs(sintheta) > .0001) cosphi   = t0.x()/sintheta;
  Vector p(costheta/cosphi, 0., -sintheta);
  p  = p.unit();
  
  // Set the number of inner samples 
  int numInner = (int)(numSamples/2.) ; 
  
  // Loop over samples
  int is = 0;
  for(is = 0; is < numSamples+numInner; is++)
    {
      // Set starting point from this sample trajectory
      Point x0 = cal_top;
      // Note: the inner samples are done at radius/4. and the inner and outer 
      //       samples are rotated by  M_PI/numSamples w.r.t. each other
      if(is <numInner)
        {
          double rotAng = (is-1)*2.*M_PI/numInner; 
          CLHEP::HepRotation rot(t0, rotAng);
          Vector delta = rot*p;
          Point xI = x0 + .25*radius*delta;
          double s = (cal_top.z() - xI.z())/costheta;
          Ray segmt( xI, t0); 
          x0 = segmt.position(s);
        }
      else
        {
          double rotAng = (is-1)*2.*M_PI/numSamples + M_PI/numSamples; 
          CLHEP::HepRotation rot(t0, rotAng);
          Vector delta = rot*p;
          Point xI = x0 + radius*delta;
          double s = (cal_top.z() - xI.z())/costheta;
          Ray segmt( xI, t0); 
          x0 = segmt.position(s);
        } 
      // Check if the start is inside LAT
      if (x0.x()<xLo || x0.x()>xHi || x0.y()<yLo || x0.y()>yHi) continue; 
      
      // Compute the arclength through the CAL
      double s_xp   = (-xHi + x0.x())/t0.x();
      double s_xm   = (-xLo + x0.x())/t0.x();
      double s_minx = (s_xp > s_xm) ? s_xp:s_xm; // Choose soln > 0. 
      
      double s_yp   = (-yHi + x0.y())/t0.y();
      double s_ym   = (-yLo + x0.y())/t0.y();
      double s_miny = (s_yp > s_ym) ? s_yp:s_ym; // Choose soln > 0. 
      
      double s_minz = -(m_calZTop - m_calZBot)/t0.z();
      // Now pick min. soln. of the x, y, and z sides 
      double s_min  = (s_minx < s_miny) ? s_minx:s_miny;  
      s_min         = (s_min  < s_minz) ? s_min :s_minz;
      
      // Set up a propagator to calc. rad. lens. 

      // might as well leave this in, probably will need it soon!
      //std::cout << "EvtValsTool propagator " << x0 << " " << t0 << std::endl;

      m_G4PropTool->setStepStart(x0, t0);
      m_G4PropTool->step(s_min);  
      
      // Loop over the propagator steps to extract the materials
      int numSteps = m_G4PropTool->getNumberSteps();
      double rl_CsI       = 0.;
      double rl_CsICntr   = 0.; 
      double rl_Stuff     = 0.;
      double rl_StuffCntr = 0.;
      idents::VolumeIdentifier volId;
      idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();
      int istep  = 0;
      double last_step_z; 
      bool centroid = true;
      for(; istep < numSteps; ++istep)
        {
          volId = m_G4PropTool->getStepVolumeId(istep);
          volId.prepend(prefix);
          //std::cout << istep << " " << volId.name() << std::endl;
          bool inXtal = ( volId.size()>7 && volId[0]==0 
                          && volId[3]==0 && volId[7]==0 ? true : false );
          double radLen_step = m_G4PropTool->getStepRadLength(istep);
          double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
          Point x_step       = m_G4PropTool->getStepPosition(istep);
          if(istep == 0) last_step_z = x_step.z(); 
          if(inXtal)
            {
              //std::cout << "inXtal " << volId.name() << " " << arcLen_step << " " << radLen_step << std::endl;
              rl_CsI  += radLen_step;
              if(is < numInner) m_arcLen_CsI  += arcLen_step;
            }
          else
            {
              rl_Stuff += radLen_step;
              if(is < numInner) m_arcLen_Stuff += arcLen_step;
            }
          if(x_step.z() >= calcentroidposz || centroid) 
            {
              double step_frac = 1.; 
              if(x_step.z() <= calcentroidposz) 
                {
                  double denominator = last_step_z - x_step.z();
                  
                  centroid = false;
                  
                  // Protect against the case where last_step_z and x_step.z() are the same
                  // (see above where the two are set equal for the first step, it can happen
                  //  that this code is executed on the first step...)
                  if (denominator < 0.001) step_frac = 0.;
                  else                     step_frac = (last_step_z - calcentroidposz) / denominator;
                }
              if(is < numInner) m_arcLen_Cntr += step_frac*arcLen_step;
              if(!inXtal) rl_StuffCntr  += step_frac*radLen_step;
              else        rl_CsICntr    += step_frac*radLen_step;
              last_step_z = x_step.z();
            }
        }
      // Compute sample weighting factor - using Simpson's 1:4:1 Rule
      double sample_factor = (is < numInner) ? 4.:1.; 
      weights += sample_factor;
      
      // Increment accumlation variables
      m_radLen_CsI   += rl_CsI*sample_factor;
      m_rms_RL_CsI   += rl_CsI*rl_CsI*sample_factor;
      m_radLen_Stuff += rl_Stuff*sample_factor;
      m_rms_RL_Stuff += rl_Stuff*rl_Stuff*sample_factor;
      m_radLen_Cntr  += (rl_CsICntr+rl_StuffCntr)*sample_factor;
      m_rms_RL_Cntr  += (rl_CsICntr+rl_StuffCntr)*(rl_CsICntr+rl_StuffCntr)*sample_factor; 
      m_radLen_CntrStuff += rl_StuffCntr*sample_factor;
      m_rms_RL_CntrStuff += rl_StuffCntr*rl_StuffCntr*sample_factor;
    }
  
  // Form the results
  if(weights < 1.) return 0;
  m_radLen_CsI   /= weights;
//   m_rms_RL_CsI    = sqrt(m_rms_RL_CsI/weights -m_radLen_CsI*m_radLen_CsI);
//   m_radLen_Stuff /= weights;
//   m_rms_RL_Stuff  = sqrt(m_rms_RL_Stuff/weights - m_radLen_Stuff*m_radLen_Stuff);
//   m_radLen_Cntr  /= weights;
//   m_rms_RL_Cntr   = sqrt(m_rms_RL_Cntr/weights - m_radLen_Cntr*m_radLen_Cntr);
//   m_radLen_CntrStuff  /= weights;
//   m_rms_RL_CntrStuff   = sqrt(m_rms_RL_CntrStuff/weights - m_radLen_CntrStuff*m_radLen_CntrStuff);
  
//   double innerNo = std::max(1., 1.*numInner);
//   m_arcLen_Stuff /= innerNo;
//   m_arcLen_CsI   /= innerNo;
//   m_arcLen_Cntr  /= innerNo;
  
  return m_radLen_CsI;
}
