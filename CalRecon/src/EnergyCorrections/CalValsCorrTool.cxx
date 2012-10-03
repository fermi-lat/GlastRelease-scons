/** @file CalValsCorrTool.cxx
@brief implementation of the class CalValsCorrTool

$Header$

*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalEventEnergy.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

#include <CalRecon/ICalEnergyCorr.h>
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrUtil/ITkrGeometrySvc.h"
#include "geometry/Ray.h"

#include "CLHEP/Vector/Rotation.h"

#include "TMath.h"
#include <stdexcept>

/**   
* @class CalValsCorrTool
* @author Bill Atwood
*
* Tool to corrected energy for cracks and leakage.
*
* Copied by THB from AnalysisNtuple::CalValsTool.cxx revision 1.43
*
* $Header$
*/

class CalValsCorrTool : public AlgTool, virtual public ICalEnergyCorr
{
public:

    //! destructor
    CalValsCorrTool( const std::string& type, const std::string& name, const IInterface* parent);
    ~CalValsCorrTool() {}; 

    StatusCode initialize();

    // worker function to get the corrected energy      
    Event::CalCorToolResult* doEnergyCorr(Event::CalCluster*, Event::TkrTree* );



private:

    /// Bill's calculation here
    void  calculate(Point x0, Vector t0, double tkr_RLn, double tkr_Energy);
    double containedFraction(Point  pos, double gap, double r, double costh, double phi) const;
    StatusCode aveRadLens(Point x0, Vector t0, double radius, int numSamples);
        Event::CalCorToolResult* loadResults(); 

    /// gets the CAL info from detModel
    StatusCode getCalInfo();

    /// TkrGeometrySvc used for access to tracker geometry info
    ITkrGeometrySvc*   m_tkrGeom;

    /// Pointer to the Gaudi data provider service
    IDataProviderSvc*  m_dataSvc;

    /// G4 Propagator tool
    IPropagator *      m_G4PropTool; 

    /// Detector Service
    IGlastDetSvc *     m_detSvc; 

    /// Internal Cluster and vertex data
    Event::CalCluster* m_cluster;
    Event::TkrTree*    m_tree; 

    /// some Geometry
    double             m_towerPitch;
    int                m_xNum;
    int                m_yNum;
    
    /// CAL vars
    double             m_calXWidth;
    double             m_calYWidth;
    double             m_calZTop;
    double             m_calZBot;
    Point              m_cal_pos;
    Point              m_cal_top; 
    Vector             m_cal_dir;

    // Internal Variables
    double             m_radLen_CsI, m_rms_RL_CsI; //Rad. lengths along axis in CsI Xtals (and rms)
    double             m_radLen_Stuff, m_rms_RL_Stuff; // Rad. lengths of other material along axis (and rms)
    double             m_radLen_Cntr, m_rms_RL_Cntr; // Rad. lengths along axis to centroid (and rms)
    double             m_radLen_CntrStuff, m_rms_RL_CntrStuff; // Rad. length of other material (not CsI) to centroid

    double             m_arcLen_CsI;      // Path length along shower axis in CsI Xtals
    double             m_arcLen_Stuff;    // Path length along shower axis for material other then CsI
    double             m_arcLen_Cntr;     // Path length along shower axis to energy centroid

    double             m_gap_fraction;    // Fraction of edge cylinder that falls inside a gap
    double             m_edge_correction; // Edge-Gap multiplicative corrections factor
    double             m_leakage_correction; // Containment fraction for shower
    double             m_total_correction;// Total multiplicative correction - includes ad-hoc piece - see code
    double             m_raw_energy;      // Raw summed energy from Xtals
    double             m_corr_energy;     // Fully corrected energy 
    double             m_deltaT;          // Difference between predicted energy centriod and meas. centroid
    double             m_t_Pred;          // Predicted energy centriod
    double             m_t;               // Location of measured energy centroid in rad. len.
    double             m_t_total;         // Total rad. len. along event axis used in leakage calc. 
    unsigned int       m_status_bits; // Status bits to be set in CalCorResults

   /// Control Parameters set via JobOptions parameters
    double             m_minEnergy;      // Min. energy required to due the corrections
    double             m_maxEdgeCorr;    // Max. allowed edge corretion factor
    double             m_edgeFracCutOff; // Contained fraction (edges) min. (or cutoff)
    double             m_minCorrEnergy;  // Min Energy for which to make leakage correction
    double             m_leakConvergence;// Leakage convergence fraction 
    double             m_minLeakFrac;    // Min allowed leakage fraction 
    double             m_minCsIRLn;      // Min. rad. len. of CsI for leakage correction
};

#include "GaudiKernel/DeclareFactoryEntries.h"
DECLARE_TOOL_FACTORY(CalValsCorrTool) ;

namespace {
  /// Shower model parameters a & b (alpha & beta)
   double beta(double E) {
                double b       = .44+.029*log(E/1000.)/2.3026; //Wallet card fit post 100GeV run 29-apr-05
        return b;
    }  
    double alpha(double E) {
                double a       = 2.496+1.191*log(E/1000.)/2.3026; //Post 100 GeV run 18-apr-05
        return a;
    }
 
    /// A local erf function: alg. from Receipts in C
    double erf_cal(double x) {
        double t = 1./(1.+.47047*x);
        double results = 1. -(.34802 - (.09587 -.74785*t)*t)*t*exp(-x*x);
        return results;
    }
    double cal_trans(double x) {
        if(x < 0) return (.5*(1. + erf_cal(-x)));
        else      return (.5*(1. - erf_cal( x)));
    }

    /// sign of a number
    double sign(double x) { return x>0 ? 1.: -1. ;}


    /// turn a global coordinate (tower, ladder, wafer) roughly into a local one
    double globalToLocal(double x, double pitch, int n) {
        double xNorm = x/pitch + 0.5*n;
        return sign(x)*(fmod(fabs(xNorm),1.0) - 0.5)*pitch ;
    }

    double circleFraction(double r) {
        double rl = (fabs(r) < 1.) ? fabs(r):1.; 
        double a_slice = 2.*(M_PI/4. - rl*sqrt(1.-rl*rl)/2. - asin(rl)/2.);
        double in_frac = 1.-a_slice/M_PI;
        if(r < 0.) in_frac = a_slice/M_PI;
        return in_frac;
    }

    double circleFractionSimpson(double r, double angle_factor) {
        double slice_0 = circleFraction(r);
        double slice_p = circleFraction(r+angle_factor);
        double slice_m = circleFraction(r-angle_factor);
        return (slice_p + 4.*slice_0 + slice_m)/6.;
    }
}

CalValsCorrTool::CalValsCorrTool( const std::string & type, 
                                  const std::string & name, 
                                  const IInterface * parent )
                                : AlgTool(type,name,parent)
{ 
    declareInterface<ICalEnergyCorr>(this) ;

        //Declare the control parameters. Defaults appear here
    declareProperty("MinEnergy",          m_minEnergy      = 10.);
        declareProperty("MaxEdgeCorr",        m_maxEdgeCorr    = 2.5);
        declareProperty("EdgeFractionCutOff", m_edgeFracCutOff = .5);
        declareProperty("MinLeakCorrEnergy",  m_minCorrEnergy  = 50.);
        declareProperty("LeakConvergenceFrac",m_leakConvergence= .005); 
        declareProperty("MinLeakFraction",    m_minLeakFrac    = .10); 
        declareProperty("MinCsIRadLens",      m_minCsIRLn      = 2.0); 

}

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc
StatusCode CalValsCorrTool::initialize()
{
//    if (EnergyCorr::initialize().isFailure())
//    { 
//        return StatusCode::FAILURE ; 
//    }

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing CalValsCorrTool" <<endreq;

    // find TkrGeometrySvc service
    if ((sc = service("TkrGeometrySvc", m_tkrGeom, true)).isFailure())
    {
        throw GaudiException("Service [TkrGeometrySvc] not found", name(), sc);
    }

    m_towerPitch = m_tkrGeom->towerPitch();
    m_xNum       = m_tkrGeom->numXTowers();
    m_yNum       = m_tkrGeom->numYTowers();

    if ((sc = getCalInfo()).isFailure()) 
    {
        throw GaudiException("Could not initialize the CAL constants", name(), sc);
    }

    if ((sc = service("GlastDetSvc", m_detSvc, true)).isFailure())
    { 
        throw GaudiException("Service [GlastDetSvc] not found", name(), sc);
    }

    //Locate and store a pointer to the data service which allows access to the TDS
    if ((sc = service("EventDataSvc", m_dataSvc)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
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

    return sc;
}


Event::CalCorToolResult* CalValsCorrTool::doEnergyCorr(Event::CalCluster* cluster, Event::TkrTree* tree)
{
    //Purpose and method:
    //
    //   This function calls CalValsTool and extracts the crack/leakage corrected
    //   energy, storing it in CalCluster.
    // 
    // TDS input: none
    // TDS output: CalClusters

    Event::CalCorToolResult* corResult = 0;
    MsgStream lm(msgSvc(), name());

    if (!cluster)
    {
        lm << MSG::DEBUG << "Ending doEnergyCorr: No Cluster" 
            << endreq;
        return corResult;
    }
    
    m_cluster = cluster;
    m_tree    = tree; 

        //Make sure we have valid cluster data
    if (!m_cluster) return corResult;

    // Put here a place holder for Event Axis Calculation!!!!!!!!!!!!!
    if(!m_cluster->checkStatusBit(Event::CalCluster::CENTROID)) return corResult;

    m_cal_pos = m_cluster->getPosition();
    Point x0  = m_cal_pos;
    Vector t0 = m_cluster->getDirection();

    if( !m_cluster->checkStatusBit(Event::CalCluster::MOMENTS) ||
              m_cluster->getRmsLong() < .1) { // Trap NaN condition caused by Moments failure
        t0 = Vector(0., 0., 1.);
    }

    double tkr_Energy = 0.; 
    double tkr_RLn    = 0.;

    if (m_tree != 0) 
    {
        x0 = m_tree->getAxisParams()->getEventPosition();
        t0 = m_tree->getAxisParams()->getEventAxis();    

        // Start at the top most point (ie in the top silicon layer) 
        int    topLayer = m_tree->getHeadNode()->front()->getTreeStartLayer();
        double zTopLyr  = std::max(m_tkrGeom->getLayerZ(topLayer, 0), m_tkrGeom->getLayerZ(topLayer, 1));

        // Translate x0 to this z position
        double arcLen   = t0.z() > 0. ? (zTopLyr - x0.z()) / t0.z() : 0.;

        // Translate the tree start position to middle of top silicon layer
        x0 = x0 + arcLen * t0;

        // Now set up and call propagator to get the radiation lengths to calorimeter
        arcLen = t0.z() > 0. ? (x0.z() - m_tkrGeom->calZTop()) / t0.z() : 0.;
        m_G4PropTool->setStepStart(x0, -t0);
        m_G4PropTool->step(arcLen);
        tkr_RLn = m_G4PropTool->getRadLength(); 

        // Patch for error in KalFitTrack: 1/2 of first radiator left out
        int topPlane = m_tkrGeom->getPlane(zTopLyr);

        if (m_tkrGeom->isTopPlaneInLayer(topPlane)) 
        {
            tkr_RLn += 0.5*m_tkrGeom->getRadLenConv(topLayer) / t0.z();
        }

        // add up the rad lens; this could be a local array if you're bothered by the overhead
        //   but hey, compared to the propagator...
        double tkr_radLen_nom = 0.; 
        int layerCount = topLayer;
        for(; layerCount>=0; --layerCount) 
        {
            tkr_radLen_nom += m_tkrGeom->getRadLenConv(layerCount) 
                             + m_tkrGeom->getRadLenRest(layerCount);
        }
        tkr_radLen_nom /= t0.z();
        if      (tkr_RLn > tkr_radLen_nom * 1.5) {tkr_RLn = tkr_radLen_nom * 1.5;}
        else if (tkr_RLn < tkr_radLen_nom * .5)  {tkr_RLn  = tkr_radLen_nom * .5;}

       // First Tree Based Tracker energy
                double tkrThinClusters = m_tree->getHeadNode()->getNumThinNodesInTree();
                double tkrThickClusters = m_tree->getHeadNode()->getNumThickNodesInTree();
                double tkrBlankClusters = m_tree->getHeadNode()->getNumBlankNodesInTree();
                double tkrThickLeaves = m_tree->getHeadNode()->getNumThickLeavesInTree();
                double tkrThinClsCosTheta = tkrThinClusters/fabs(t0.z());
                double tkrThickClsCosTheta = tkrThickClusters/fabs(t0.z());
                double tkrBlankClsCosTheta = tkrBlankClusters/fabs(t0.z());
                double tkrThinClsCosThetaSq = tkrThinClsCosTheta/fabs(t0.z());
                double tkrThickClsCosThetaSq = tkrThickClsCosTheta/fabs(t0.z()); 

                // Coeffecients from linear regression:  2 fits  - Thin and Thick conversions
                if(topLayer > 5) {
                        tkr_Energy = -0.46*tkrThickLeaves + 1.09*tkrThinClsCosTheta + 4.37*tkrThickClsCosTheta
                                         -1.25*tkrBlankClsCosTheta + .53*tkrThinClsCosThetaSq -.44*tkrThickClsCosThetaSq;
                }
                else {
                        tkr_Energy = -0.97*tkrThickLeaves + 4.39*tkrThickClsCosTheta +
                                          1.37*tkrBlankClsCosTheta  -.28*tkrThickClsCosThetaSq;
                }
                // Layer normalizations - first 2 layer shouldn't happen
                double layerCoefs [18] = { 1., 1., 1.022,.981,1.010,1.018,   
                                               .896, .894, .888, .904, .925, .917, .921, .940, .934, .983, .973, .967};
                tkr_Energy /= layerCoefs[topLayer];
    }

        // Now do the energy correction and calculation of several vars. used in bkg. rejection
    calculate(x0, t0, tkr_RLn, tkr_Energy);

    if (m_status_bits != Event::CalCorToolResult::ZERO) corResult = loadResults();

    return corResult;
}

void CalValsCorrTool::calculate(Point x0, Vector t0, double t_tracker, double tkr_energy)
{
        // Temporary Location for CalValsCorrTools status Bits
         enum statusBits {VALIDEDGECORR  = 0x01000000,  // Edge Correction completed
                          VALIDAVERLN    = 0x02000000,  // Ave. Rad. lengths Calcc. o.k.
                                          MINCSIRLN      = 0x04000000,  // >Min. CsI RLn present to do leakage
                              MAXLEAKCORR    = 0x08000000,  // Leakage exceeds max. allowed correction
                                          MAXLEAKITER    = 0x10000000   // Leakage Corr. failed to converge on iterations
         };

    //Initializations - need to insure a sensible CalCorToolResult
        m_edge_correction    = 1.;
        m_leakage_correction = 1.;
    m_radLen_CsI         = 0.;
    m_radLen_Stuff       = 0.;
    m_radLen_Cntr        = 0.;
    m_radLen_CntrStuff   = 0.;

    m_arcLen_CsI         = 0.;
    m_arcLen_Stuff       = 0.;
    m_arcLen_Cntr        = 0.;

        m_gap_fraction       = 0.;
        m_total_correction   = 1.;
        m_deltaT             = 0.;
        m_t_Pred             = 0.;
        m_t                  = 0.;
        m_t_total            = 0.;
        m_status_bits        = 0;
    m_raw_energy   = m_cluster->getMomParams().getEnergy();
    m_corr_energy = m_raw_energy;
        m_cal_pos  = m_cluster->getPosition();
    m_cal_dir  = m_cluster->getDirection();

// DC: redundant with correctionName
//    m_status_bits |=  Event::CalCorToolResult::CALVALS;

    // Construct Event Axis along which the shower will be evaluated 
    Ray axis(x0, t0); 
    double arc_len = (x0.z()- m_calZTop)/t0.z(); 
    m_cal_top = axis.position(-arc_len);   // Event axis entry point to top of Cal Stack 
    axis      = Ray(m_cal_top, t0); 

    // Note: this method fills internal variables such as m_radLen_CsI & m_radLen_Stuff
    double rm_hard   = 40. + 36.*cal_trans((m_raw_energy-250)/200.)* arc_len/500.;

    if(aveRadLens(m_cal_top, -t0, rm_hard/4., 6) == StatusCode::FAILURE) return; 
        m_status_bits |= VALIDAVERLN; 

    double t_cal_tot  = m_radLen_CsI + m_radLen_Stuff;// rad. len. in Cal
    m_t_total         = t_tracker + t_cal_tot;        // Total rad. len. in LAT
    m_t               = t_tracker + m_radLen_Cntr;    // Energy centroid in rad. len.

    if(m_raw_energy < m_minEnergy) return;  
        m_status_bits |=  Event::CalCorToolResult::VALIDPARAMS;

    // this "cos(theta)" doesn't distinguish between up and down
    double costh  = fabs(t0.z()); 
    // This "phi" is restricted to the range 0 to pi/2
    // protect against zero denominator
    double phi_90    = (fabs(t0.x())<1.e-7) ? 0.5*M_PI : fabs(atan(-t0.y()/t0.x()));

        //        Edge Correction Calculation - done layer-by-layer
    // Factors controlling edge correction: core radius, fringe radius, core fraction
    // and size of inter-tower gap (active-to-active area) 
    // What counts is radius of the shower at the CAL entrance.  At low energy, for early
    // conversions there should be a large radius.  
         double rm_soft   = 80. + 65.*cal_trans((m_raw_energy-250)/200.)* arc_len/500.;
        double hard_frac =  .5  + .5*cal_trans((670.- m_raw_energy)/300.);   
 
    double test_energy = 0;
    double arg1 = -1.25;
    double arg2 = 670./300.;
    //std::cout << "vals = " << cal_trans(arg1) << ", " << cal_trans(arg2) << std::endl;
    double rm_h1 = 40. + 36.*cal_trans(arg1)*250/500.;
    double rm_hf1 = .5 + .5*cal_trans(arg2);
    //std::cout << "rm_h1/hf1 = " << rm_h1 << ", " << rm_hf1 << std::endl;  
     
    // This is the effective inter-tower CalModule gap: approx = towerpitch - calwidth = 44 - 46 mm
        double gap       = 45.; //Physical gap

        // Gap loss factor - for the fraction of the shower in this layer which is in an inter-tower gap
        // this is the fraction that is just lost and not passed on to the next layer. For the last CAL
        // Layer - this factor is set to 1.0
    double gap_loss_factor = .30 + .45*costh;

    // Get the lower and upper limits for the CAL in the installed towers
//    double deltaX = 0.5*(m_xNum*m_towerPitch - m_calXWidth);
//    double deltaY = 0.5*(m_yNum*m_towerPitch - m_calYWidth);
//    double calXLo = m_tkrGeom->getLATLimit(0,LOW)  + deltaX;
//    double calXHi = m_tkrGeom->getLATLimit(0,HIGH) - deltaX;
//    double calYLo = m_tkrGeom->getLATLimit(1,LOW)  + deltaY;
//    double calYHi = m_tkrGeom->getLATLimit(1,HIGH) - deltaY;

    // Apply circle correction layer by layer in the calorimeter
    // This is a essentially just a geometric correction which is 
    // "integrated" layer by layer through the the Calorimeter. 
    // Initialize the various sums (computing <t> and ene_sum with edge corrections)
    Event::CalClusterLayerDataVec& lyrDataVec = *m_cluster;
        m_corr_energy = 0.; 
        m_gap_fraction = 0.; 
//    double edge_corr = 0.; 
    double good_layers = 0.;
        double layer_energy_sum = 0.; 
    for(int i=0; i<8; i++){
           layer_energy_sum += lyrDataVec[i].getEnergy();
       if(lyrDataVec[i].getEnergy() < m_minEnergy/2.) {
            m_corr_energy += lyrDataVec[i].getEnergy();
            continue; 
        }
        double arc_len = (lyrDataVec[i].getPosition().z() - m_cal_top.z())/axis.direction().z();
        Point xyz_layer = axis.position(arc_len);
        
        double in_frac_soft = containedFraction(xyz_layer, gap, rm_soft, costh, phi_90);
        
        // Cut off correction upon leaving (through a side)
        if(in_frac_soft < m_edgeFracCutOff) in_frac_soft = m_edgeFracCutOff;
        
        double in_frac_hard = containedFraction(xyz_layer, gap, rm_hard, costh, phi_90);
        
                if(i == 7) gap_loss_factor = 1.; //NOTE: HARDWIRED IN LAST CAL LAYER 7
        double corr_factor 
            = 1./((1.-hard_frac)*(1.-gap_loss_factor*(1.-in_frac_soft)) + 
                                   hard_frac*(1.-gap_loss_factor*(1.-in_frac_hard)));

        if(corr_factor > m_maxEdgeCorr) corr_factor = m_maxEdgeCorr;  
        double ene_corr = lyrDataVec[i].getEnergy()*corr_factor;
        m_corr_energy += ene_corr;
        m_gap_fraction += (1.-hard_frac)*(1.-in_frac_soft) + 
                                           hard_frac*(1.-in_frac_hard);
        good_layers  += 1.; 
    }
    if(m_corr_energy < m_minEnergy) return;
        m_status_bits |= VALIDEDGECORR; 

        double zero_suppression_energy = m_raw_energy - layer_energy_sum;
        m_corr_energy += zero_suppression_energy; 
    m_edge_correction = m_corr_energy/m_raw_energy;
    if (good_layers>0) m_gap_fraction /= good_layers;

    //           Average Radiation Lengths
    // The averaging is set for 6 + 3 samples at a radius of rm_hard/4

    //// Note: this method fills internal variables such as m_radLen_CsI & m_radLen_Stuff
    //if(aveRadLens(m_cal_top, -t0, rm_hard/4., 6) == StatusCode::FAILURE) return; 
        //m_status_bits |= VALIDAVERLN; 

        if(m_radLen_CsI < m_minCsIRLn) return;
        m_status_bits |= MINCSIRLN; 

    m_corr_energy    *= t_cal_tot/m_radLen_CsI;       // Correction for non-CsI shower

        //                    Energy Leakage Correction
    // Note: the TMath incomplete gamma function is really the fractional inner incomplete
    //       gamma function - i.e. its just what we need to do shower models. 
    // 
        // Algorithm:
        //       The givens are the total rad. len., the location of the observed energy centroid
        //       (in rad. len.), and the (now edge corrected) observed energy.  The true energy 
        //       set the overall shape of the shower model albeit the shape varys slowly
        //       with energy (like log(E)). Given the shower shape, the leakage can be estimated
        //       and the observed energy corrected.  In the iterative proceed below we are effectively
        //       fitting the location of the energy centroid, given the constraints of total rad. len.
        //       and observed energy.  Specifically:
        //          1) the present estimate of the corrected energy givens the shower model b parameter
        //             (the b parameter simply is a scale factor for the rad. len.)
        //          2) the a parameter (location of centroid * b) can then be calculated with a correction
        //             for the finite length (Note: to boot-strap - no correction is made for the zeroth
        //             iteration)
        //          3) with a & b and the total rad. len. the contaiment fraction is estimated and a 
        //             new corrected energy estimated
        //          4) loop back to 1) until the corrected energy changes by less then the convergence 
        //             criteria
        //
        //   Note that in principle we need the FULL observed energy and the FULL rad. len. to do this
        //   This requires the Tracker pieces.  In practice the Tracker rad. len. are vital, however the 
        //   energy from the Tracker only starts to become important when it constitues a large fraction of 
        //   the total energy (e.g. at 100 MeV the Tracker energy is ~ 50% while at a GeV its ~ 10%).  
        //       
    double b1 =  beta(m_corr_energy);
        double a1 = m_t*b1; 
    m_leakage_correction   = TMath::Gamma(a1, b1*m_t_total);
        double e_cal_corr      = m_corr_energy + tkr_energy;
        double e_cal_corr_next = (m_corr_energy + tkr_energy)/m_leakage_correction;
        int counter = 0;
        if(m_raw_energy > m_minCorrEnergy) { //There are convergence issues for small energies
            while ((fabs((e_cal_corr_next-e_cal_corr)/e_cal_corr) > m_leakConvergence) && counter < 20) {
                    e_cal_corr = e_cal_corr_next;
            b1 =  beta(e_cal_corr_next); 
                    a1 = m_t*b1*m_leakage_correction/TMath::Gamma(a1+1.,b1*m_t_total); 
            m_leakage_correction = TMath::Gamma(a1, b1*m_t_total);
                e_cal_corr_next = (m_corr_energy + tkr_energy)/m_leakage_correction;
                        counter++;
            }
        }
        // Limit the leakage correction factor
        if(m_leakage_correction < m_minLeakFrac) {
                m_status_bits |= MAXLEAKCORR; 
                m_leakage_correction = m_minLeakFrac;
        }
        if(counter > 20) m_status_bits |= MAXLEAKITER; 
   
        // Lump the entire correction into the CAL piece (that's why the addition/subtraction of tkr_energy)
    m_corr_energy = (m_corr_energy + tkr_energy)/m_leakage_correction - tkr_energy;
 
    m_t_Pred      = a1/b1*(TMath::Gamma(a1+1.,b1*m_t_total)/
                             TMath::Gamma(a1,b1*m_t_total)); 
    m_deltaT      = m_t - m_t_Pred;


    // The "final" correction derived empirically from analyizing and flattening the 
    // resultant energy in cos(theta) and log10(E)   
    double ad_hoc_factor = (1.0+.20*(1.-costh));

    // Apply final correction 
    m_corr_energy = m_corr_energy * ad_hoc_factor;
    m_total_correction = m_corr_energy/m_raw_energy;

        // NOTE:  BIG Change:  leaving out the leakage correction.  It cause too much dispersion below 1 GeV
        m_corr_energy *= m_leakage_correction; 
    
        return;
}

Event::CalCorToolResult* CalValsCorrTool::loadResults()
{
        // Create a Results Object
        Event::CalCorToolResult *corResult = new Event::CalCorToolResult();

    // Fill in the corrected information and exit
    Event::CalParams params(m_corr_energy, .1*m_corr_energy, 
                            m_cal_pos.x(), m_cal_pos.y(), m_cal_pos.z(), 1., 0., 0., 1., 0., 1.,
                            m_cal_dir.x(), m_cal_dir.y(), m_cal_dir.z(), 1., 0., 0., 1., 0., 1.);

    corResult->setStatusBit(Event::CalCorToolResult::VALIDPARAMS);
    corResult->setCorrectionName(type());
    corResult->setParams(params);
    corResult->setChiSquare(1.);
    (*corResult)["CorrectedEnergy"] = m_corr_energy ;
    (*corResult)["CalTopX0"]        = m_cal_top.x() ;
        (*corResult)["CalTopY0"]        = m_cal_top.y() ;
        (*corResult)["CsIRLn"]          = m_radLen_CsI ;
        (*corResult)["CALRLn"]          = m_radLen_CsI + m_radLen_Stuff ;
    (*corResult)["LATRLn"]          = m_t_total ;
    (*corResult)["StuffRLn"]        = m_radLen_Stuff ;
    (*corResult)["CntrRLn"]         = m_t ;
        (*corResult)["CntrRLnStuff"]    = m_radLen_CntrStuff ;
        (*corResult)["CsIArcLen"]       = m_arcLen_CsI ;
        (*corResult)["GapFraction"]     = m_gap_fraction ;
    (*corResult)["EdgeCorrection"]  = m_edge_correction ;
    (*corResult)["LeakCorrection"]  = m_leakage_correction ;
        (*corResult)["TotalCorrection"] = m_total_correction ;
        (*corResult)["PredCntr"]        = m_t_Pred ;        
    (*corResult)["DeltaCntr"]       = m_deltaT ;        
 
    return corResult;
}

StatusCode CalValsCorrTool::getCalInfo()
{
    m_calZTop = m_tkrGeom->calZTop();
    m_calZBot = m_tkrGeom->calZBot();
    m_calXWidth = m_tkrGeom->calXWidth();
    m_calYWidth = m_tkrGeom->calYWidth();

    return StatusCode::SUCCESS;
}

double CalValsCorrTool::containedFraction(Point pos, double gap, 
                                          double r, double costh, double phi) const
{
    double x = pos.x();
    double y = pos.y();
    // Get the projected angles for the gap
    double tanth = sqrt(1.-costh*costh)/costh;
    double gap_x = gap - 40.*sin(phi)*tanth; 
    if(gap_x < 5.) gap_x = 5.;
    double gap_y = gap - 40.*cos(phi)*tanth;
    if(gap_y < 5.) gap_y = 5.; 

    // X Edges
    double x_twr = globalToLocal(x, m_towerPitch, m_xNum); // member function of ValBase
    double edge = m_towerPitch/2. - fabs(x_twr);
    double r_frac_plus = (edge-gap_x/2.)/r; 
    double angle_factor = sin(phi)*(1./costh - 1.);
    double in_frac_x  =  circleFractionSimpson(r_frac_plus, angle_factor);
    // This should work even for missing towers, 
    // because the radlens, etc., come from the propagator
    if (x>m_tkrGeom->getLATLimit(0,LOW)+0.5*m_towerPitch
        && x<m_tkrGeom->getLATLimit(0,HIGH)-0.5*m_towerPitch) 
    {

        double r_frac_minus = (edge + gap_x/2.)/r;
        in_frac_x += circleFractionSimpson(-r_frac_minus, angle_factor);
    }

    // Y Edges
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);
    edge = m_towerPitch/2. - fabs(y_twr);
    r_frac_plus = (edge-gap_y/2.)/r; 
    angle_factor = cos(phi)*(1./costh - 1.);
    double in_frac_y  =  circleFractionSimpson(r_frac_plus, angle_factor);
    if (y>m_tkrGeom->getLATLimit(1,LOW)+0.5*m_towerPitch
        && y<m_tkrGeom->getLATLimit(1,HIGH)-0.5*m_towerPitch) 
    {// X edge is not outside limit of LAT
        double r_frac_minus = (edge + gap_y/2.)/r;
        in_frac_y += circleFractionSimpson(-r_frac_minus, angle_factor);
    }

    // Cross term assumes x and y are independent 
    double in_frac = 1.;
    if(in_frac_x > .999) in_frac = in_frac_y;
    else if(in_frac_y > .999) in_frac = in_frac_x; 
    else in_frac = 1. - (1.-in_frac_x) - (1.-in_frac_y) + 
        (1.-in_frac_x)*(1.-in_frac_y);  //Cross Term Correction?  
    if(in_frac < .01) in_frac = .01;

    return in_frac;
}

StatusCode CalValsCorrTool::aveRadLens(Point /* x0 */, Vector t0, double radius, int numSamples)
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
        
    double xLo = m_tkrGeom->getLATLimit(0, LOW);
    double xHi = m_tkrGeom->getLATLimit(0, HIGH);
    double yLo = m_tkrGeom->getLATLimit(1, LOW);
    double yHi = m_tkrGeom->getLATLimit(1, HIGH);
    // Ph. Bruel : hardcoded modification in order to take into account the CU geometry (tower 1 without tracker)
    if(xLo==0) xLo = -m_towerPitch;

    // Only do leakage correction for tracks which "hit" the Calorimeter
        if (m_cal_top.x()<xLo || m_cal_top.x()>xHi || m_cal_top.y()<yLo || m_cal_top.y()>yHi) 
                return StatusCode::FAILURE;
    
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
        for(is = 0; is < numSamples+numInner; is++) {

                // Set starting point from this sample trajectory
                Point x0 = m_cal_top;
                // Note: the inner samples are done at radius/4. and the inner and outer 
                //       samples are rotated by  M_PI/numSamples w.r.t. each other
                if(is <numInner) {
                        double rotAng = (is-1)*2.*M_PI/numInner; 
            CLHEP::HepRotation rot(t0, rotAng);
                        Vector delta = rot*p;
                        Point xI = x0 + .25*radius*delta;
                        double s = (m_cal_top.z() - xI.z())/costheta;
                        Ray segmt( xI, t0); 
                        x0 = segmt.position(s);
                }
                else {
                        double rotAng = (is-1)*2.*M_PI/numSamples + M_PI/numSamples; 
                        CLHEP::HepRotation rot(t0, rotAng);
                        Vector delta = rot*p;
                        Point xI = x0 + radius*delta;
                        double s = (m_cal_top.z() - xI.z())/costheta;
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
       for(; istep < numSteps; ++istep) {
           volId = m_G4PropTool->getStepVolumeId(istep);
           volId.prepend(prefix);
           //std::cout << istep << " " << volId.name() << std::endl;
           bool inXtal = ( volId.size()>7 && volId[0]==0 
                          && volId[3]==0 && volId[7]==0 ? true : false );
           double radLen_step = m_G4PropTool->getStepRadLength(istep);
           double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
           Point x_step       = m_G4PropTool->getStepPosition(istep);
                   if(istep == 0) last_step_z = x_step.z(); 
           if(inXtal) {
            //std::cout << "inXtal " << volId.name() << " " << arcLen_step << " " << radLen_step << std::endl;
               rl_CsI  += radLen_step;
               if(is < numInner) m_arcLen_CsI  += arcLen_step;
           }
           else {
               rl_Stuff += radLen_step;
               if(is < numInner) m_arcLen_Stuff += arcLen_step;
           }
           if(x_step.z() >= m_cal_pos.z() || centroid) {
               double step_frac = 1.; 
                           if(x_step.z() <= m_cal_pos.z()) {
                   double denominator = last_step_z - x_step.z();

                   centroid = false;

                   // Protect against the case where last_step_z and x_step.z() are the same
                   // (see above where the two are set equal for the first step, it can happen
                   //  that this code is executed on the first step...)
                   if (denominator < 0.001) step_frac = 0.;
                   else                     step_frac = (last_step_z - m_cal_pos.z()) / denominator;
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
        if(weights < 1.) return StatusCode::FAILURE; 
        m_radLen_CsI   /= weights;
    m_rms_RL_CsI    = sqrt(m_rms_RL_CsI/weights -m_radLen_CsI*m_radLen_CsI);
    m_radLen_Stuff /= weights;
        m_rms_RL_Stuff  = sqrt(m_rms_RL_Stuff/weights - m_radLen_Stuff*m_radLen_Stuff);
        m_radLen_Cntr  /= weights;
        m_rms_RL_Cntr   = sqrt(m_rms_RL_Cntr/weights - m_radLen_Cntr*m_radLen_Cntr);
        m_radLen_CntrStuff  /= weights;
        m_rms_RL_CntrStuff   = sqrt(m_rms_RL_CntrStuff/weights - m_radLen_CntrStuff*m_radLen_CntrStuff);

        double innerNo = std::max(1., 1.*numInner);
        m_arcLen_Stuff /= innerNo;
    m_arcLen_CsI   /= innerNo;
        m_arcLen_Cntr  /= innerNo;

        return StatusCode::SUCCESS;
}
