/** @file CalValsCorrTool.cxx
@brief implementation of the class CalValsCorrTool

$Header$

*/

#include "CalValsCorrTool.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ToolFactory.h"

#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrUtil/ITkrGeometrySvc.h"


static const ToolFactory<CalValsCorrTool>  s_factory;
const IToolFactory& CalValsCorrToolFactory = s_factory;

#include "TMath.h"
#include <stdexcept>

namespace {


    double beta(double E) {
        // Following expression for b (beta) comes from CalRecon
        // double beta      = exp(0.031*log(CAL_EnergySum/1000.))/2.29;
        // The above expression gives too small values
        double b      = 1.2*exp(0.031*log(E/1000.))/2.3;   //July 16, 2003 * 1.2
        //double b      = .51*exp(0.02*log(E/1000.));  // 17-July
        return b;
    }

    double alpha(double E) {
        // Following expression for a (alpha) comes from CalRecon
        // double alpha      = 2.65*exp(0.15*log(E/1000.));
        // The above expression gives too large values
        //double a = 2.50*exp(0.19*log(E/1000.));    //RR182
        // Found that these power law models diverge as E increases.  
        //double a = .98*log10(E) - .23; 
        double a = .9*log10(E); 
        return a;
    }

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


CalValsCorrTool::CalValsCorrTool( const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent)
                                 : EnergyCorr(type,name,parent)
{
    // declare base interface for all consecutive concrete classes
    declareInterface<IEnergyCorr>(this);
}



StatusCode CalValsCorrTool::initialize()

// This function does following initialization actions:
//    - extracts geometry constants from xml file using GlastDetSvc

{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log << MSG::INFO << "Initializing CalValsCorrTool" <<endreq;

    sc = serviceLocator()->service( "EventDataSvc", m_pEventSvc, true );
    if(sc.isFailure()){
        log << MSG::ERROR << "Could not find EventDataSvc" << std::endl;
        return sc;
    }
    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
        return StatusCode::FAILURE;
    }


    // find TkrGeometrySvc service
    if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
        log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
        return StatusCode::FAILURE;
    }

    m_towerPitch = m_geoSvc->towerPitch();
    m_xNum = m_geoSvc->numXTowers();
    m_yNum = m_geoSvc->numYTowers();

    if (getCalInfo().isFailure()) {
        log << MSG::ERROR << "Couldn't initialize the CAL constants" << endreq;
        return StatusCode::FAILURE;
    }

    // pick up the chosen propagator
    if (service("GlastPropagatorSvc", m_propSvc, true).isFailure()) {
        log << MSG::ERROR << "Couldn't find the GlastPropagatorSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    IToolSvc* toolSvc = 0;
    if(service("ToolSvc", toolSvc, true).isFailure()) {
        log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
        return StatusCode::FAILURE;
    }
    if(!toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
        log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
        return StatusCode::FAILURE;
    }

    return sc;
}


StatusCode CalValsCorrTool::doEnergyCorr(double, Event::CalCluster* cluster)

//Purpose and method:
//
//   This function calls CalValsTool and extracts the crack/leakage corrected
//   energy, storing it in CalCluster.
// 
// TDS input: none
// TDS output: CalClusters


{
    MsgStream lm(msgSvc(), name());
    StatusCode sc =  calculate();
    double correctedEnergy = CAL_Energy_Corr;
    if (sc == StatusCode::SUCCESS) cluster->setEnergyCorrected(correctedEnergy);

    return sc;
}

StatusCode CalValsCorrTool::calculate()
{
    StatusCode sc = StatusCode::SUCCESS;

    //MsgStream logstream(msgSvc(), name());

    // Recover Track associated info. 
    SmartDataPtr<Event::TkrFitTrackCol>  
        pTracks(m_pEventSvc,EventModel::TkrRecon::TkrFitTrackCol);
    SmartDataPtr<Event::TkrVertexCol>     
        pVerts(m_pEventSvc,EventModel::TkrRecon::TkrVertexCol);
    //SmartDataPtr<Event::TkrClusterCol> 
    //    pClusters(m_pEventSvc,EventModel::TkrRecon::TkrClusterCol);
    SmartDataPtr<Event::CalClusterCol>     
        pCals(m_pEventSvc,EventModel::CalRecon::CalClusterCol);
    SmartDataPtr<Event::CalXtalRecCol> 
        pxtalrecs(m_pEventSvc,EventModel::CalRecon::CalXtalRecCol);

    //Do some vital initializations
    CAL_EnergySum     = 0.;
    CAL_Energy_LLCorr = 0.;
    CAL_EneSum_Corr   = 0.;
    CAL_Energy_Corr   = 0.;

    //Make sure we have valid cluster data
    if (!pCals) return sc;

    double z0 = 0.0; // placeholder for offset

    Event::CalCluster* calCluster = pCals->front();

    CAL_EnergySum   = calCluster->getEnergySum();
    for(int i = 0; i<8; i++) CAL_eLayer[i] = calCluster->getEneLayer(i);
    CAL_Long_Rms      = calCluster->getRmsLong();
    CAL_Trans_Rms     = calCluster->getRmsTrans();

    if(CAL_EnergySum>0.0) {
        CAL_Lyr0_Ratio  = CAL_eLayer[0]/CAL_EnergySum;
        CAL_Lyr7_Ratio  = CAL_eLayer[7]/CAL_EnergySum;
        CAL_BkHalf_Ratio = (CAL_eLayer[4]+CAL_eLayer[5]+
            CAL_eLayer[6]+CAL_eLayer[7])/CAL_EnergySum;
        CAL_LRms_Ratio    = CAL_Long_Rms / CAL_EnergySum;
    }

    CAL_Energy_LLCorr = calCluster->getEnergyLeak();

    //Code from meritAlg
    int no_xtals=0;
    int no_xtals_trunc=0;
    double max_log_energy = 0.; 
    for( Event::CalXtalRecCol::const_iterator jlog=pxtalrecs->begin(); jlog != pxtalrecs->end(); ++jlog){

        const Event::CalXtalRecData& recLog = **jlog;

        double eneLog = recLog.getEnergy();
        if(eneLog > max_log_energy) max_log_energy = eneLog; 
        if(eneLog>0)no_xtals++;
        if(eneLog>0.01*CAL_EnergySum)no_xtals_trunc++;
    }
    CAL_Xtal_Ratio= (no_xtals>0) ? float(no_xtals_trunc)/no_xtals : 0;
    CAL_No_Xtals_Trunc = float(no_xtals_trunc); 
    CAL_Xtal_maxEne = max_log_energy; 

    //Overwritten below, if there are tracks
    CAL_Energy_Corr = calCluster->getEnergyCorrected(); 
    if(CAL_EnergySum < 5.) return sc;  

    Point  cal_pos  = calCluster->getPosition();
    Vector cal_dir  = calCluster->getDirection();

    CAL_xEcntr      = cal_pos.x();
    CAL_yEcntr      = cal_pos.y();
    CAL_zEcntr      = cal_pos.z();
    CAL_xdir        = cal_dir.x();
    CAL_ydir        = cal_dir.y();
    CAL_zdir        = cal_dir.z();

    int view = 0;
    Point pos(CAL_xEcntr, CAL_yEcntr, 0);
    CAL_TwrEdgeCntr = activeDist(pos, view);

    CAL_TE_Nrm = ( view==0 ? cal_dir.x() : cal_dir.y() );

    // with no tracks we've done what we can do.
    if(!pTracks) return sc; 
    int num_tracks = pTracks->size(); 
    if(num_tracks <= 0 ) return sc;

    // Get the first track
    Event::TkrFitConPtr pTrack1 = pTracks->begin();
    Event::TkrKalFitTrack* track_1 = dynamic_cast<Event::TkrKalFitTrack*>(*pTrack1);

    // Get the start and direction 
    Point  x0 = track_1->getPosition();
    Vector t0 = track_1->getDirection();

    // If vertexed - use first vertex
    if(pVerts) {
        Event::TkrVertexColPtr gammaPtr =  pVerts->begin(); 
        Event::TkrVertex *gamma = *gammaPtr; 
        x0 = gamma->getPosition();
        t0 = gamma->getDirection();
    }

    // this "cos(theta)" doesn't distinguish between up and down
    double costh  = fabs(t0.z()); 
    // This "phi" is restricted to the range 0 to pi/2
    // protect against zero denominator
    double phi_90    = (fabs(t0.x())<1.e-7) ? 0.5*M_PI : fabs(atan(-t0.y()/t0.x()));

    // Find the distance for energy centroid to track axis
    Vector x_diff = x0 - cal_pos;
    double x_diff_sq = x_diff*x_diff;
    double x_diff_t0 = x_diff*t0;
    CAL_Track_DOCA = sqrt(x_diff_sq - x_diff_t0*x_diff_t0);

    // Compare Tkr and Cal directions.  Note: the direction in Cal is opposite to tracking!
    if(fabs(cal_dir.x()) < 1.) {
        double cosCalt0 = -t0*cal_dir; 
        CAL_Track_Angle = acos(cosCalt0);
    }
    else CAL_Track_Angle = -.1; 

    // Section to apply edge and leakage correction to Cal data. 
    // First apply layer by layer edge correction 
    Vector t_axis = x_diff.unit();  // Using "event" axis. Alternative: t_axis = - t0
    Ray axis(x0, t_axis); 
    double arc_len = (x0.z()- m_calZTop)/t_axis.z(); 
    Point cal_top = axis.position(-arc_len); // Event axis entry point to top of Cal Stack 
    axis      = Ray(cal_top, t_axis); // alternative: Ray(cal_top, -t0); //
    CAL_x0    = cal_top.x();
    CAL_y0    = cal_top.y();
    CAL_z0    = cal_top.z();
    std::vector<Vector> pos_Layer = calCluster->getPosLayer();
    std::vector<double> ene_Layer = calCluster->getEneLayer();

    double ene_sum_corr = 0.;

    // Factors controlling edge correction: core radius, fringe radius, core fraction
    // and size of inter-tower gap (active-to-active area) 
    // What counts is radius of the shower at the CAL entrance.  At low energy, for early
    // conversions there should be a large radius.  
    double rm_hard   = 36. + 36.*cal_trans((CAL_EnergySum-250)/200.)* arc_len/500.;  //was 18 - 36mm is ~Rm
    double rm_soft   = 65. + 65.*cal_trans((CAL_EnergySum-250)/200.)* arc_len/500.;  
    double hard_frac =  .5  + .5*cal_trans((670.- CAL_EnergySum)/300.);   

    // This is the effective inter-tower CalModule gap: approx = towerpitch - calwidth = 44 - 46 mm
    double gap       = 36.;

    // Find the distance from the LAT edge for the leading track
    Event::TkrFitTrackBase::TrackEnd end = Event::TkrFitTrackBase::End;
    Ray trj_1 = Ray(track_1->getPosition(end), track_1->getDirection(end));
    double delta_z = trj_1.position().z() - m_calZTop;
    arc_len = delta_z/fabs(trj_1.direction().z());

    Point cal_1 = trj_1.position(arc_len);

    // Get the lower and upper limits for the CAL in the installed towers
    double deltaX = 0.5*(m_xNum*m_towerPitch - m_calXWidth);
    double deltaY = 0.5*(m_yNum*m_towerPitch - m_calYWidth);
    double calXLo = m_geoSvc->getLATLimit(0,LOW)  + deltaX;
    double calXHi = m_geoSvc->getLATLimit(0,HIGH) - deltaX;
    double calYLo = m_geoSvc->getLATLimit(1,LOW)  + deltaY;
    double calYHi = m_geoSvc->getLATLimit(1,HIGH) - deltaY;

    // Find the distance closest to an edge
    double calX = cal_1.x(); double calY = cal_1.y();
    double dX = std::max(calXLo-calX, calX-calXHi);
    double dY = std::max(calYLo-calY, calY-calYHi);
    CAL_LATEdge = -std::max(dX, dY);
    // Apply circle correction layer by layer in the calorimeter
    // This is a essentially just a geometric correction which is 
    // "integrated" layer by layer through the the Calorimeter. 
    double max_corr = 2.5; //This limits the size of the edge correction (was 1.4)
    double edge_corr = 0.; 
    double good_layers = 0.; 
    for(int i=0; i<8; i++){
        if(ene_Layer[i] < 5.) {
            ene_sum_corr += ene_Layer[i];
            continue; 
        }
        Vector pos = cal_top;
        double arc_len = (pos_Layer[i] - pos).magnitude();
        Point xyz_layer = axis.position(-arc_len);

        double in_frac_soft = containedFraction(xyz_layer, gap, rm_soft, costh, phi_90);

        // Cut off correction upon leaving (through a side)
        if(in_frac_soft < .5) in_frac_soft = .5;

        double in_frac_hard = containedFraction(xyz_layer, gap, rm_hard, costh, phi_90);

        double corr_factor 
            = 1./((1.-hard_frac)*in_frac_soft + hard_frac*in_frac_hard);
        if(corr_factor > max_corr) corr_factor = max_corr;  
        double ene_corr = ene_Layer[i]*corr_factor;
        ene_sum_corr += ene_corr;
        edge_corr    += corr_factor;
        good_layers  += 1.; 
    }
    if(ene_sum_corr < 1.) return sc;
    CAL_EdgeSum_Corr = ene_sum_corr/CAL_EnergySum;
    if (good_layers>0) edge_corr /= good_layers;
    CAL_Edge_Corr = edge_corr; 
    double cal_half_width = 0.5*std::max(m_calXWidth, m_calYWidth);

    //       Leakage Correction  
    // First: get the rad.lens. in the tracker 
    double t_tracker = track_1->getTkrCalRadlen();
    // Patch for error in KalFitTrack: 1/2 of first radiator left out
    int layer = track_1->getLayer();
    t_tracker += 0.5*m_geoSvc->getReconRadLenConv(layer)/costh;
    // Need to fix a problem here.  There can be large fluctuations on single
    // trajectories.  This should be fixed in TkrValsTool probably by averaging
    // over a cylinder as we do in CalValsTool.
    double tkr_radLen_nom = 0.; 
    if(layer < 12) tkr_radLen_nom = (11.5 - layer)*.045 + 4*.195 + .03;
    else           tkr_radLen_nom =  (15.5 - layer)*.195 + .03;
    tkr_radLen_nom /= costh;
    if(t_tracker > tkr_radLen_nom * 1.5) t_tracker = tkr_radLen_nom * 1.5;
    if(t_tracker < tkr_radLen_nom * .5) t_tracker  = tkr_radLen_nom * .5;

    // Find the distance in Cal to nearest edge along shower axis
    //       Need to check sides as well as the back
    // double s_exit = -(axis.position().z() - z0)/t_axis.z();//  Not presently used
    // Point tkr_exit= axis.position(s_exit);    // Exit point from tracker

    Point pos0(CAL_x0, CAL_y0, 0);
    CAL_TwrEdge0 = activeDist(pos0, view);
    // Now get averaged radiation lengths
    // The averaging is set for 6 + 1 samples at a radius of rm_hard/4
    // Note: this method fills internal variables such as m_radLen_CsI & m_radLen_Stuff
    if(aveRadLens(cal_top, -t_axis, rm_hard/4., 6) == StatusCode::FAILURE) return sc; 

    double t_cal_tot = m_radLen_CsI + m_radLen_Stuff; //m_radLen_CsI; //s_min/20.; // 
    ene_sum_corr    *=  t_cal_tot/m_radLen_CsI;       // Correct for unseen stuff
    double t_total   =  t_tracker + t_cal_tot;

    // Energy centroid in radiation lengths for this event.
    double t         = t_tracker + m_radLen_Cntr;

    // 2 Iterations to find energy
    // Note: the TMath incomplete gamma functions is really the fractional inner incomplete
    //       gamma function - i.e. its just what we need to do shower models.  
    //Leon: this is the algorithm that gave the best results for estimating
    //      leakage:
    //        a) In the shower model (wallet card) the ratio of the a & b paremters
    //           is the energy centroid (check it for yourself).  Also b is all but
    //           energy independent (not ours ... problem here -see below ... but the wallet card's)
    //           So given our event centroid and b - > make first guess for a.  But the ratio of 
    //           a/b = <t> only in an infinitely deep cal.  
    //        b) This gives first approx. for energy leakage fraction from model.
    //        c) Now find a & b from parametric models using "corrected energy".
    //      I don't think this is great. The reality is that first estimation correction raises
    //      the energy by ~ 10% at 10 GeV.  The second iteration using the a & b models gives
    //      a total correction of ~ + 33%. 
    //      The powerlaw parameterization of b and a are funny. Particularly b - our b for 
    //      CsI is ~ .46 and raising with E (i.e. too small - like smaller then uranium!  - and 
    //      changing too fast).  These are vistiges from the code in CalRecon. 
    //      Overall it seems to work albeit for heuristic methodology. Note however we may get strange
    //      results when used outside the energy range 50 MeV - 18 GeV.  I hope here that 
    //      the stuff Giebells is doing will come to the rescue - transistion from one to the other??
    //   
    //  double a1 = alpha(ene_sum_corr);
    double b1 =  beta(ene_sum_corr); 
    double a1 = b1*t + 1.; //  + 1 is to compensates for finite cal effect.  This would better 
    // match it  to the final answer.  In reality it has a minor effect since all this 
    // typically affects things as log(E). 
    double in_frac_1 = TMath::Gamma(a1, b1*t_total);
    double a2 = alpha(ene_sum_corr/in_frac_1);
    double b2 =  beta(ene_sum_corr/in_frac_1); 
    double in_frac_2 = TMath::Gamma(a2, b2*t_total);

    ene_sum_corr /= in_frac_2; //Using 2nd order solution

    // Final Delta t factor (?).  Again - present in sumlated data - Needs to be understood 
    // Perhaps this is the event to event compensation for variations in the initiation of showers. 
    CAL_t_Pred      = a2/b2*(TMath::Gamma(a2+1.,b2*t_total)/
        TMath::Gamma(a2,b2*t_total)); 
    CAL_deltaT      = t - CAL_t_Pred;

    // The "final" correction derived empirically from analyizing and flattening the 
    // resultant energy in cos(theta) and log10(E) 
    double logEsum = log10(ene_sum_corr); 
    //double ad_hoc_factor = (1.35-.15*logEsum )/(1.+.9*(.6 - costh)*(.6 - costh))/
    // (1.-.125*logEsum + (.125*logEsum)*costh);  
    //double ad_hoc_factor = (1.23 - .065*logEsum)/(1.+.14*(logEsum-1.)*(costh-.74));  
    double ad_hoc_factor = 1.07; //(1.08 + .02*logEsum);  
    // Store the results away 
    CAL_EneSum_Corr = ene_sum_corr * ad_hoc_factor;
    CAL_Energy_Corr = CAL_EnergySum*edge_corr * ad_hoc_factor/in_frac_2; 
    CAL_TotSum_Corr = CAL_EneSum_Corr/CAL_EnergySum;
    //CAL_Total_Corr  = CAL_Energy_Corr/CAL_EnergySum; 
    CAL_CsI_RLn     = m_radLen_CsI;
    CAL_MIP_Diff    = CAL_EnergySum - 12.07*m_radLen_CsI;
    const double minRadLen = 0.1;
    CAL_MIP_Ratio   = CAL_EnergySum /(12.07*std::max(m_radLen_CsI, minRadLen));
    CAL_Tot_RLn     = t_cal_tot;
    CAL_LAT_RLn     = t_total;
    CAL_Cnt_RLn     = t;
    CAL_DeadTot_Rat = m_radLen_Stuff/std::max(minRadLen, t_total);
    CAL_DeadCnt_Rat = m_radLen_CntrStuff/std::max(minRadLen, t);
    //CAL_a_Parm      = a2;
    //CAL_b_Parm      = b2; 

    CAL_Leak_Corr  = in_frac_2;   
    //    CAL_TwrGap      = m_arcLen_Stuff - m_arcLen_Gap;

    return sc;
}

StatusCode CalValsCorrTool::getCalInfo()
{
    m_calZTop = m_geoSvc->calZTop();
    m_calZBot = m_geoSvc->calZBot();
    m_calXWidth = m_geoSvc->calXWidth();
    m_calYWidth = m_geoSvc->calYWidth();

    return StatusCode::SUCCESS;
}

double CalValsCorrTool::activeDist(Point pos, int &view) const
{
    double edge = 0.;
    double x = pos.x();
    double y = pos.y();
    double x_twr = globalToLocal(x, m_towerPitch, m_xNum);
    double y_twr = globalToLocal(y, m_towerPitch, m_yNum);

    if( fabs(x_twr) > fabs(y_twr) ) {
        edge = m_towerPitch/2. - fabs(x_twr);
        view = 0; 
    }
    else {
        edge = m_towerPitch/2. - fabs(y_twr);
        view = 1;
    }
    return edge;
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
    if (x>m_geoSvc->getLATLimit(0,LOW)+0.5*m_towerPitch
        && x<m_geoSvc->getLATLimit(0,HIGH)-0.5*m_towerPitch) 
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
    if (y>m_geoSvc->getLATLimit(1,LOW)+0.5*m_towerPitch
        && y<m_geoSvc->getLATLimit(1,HIGH)-0.5*m_towerPitch) 
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

StatusCode CalValsCorrTool::aveRadLens(Point x0, Vector t0, double radius, int numSamples)
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


    double xLo = m_geoSvc->getLATLimit(0, LOW);
    double xHi = m_geoSvc->getLATLimit(0, HIGH);
    double yLo = m_geoSvc->getLATLimit(1, LOW);
    double yHi = m_geoSvc->getLATLimit(1, HIGH);

    // Only do leakage correction for tracks which "hit" the Calorimeter
    if (CAL_x0<xLo || CAL_x0>xHi || CAL_y0<yLo || CAL_y0>yHi) return StatusCode::FAILURE;

    double costheta = t0.z();
    double sintheta = sqrt(1.-costheta*costheta);
    double cosphi   = t0.x()/sintheta;
    Vector p(costheta/cosphi, 0., -sintheta);
    p  = p.unit();
    int num_good = 0; 
    int is = 0;
    for(is = 0; is < numSamples+1; is++) {
        Point x0(CAL_x0, CAL_y0, CAL_z0);
        if(is > 0) { // Compute new starting point onto of Cal
            double rotAng = (is-1)*2.*M_PI/numSamples; 
            HepRotation rot(t0, rotAng);
            Vector delta = rot*p;
            Point xI = x0 + radius*delta;
            double s = (CAL_z0 - xI.z())/costheta;
            Ray segmt( xI, t0); 
            x0 = segmt.position(s);
        }
        // Check if the start is inside LAT
        if (x0.x()<xLo || x0.x()>xHi || x0.y()<yLo || x0.y()>yHi) continue; 
        num_good += 1; 

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

        // Now loop over the steps to extract the materials
        int numSteps = m_G4PropTool->getNumberSteps();
        double rl_CsI       = 0.;
        double rl_CsICntr   = 0.; 
        double rl_Stuff     = 0.;
        double rl_StuffCntr = 0.;
        idents::VolumeIdentifier volId;
        idents::VolumeIdentifier prefix = m_detSvc->getIDPrefix();
        int istep  = 0;
        for(; istep < numSteps; ++istep) {
            volId = m_G4PropTool->getStepVolumeId(istep);
            volId.prepend(prefix);
            //std::cout << istep << " " << volId.name() << std::endl;
            bool inXtal = ( volId.size()>7 && volId[0]==0 
                && volId[3]==0 && volId[7]==0 ? true : false );
            double radLen_step = m_G4PropTool->getStepRadLength(istep);
            double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
            Point x_step       = m_G4PropTool->getStepPosition(istep);
            if(inXtal) {
                //std::cout << "inXtal " << volId.name() << " " << arcLen_step << " " << radLen_step << std::endl;
                rl_CsI  += radLen_step;
                if(is == 0) m_arcLen_CsI  += arcLen_step;
            }
            else {
                rl_Stuff += radLen_step;
                if(is == 0) m_arcLen_Stuff += arcLen_step;
            }
            if(x_step.z() >= CAL_zEcntr) {  
                if(is == 0) m_arcLen_Cntr  += arcLen_step;
                if(!inXtal) rl_StuffCntr += radLen_step;
                else        rl_CsICntr   += radLen_step;
            }
        }
        // Increment accumlation variables
        m_radLen_CsI   += rl_CsI;
        m_rms_RL_CsI   += rl_CsI*rl_CsI;
        m_radLen_Stuff += rl_Stuff;
        m_rms_RL_Stuff += rl_Stuff*rl_Stuff;
        m_radLen_Cntr  += rl_CsICntr+rl_StuffCntr;
        m_rms_RL_Cntr  += (rl_CsICntr+rl_StuffCntr)*(rl_CsICntr+rl_StuffCntr); 
        m_radLen_CntrStuff += rl_StuffCntr;
        m_rms_RL_CntrStuff += rl_StuffCntr*rl_StuffCntr;
    }
    m_radLen_CsI   /= num_good;
    m_rms_RL_CsI    = sqrt(m_rms_RL_CsI/num_good -m_radLen_CsI*m_radLen_CsI);
    m_radLen_Stuff /= num_good;
    m_rms_RL_Stuff  = sqrt(m_rms_RL_Stuff/num_good - m_radLen_Stuff*m_radLen_Stuff);
    m_radLen_Cntr  /= num_good;
    m_rms_RL_Cntr   = sqrt(m_rms_RL_Cntr/num_good - m_radLen_Cntr*m_radLen_Cntr);
    m_radLen_CntrStuff  /= num_good;
    m_rms_RL_CntrStuff   = sqrt(m_rms_RL_CntrStuff/num_good - m_radLen_CntrStuff*m_radLen_CntrStuff);
    return StatusCode::SUCCESS;
}

StatusCode CalValsCorrTool::finalize()
{
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

StatusCode CalValsCorrTool::execute()
{
    StatusCode sc = StatusCode::SUCCESS;

    return sc;
}

