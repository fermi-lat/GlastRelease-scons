//------------------------------------------------------------------------------
//
//     TkrKalFitTrack
//
//     Implementation of the Tracker Fit Track output class
//				
//------------------------------------------------------------------------------


#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"

using namespace Event;

// Constructor takes no arguments and basically just initializes to
// a null state. In the class' current incarnation, it is expected to
// be inherited from a class which actually does the track fit. 
TkrKalFitTrack::TkrKalFitTrack()
{
    m_status           = EMPTY;
    m_type             = 0;

    m_energy0          = 0.;
    m_x0               = Point(0.,0.,0.);
    m_dir              = Vector(0.,0.,0.);

    m_chisq            = -1e6;
    m_chisqSmooth      = -1e6;
    m_rmsResid         = 0.;
    m_Q                = -1e6;

    m_KalEnergy        = 0.;
    m_KalThetaMS       = 0.;
    m_Xgaps            = 0;
    m_Ygaps            = 0;
    m_XistGaps         = 0;
    m_YistGaps         = 0;

    m_numSegmentPoints = 0;
    m_chisqSegment     = 0.;
    m_nxHits           = 0;
    m_nyHits           = 0;
    m_KalEnergyErr     = 0.;
    m_TkrCal_radlen    = 0.; 


    clear();
}


void TkrKalFitTrack::initializeBase(Status stat, int type, double energy0, Point& x0, Vector& dir)
{
    m_status   = stat;
    m_type     = type;
    m_energy0  = energy0;
    m_x0       = x0;
    m_dir      = dir;

    clear();
}

void TkrKalFitTrack::initializeQual(double chiSq, double ChiSqSmooth, double rms, double quality, double e, double e_err, double ms)
{
    m_chisq        = chiSq;
    m_chisqSmooth  = ChiSqSmooth;
    m_rmsResid     = rms;
    m_Q            = quality;
    m_KalEnergy    = e;
    m_KalEnergyErr = e_err;
    m_KalThetaMS   = ms;
}

void TkrKalFitTrack::initializeGaps(int xgaps, int ygaps, int x1st, int y1st)
{
    m_Xgaps    = xgaps;
    m_Ygaps    = ygaps;
    m_XistGaps = x1st;
    m_YistGaps = y1st;
}

void TkrKalFitTrack::initializeKal(int nSegPoints, double chisqSeg, int nxHits, int nyHits,
                                   double radlen)
{
    m_numSegmentPoints = nSegPoints;
    m_chisqSegment     = chisqSeg;
    m_nxHits           = nxHits;
    m_nyHits           = nyHits;
    m_TkrCal_radlen    = radlen;
}


TkrKalFitTrack::~TkrKalFitTrack()
{
    clear();
}

bool TkrKalFitTrack::empty(int numHits) const
{
    bool empty = false;
    if (getLayer()     < 0)       empty = true;
    if (getNumHits()   < numHits) empty = true;
    if (getChiSquare() < 0.)      empty = true;

    return empty;
}

Vector TkrKalFitTrack::getDirection(TrackEnd end) const
{
    TkrFitPar trkPar = getFoLPlane(end).getHit(TkrFitHit::SMOOTH).getPar();

    return Vector(-trkPar.getXSlope(),-trkPar.getYSlope(),-1.).unit();
}

void TkrKalFitTrack::writeOut(MsgStream& log) const
{
    
    log << MSG::DEBUG;
    if (log.isActive() ) {
        log << " --- TkrKalFitTrack::writeOut --- "            << endreq
            << " Position      = " << getPosition().x()  << " " 
            << getPosition().y()  << " " << getPosition().z()  << endreq
            << " Direction     = " << getDirection().x() << " " 
            << getDirection().y() << " " << getDirection().z() << endreq
            << " Energy        = " << getEnergy() << endreq
            << " first Layer   = " << getLayer() << endreq
            << " Tower         = " << getTower() << endreq
            << " quality       = " << getQuality()       << endreq
            << " num m_hits    = " << getNumHits();
    }
    log << endreq;
}


TkrFitPlane TkrKalFitTrack::getFoLPlane(TrackEnd end) const
{
    if (size() == 0) 
    {
        return TkrFitPlane();
    }
    else
    {
        if (end == TkrFitTrackBase::Start) return front();
        else                               return back();
    }
}

Ray TkrKalFitTrack::getRay(TrackEnd) const 
{
    return Ray(getPosition(),getDirection());
}


double TkrKalFitTrack::getKink(int kplane) const
{
  // Purpose and Method: Computes the 3D angle between two
  //         track segments (3 pairs of (x,y) hits)
  // Inputs: plane where the kink is located
  // Outputs: the kink angle
  // Dependencies: None
  // Restrictions and Caveats:  None

    //Only valid past first 2 planes 
    if(kplane < 2 )
        return 10.; // Rogue value  - Invalid arg.

    int nplanes = size();
    if(kplane + 2 > nplanes) //Same deal - last planes not valid
        return 10.; 

    TkrFitPlaneConPtr planeIter = begin();

    Vector t0(0.,0.,0.);
    Vector t1(0.,0.,0.);
    int old_Plane_Id = planeIter[kplane-2].getIDPlane(); 
    double sX = 0.;
    double sY = 0.; 
    double count = 0.;
    
    for (int iplane = kplane-2; iplane < nplanes; iplane++) {
          
        if(planeIter[iplane].getIDPlane() == old_Plane_Id) {
            sX += planeIter[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY += planeIter[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            count += 1.; 
            continue; 
        }    
        t1 = Vector(-sX/count, -sY/count, -1.).unit();
        
        if(t0.magnitude() < .1) { // Trap for first plane
            t0 = t1; 
            sX = planeIter[iplane].getHit(TkrFitHit::SMOOTH).getPar().getXSlope();
            sY = planeIter[iplane].getHit(TkrFitHit::SMOOTH).getPar().getYSlope();
            count =1.;
            old_Plane_Id = planeIter[iplane].getIDPlane();           
            continue;
        }
        else break; 
    }
        
    double t0t1 = t0*t1;
    double kink_theta = acos(t0t1);
    return kink_theta;
}

double TkrKalFitTrack::getKinkNorma(int kplane) const
{
  // Purpose and Method: Computes the normalized 3D kink angle.
  //          (see KalFitTrack::kink()). Normalization is based on
  //          expected multiple scattering
  // Inputs: plane where kink is located
  // Outputs: no. of sigmas for the measured kink
  // Dependencies: None
  // Restrictions and Caveats:  None

    double kink_angle = getKink(kplane); 
    if(kink_angle > 3.) return kink_angle;

    TkrFitPlaneConPtr planeIter = begin();

    double rad_len = planeIter[kplane].getRadLen();
    if(planeIter[kplane].getIDPlane() == planeIter[kplane-1].getIDPlane()){
        rad_len += planeIter[kplane-1].getRadLen();
    } //Assume its the next one
    else{
        rad_len += planeIter[kplane+1].getRadLen();
    }
    double thetaMS_prj =  13.6/getEnergy() * sqrt(rad_len) * 
                           (1. + .038*log(rad_len)); 
    double sigma_kink = kink_angle /(1.414*thetaMS_prj); 
    return sigma_kink;
}
