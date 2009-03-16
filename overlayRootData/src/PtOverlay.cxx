#include "overlayRootData/PtOverlay.h"
#include "Riostream.h"

ClassImp(PtOverlay)

//
// Class:  PtOverlay        
// Description: Digitizations for Cal
// Based on the PtOverlay TDS class in GlastEvent written by J. Eric Grove.
//              
// Author:  Heather Kelly
//

PtOverlay::PtOverlay() 
{
    Clear();

    return;
}

PtOverlay::~PtOverlay() 
{ 
    return;
}

void PtOverlay::initialize(const Double_t  start,
                           const TVector3& sc_position,
                           const Float_t   lat_geo,
                           const Float_t   lon_geo,
                           const Float_t   lat_mag,
                           const Float_t   rad_geo,
                           const Float_t   ra_scz,
                           const Float_t   dec_scz,
                           const Float_t   ra_scx,
                           const Float_t   dec_scx,
                           const Float_t   zenith_scz,
                           const Float_t   B,
                           const Float_t   L,
                           const Float_t   lambda,
                           const Float_t   R,
                           const Float_t   bEast,
                           const Float_t   bNorth,
                           const Float_t   bUp) 
{
    m_start       = start;
    m_sc_position = sc_position;
    m_lat_geo     = lat_geo;
    m_lon_geo     = lon_geo;
    m_lat_mag     = lat_mag;
    m_rad_geo     = rad_geo;
    m_ra_scz      = ra_scz;
    m_dec_scz     = dec_scz;
    m_zenith_scz  = zenith_scz;
    m_B           = B;
    m_L           = L;
    m_lambda      = lambda;
    m_R           = R;
    m_bEast       = bEast;
    m_bNorth      = bNorth;
    m_bUp         = bUp;
}

void PtOverlay::Clear(Option_t *option) {
    m_start       = 0.;
    m_sc_position = TVector3(0.,0.,0.);
    m_lat_geo     = 0.;
    m_lon_geo     = 0.;
    m_lat_mag     = 0.;
    m_rad_geo     = 0.;
    m_ra_scz      = 0.;
    m_dec_scz     = 0.;
    m_zenith_scz  = 0.;
    m_B           = 0.;
    m_L           = 0.;
    m_lambda      = 0.;
    m_R           = 0.;
    m_bEast       = 0.;
    m_bNorth      = 0.;
    m_bUp         = 0.;
}

void PtOverlay::Print(Option_t *option) const 
{
    TObject::Print(option);
    cout.precision(2);
    std::cout << "Start time: " << m_start << endl;

    std::cout << "SC position: ";
    m_sc_position.Print();
    std::cout << std::endl;
}
