
#ifndef overlayRootData_PtOverlay_H
#define overlayRootData_PtOverlay_H 1

#include "commonRootData/idents/CalXtalId.h"
#include "TObject.h"
#include "TVector3.h"

/** @class PtOverlay
 * @brief Root CAL data object for overlays
 *
 * @author Tracy Usher
 * $Header$
*/

class PtOverlay : public TObject  { 
    
public:
    
    PtOverlay();
    
    virtual ~PtOverlay();

    void Clear(Option_t *option ="");

    void initialize(const Double_t  start,
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
                    const Float_t   bUp);

    void Print(Option_t *option="") const;

    /** @ingroup PtOverlayGroup */
    /*@{*/
    /// Retrieve log identifier
    const Double_t  getStartTime()   const {return m_start;}
    const TVector3& getSC_Position() const {return m_sc_position;}
    const Float_t  getLatGeo()       const {return m_lat_geo;}
    const Float_t  getLonGeo()       const {return m_lon_geo;}
    const Float_t  getLatMag()       const {return m_lat_mag;}
    const Float_t  getRadGeo()       const {return m_rad_geo;}
    const Float_t  getRaScz()        const {return m_ra_scz;}
    const Float_t  getDecScz()       const {return m_dec_scz;}
    const Float_t  getRaScx()        const {return m_ra_scx;}
    const Float_t  getDecScx()       const {return m_dec_scx;}
    const Float_t  getZenithScz()    const {return m_zenith_scz;}
    const Float_t  getB()            const {return m_B;}
    const Float_t  getL()            const {return m_L;}
    const Float_t  getLambda()       const {return m_lambda;}
    const Float_t  getR()            const {return m_R;}
    const Float_t  getBEast()        const {return m_bEast;}
    const Float_t  getBNorth()       const {return m_bNorth;}
    const Float_t  getBUp()          const {return m_bUp;}

    void setStartTime(Double_t start)                {m_start       = start;}
    void setSC_Position(const TVector3& sc_position) {m_sc_position = sc_position;}
    void setLatGeo(Float_t lat_geo)                  {m_lat_geo     = lat_geo;}
    void setLonGeo(Float_t lon_geo)                  {m_lon_geo     = lon_geo;}
    void setLatMag(Float_t lat_mag)                  {m_lat_mag     = lat_mag;}
    void setRadGeo(Float_t rad_geo)                  {m_rad_geo     = rad_geo;}
    void setRaScz(Float_t ra_scz)                    {m_ra_scz      = ra_scz;}
    void setDecScz(Float_t dec_scz)                  {m_dec_scz     = dec_scz;}
    void setRaScx(Float_t ra_scx)                    {m_ra_scx      = ra_scx;}
    void setDecScx(Float_t dec_scx)                  {m_dec_scx     = dec_scx;}
    void setZenithScz(Float_t zenith_scz)            {m_zenith_scz  = zenith_scz;}
    void setB(Float_t B)                             {m_B           = B;}
    void setL(Float_t L)                             {m_L           = L;}
    void setLambda(Float_t lambda)                   {m_lambda      = lambda;}
    void setR(Float_t R)                             {m_R           = R;}
    void setBEast(Float_t bEast)                     {m_bEast       = bEast;}
    void setBNorth(Float_t bNorth)                   {m_bNorth      = bNorth;}
    void setBUp(Float_t bUp)                         {m_bUp         = bUp;}
    /*@}*/

private:
    
    // Variables to be used to calculate the Pt variables and stored into the ntuple
    Double_t m_start;
    TVector3 m_sc_position;
    Float_t  m_lat_geo;
    Float_t  m_lon_geo;
    Float_t  m_lat_mag;
    Float_t  m_rad_geo;
    Float_t  m_ra_scz;
    Float_t  m_dec_scz;
    Float_t  m_ra_scx; 
    Float_t  m_dec_scx;
    Float_t  m_zenith_scz;       ///< space craft zenith angle
    Float_t  m_B;                ///< magnetic field
    Float_t  m_L;                ///< McIllwain L parameter

    Float_t  m_lambda;
    Float_t  m_R;
    Float_t  m_bEast;
    Float_t  m_bNorth;
    Float_t  m_bUp;
    
    ClassDef(PtOverlay,1) // Pt variables Overlay Class
};

#endif
