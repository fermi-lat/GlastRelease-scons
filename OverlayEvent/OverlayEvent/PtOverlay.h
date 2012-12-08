#ifndef PtOverlay_h
#define PtOverlay_h

#include <iostream>
#include <vector>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

/**
* @class PtOverlay
* @brief Overlay TDS for GEM data
*
* $Header$
*/
static const CLID& CLID_PtOverlay = InterfaceID("PtOverlay", 1, 0);

namespace Event {


/** @class PtOverlay
  * @brief Local storage of Pt variables from input overlay event data
  * $Header$
*/
class PtOverlay : public DataObject
{
public:

    PtOverlay() { clear(); };
    virtual ~PtOverlay() { clear(); };
    virtual std::ostream& fillStream(std::ostream &s) const;
    friend std::ostream& operator << (std::ostream &s, const PtOverlay &obj);

    /// Retrieve reference to class definition structure
    virtual const CLID& clID() const   { return PtOverlay::classID(); }
    static const  CLID& classID()      { return CLID_PtOverlay; }

    void clear();

    void initPtOverlay(const double start,
                       const float* sc_position,
                       const float  lat_geo,
                       const float  lon_geo,
                       const float  lat_mag,
                       const float  rad_geo,
                       const float  ra_scz,
                       const float  dec_scz,
                       const float  ra_scx,
                       const float  dec_scx,
                       const float  zenith_scx,
                       const float  B,
                       const float  L,
                       const float  lambda,
                       const float  R,
                       const float  bEast,
                       const float  bNorth,
                       const float  bUp,

                       const int    lat_mode,
                       const int    lat_config,
                       const int    data_qual,
                       const float  rock_angle,
                       const float  livetime_frac

                      );

    const double getStartTime()   const {return m_start;}
    const float* getSC_Position() const {return &m_sc_position[0];}
    const float  getLatGeo()      const {return m_lat_geo;}
    const float  getLonGeo()      const {return m_lon_geo;}
    const float  getLatMag()      const {return m_lat_mag;}
    const float  getRadGeo()      const {return m_rad_geo;}
    const float  getRaScz()       const {return m_ra_scz;}
    const float  getDecScz()      const {return m_dec_scz;}
    const float  getRaScx()       const {return m_ra_scx;}
    const float  getDecScx()      const {return m_dec_scx;}
    const float  getZenithScz()   const {return m_zenith_scz;}
    const float  getB()           const {return m_B;}
    const float  getL()           const {return m_L;}
    const float  getLambda()      const {return m_lambda;}
    const float  getR()           const {return m_R;}
    const float  getBEast()       const {return m_bEast;}
    const float  getBNorth()      const {return m_bNorth;}
    const float  getBUp()         const {return m_bUp;}

	const int    getLATMode()      const {return m_lat_mode;}
	const int    getLATConfig()    const {return m_lat_config;}
	const int    getDataQual()     const {return m_data_qual;}
	const float  getRockAngle()    const {return m_rock_angle;}
	const float  getLivetimeFrac() const {return m_livetime_frac;}   

    void setStartTime(double start)         {m_start          = start;}
    void setSC_Position(float* sc_position) {m_sc_position[0] = sc_position[0];
                                             m_sc_position[1] = sc_position[1];
                                             m_sc_position[2] = sc_position[2];}
    void setLatGeo(float lat_geo)           {m_lat_geo        = lat_geo;}
    void setLonGeo(float lon_geo)           {m_lon_geo        = lon_geo;}
    void setLatMag(float lat_mag)           {m_lat_mag        = lat_mag;}
    void setRadGeo(float rad_geo)           {m_rad_geo        = rad_geo;}
    void setRaScz(float ra_scz)             {m_ra_scz         = ra_scz;}
    void setDecScz(float dec_scz)           {m_dec_scz        = dec_scz;}
    void setRaScx(float ra_scx)             {m_ra_scx         = ra_scx;}
    void setDecScx(float dec_scx)           {m_dec_scx        = dec_scx;}
    void setZenithScz(float zenith_scz)     {m_zenith_scz     = zenith_scz;}
    void setB(float B)                      {m_B              = B;}
    void setL(float L)                      {m_L              = L;}
    void setLambda(float lambda)            {m_lambda         = lambda;}
    void setR(float R)                      {m_R              = R;}
    void setBEast(float bEast)              {m_bEast          = bEast;}
    void setBNorth(float bNorth)            {m_bNorth         = bNorth;}
    void setBUp(float bUp)                  {m_bUp            = bUp;}

    void setLATMode(int lat_mode)           {m_lat_mode       = lat_mode;}
    void setLATConfig(int lat_config)       {m_lat_config     = lat_config;}
    void setDataQual(int data_qual)         {m_data_qual      = data_qual;}
    void setRockAngle(float rock_angle)     {m_rock_angle     = rock_angle;}
    void setLivetimeFrac(float livetime_frac) {m_livetime_frac = livetime_frac;}

 private:

    // Variables to be used to calculate the Pt variables and stored into the ntuple
    double m_start;
    float  m_sc_position[3];
    float  m_lat_geo;
    float  m_lon_geo;
    float  m_lat_mag;
    float  m_rad_geo;
    float  m_ra_scz;
    float  m_dec_scz;
    float  m_ra_scx; 
    float  m_dec_scx;
    float  m_zenith_scz;       ///< space craft zenith angle
    float  m_B;                ///< magnetic field
    float  m_L;                ///< McIllwain L parameter

    float  m_lambda;
    float  m_R;
    float  m_bEast;
    float  m_bNorth;
    float  m_bUp;

    int    m_lat_mode;
    int    m_lat_config;
    int    m_data_qual;
    float  m_rock_angle;
    float  m_livetime_frac;
 };


inline void PtOverlay::initPtOverlay(const double start,
                                     const float* sc_position,
                                     const float  lat_geo,
                                     const float  lon_geo,
                                     const float  lat_mag,
                                     const float  rad_geo,
                                     const float  ra_scz,
                                     const float  dec_scz,
                                     const float  ra_scx,
                                     const float  dec_scx,
                                     const float  zenith_scz,
                                     const float  B,
                                     const float  L,
                                     const float  lambda,
                                     const float  R,
                                     const float  bEast,
                                     const float  bNorth,
                                     const float  bUp,

                                     const int    lat_mode,
                                     const int    lat_config,
                                     const int    data_qual,
                                     const float  rock_angle,
                                     const float  livetime_frac
                                    )
{
    m_start          = start;
    m_sc_position[0] = sc_position[0];
    m_sc_position[1] = sc_position[1];
    m_sc_position[2] = sc_position[2];
    m_lat_geo        = lat_geo;
    m_lon_geo        = lon_geo;
    m_lat_mag        = lat_mag;
    m_rad_geo        = rad_geo;
    m_ra_scz         = ra_scz;
    m_dec_scz        = dec_scz;
    m_ra_scx         = ra_scx;
    m_dec_scx        = dec_scx;
    m_zenith_scz     = zenith_scz;
    m_B              = B;
    m_L              = L;
    m_lambda         = lambda;
    m_R              = R;
    m_bEast          = bEast;
    m_bNorth         = bNorth;
    m_bUp            = bUp;

    m_lat_mode       = lat_mode;
    m_lat_config     = lat_config;
    m_data_qual      = data_qual;
    m_rock_angle     = rock_angle;
    m_livetime_frac  = livetime_frac;
}

inline void PtOverlay::clear() 
{
    m_start          = 0.;
    m_sc_position[0] = 0.;
    m_sc_position[1] = 0.;
    m_sc_position[2] = 0.;
    m_lat_geo        = 0.;
    m_lon_geo        = 0.;
    m_lat_mag        = 0.;
    m_rad_geo        = 0.;
    m_ra_scz         = 0.;
    m_dec_scz        = 0.;
    m_zenith_scz     = 0.;
    m_B              = 0.;
    m_L              = 0.;
    m_lambda         = 0.;
    m_R              = 0.;
    m_bEast          = 0.;
    m_bNorth         = 0.;
    m_bUp            = 0.;
    m_lat_mode       = 0;
    m_lat_config     = 0;
    m_data_qual      = 0;
    m_rock_angle     = 0;
    m_livetime_frac  = 0;
}

inline std::ostream& PtOverlay::fillStream(std::ostream &s) const{
    s << "PtOverlay:" <<std::endl;
    s << "Start Time = " <<  getStartTime() << std::endl;
    s << "Position: " << getSC_Position()[0] << ", " << getSC_Position()[1] << ", " << getSC_Position()[2] << std::endl;
    s << "LAT Geo: " << m_lat_geo << ", LONG Geo: " << m_lon_geo << std::endl;
    s << "LAT Mag: " << m_lat_mag <<  std::endl;
    s << "Rad Geo: " << m_rad_geo << std::endl;
    s << "RA scz: " << m_ra_scz << ", DEC scz: " << m_dec_scz << ", Zenith scz: " << m_zenith_scz << std::endl;
    s << "B: " << m_B << ", L: " << m_L << ", lambda: " << m_lambda << ", R: " << m_R << std::endl;
    s << "bEast: " << m_bEast << ", bNorth: " << m_bNorth << ", bUp: " << m_bUp << std::endl;
    s << "lat_mode: " << m_lat_mode << ", lat_config: " << m_lat_config << ", data_qual: "
      << m_data_qual << ", rock_angle: " << m_rock_angle 
      << ", livetime_frac: " << m_livetime_frac << std::endl;
    return s;
}


inline std::ostream& operator<<(std::ostream &s, const PtOverlay &obj){
    return obj.fillStream(s);
}


}//namespace Event


#endif
