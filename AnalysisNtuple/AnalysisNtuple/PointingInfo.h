/** @file PointingInfo.h
// $Header$
*/

#ifndef PointingInfo_h
#define PointingInfo_h

#include <string>

class INTupleWriterSvc;
class MsgStream;

namespace astro { 
    class GPS;
    class PointingHistory;
}

/** @class PointingInfo
@brief Manage the data in the FT2 definition
*/

namespace AnalysisNtuple {

class PointingInfo {
public:
    PointingInfo(){}

    //!@brief  associate it with the the Pt part of the "merit" tuple
    void setPtTuple(INTupleWriterSvc* tuple, const std::string& tname);

    //! fill the pointing info for the current orbital status
    void execute( const astro::GPS& gps);

    void setHistoryFile( const std::string filename);

private:
    double start, stop;
    float sc_position[3];
    float lat_geo, lon_geo;
    float lat_mag;
    float rad_geo;
    float ra_zenith, dec_zenith;
    float ra_scz, dec_scz;
    float ra_scx, dec_scx;
    float in_saa;
    float livetime;
    float L; ///< McIllwain L parameter
    float B; ///< magnetic field
    float zenith_scz; ///< space craft zenith angle

    float lambda;
    float R;
    float bEast;
    float bNorth;
    float bUp;

    int   lat_mode;
    int   lat_config;
    int   data_qual;
    float rock_angle;
    float livetime_frac;
    float ft2_start;
    float ft2_stop;
    //int in_saa  // already have this?
    //astro::PointingHistory  m_history;
    std::string             m_filename;
};

}

#endif 
