/** @file PointingInfo.h
// $Header$
*/

#ifndef PointingInfo_h
#define PointingInfo_h

#include <string>

class INTupleWriterSvc;
namespace Event { class Exposure;} 

/** @class PointingInfo
Manage the data in the FT2 definition

 the root tuple, a la FT2
 http://glast.gsfc.nasa.gov/ssc/dev/fits_def/definitionFT2.html
*/

class PointingInfo {
public:
    PointingInfo(){}

    //!@brief  associate it with the the FT2 tuple
    void setFT2Tuple(INTupleWriterSvc* tuple, const std::string& tname);

    //!@brief  associate it with the the Pt part of the "merit" tuple
    void setPtTuple(INTupleWriterSvc* tuple, const std::string& tname);

    //! fill the pointing info for the given start time
    void set(double start_time);

    //! finish it.
    void finish(double stop_time, double live);

    //! accessor for time
    double start_time()const{return start;}

    // return TDS object for old scheme
    Event::Exposure* forTDS()const;
private:
    double start, stop;
    float sc_position[3];
    float lat_geo, lon_geo;
    float rad_geo;
    float ra_zenith, dec_zenith;
    float ra_scz, dec_scz;
    float ra_scx, dec_scx;
    float livetime;

};

#endif 
