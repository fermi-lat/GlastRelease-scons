/** @file IPointingInfo.h

    @brief declaration of the IPointingInfo class

$Header$

*/

#ifndef IPointingInfo_h
#define IPointingInfo_h

namespace Event { class Exposure;} 

#include "GaudiKernel/IAlgTool.h"

/** @class IPointingInfo
    @brief Provides interface to the FluxPointingInfoTool to calculate the "Pt" variables
    @author Tracy Usher
*/

static const InterfaceID IID_IPointingInfo("IPointingInfo", 1 , 0);

class IPointingInfo : virtual public IAlgTool
{
public:

    //!@brief  Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IPointingInfo; }

    //!@brief  fill the pointing info for the current orbital status
    virtual void set() = 0;

    //!@brief  finish it. 
    virtual void finish(double stop_time, double live) = 0;

    //!@brief  accessor for time
    virtual double start_time() const = 0;

    //!@brief  return TDS object for old scheme
    virtual Event::Exposure* forTDS() const = 0;

    //!@brief  associate it with the the FT2 tuple
    virtual void setFT2Tuple(const std::string& tname) = 0;

    //!@brief  associate it with the the Pt part of the "merit" tuple
    virtual void setPtTuple(const std::string& tname) = 0;

    //!@brief Provide ability to read all of PointingInfo data members
    virtual const double get_start()       const = 0;
    virtual const double get_stop()        const = 0;
    virtual const float* get_sc_position() const = 0;
    virtual const float  get_lat_geo()     const = 0;
    virtual const float  get_lon_geo()     const = 0;
    virtual const float  get_lat_mag()     const = 0;
    virtual const float  get_rad_geo()     const = 0;
    virtual const float  get_ra_zenith()   const = 0;
    virtual const float  get_dec_zenith()  const = 0;
    virtual const float  get_ra_scz()      const = 0;
    virtual const float  get_dec_scz()     const = 0;
    virtual const float  get_ra_scx()      const = 0;
    virtual const float  get_dec_scx()     const = 0;
    virtual const float  get_in_saa()      const = 0;
    virtual const float  get_livetime()    const = 0;
    virtual const float  get_L()           const = 0; 
    virtual const float  get_B()           const = 0; 
    virtual const float  get_zenith_scz()  const = 0; 

};

#endif
