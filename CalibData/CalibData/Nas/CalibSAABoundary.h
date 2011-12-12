// $Header$
#ifndef CalibData_CalibSAABoundary_h
#define CalibData_CalibSAABoundary_h

/** @class CalibTest1  

  Very simple calibration data-like class to be used for testing 
  calibration infrastructure

  @author J. Bogart
*/

#include "CalibData/CalibBase.h"
#include "CalibData/CalibModel.h"
#include <vector>

// extern const CLID& CLID_Calib_NAS_SAABoundary;


namespace CalibData {
  class CalibSAABoundary : public CalibBase {
  public:
    CalibSAABoundary(std::vector<double>& lat, std::vector<double>& lon,
                     const Gaudi::Time& since, const Gaudi::Time& till,int serNo);

    CalibSAABoundary(const CalibSAABoundary& other);

    virtual ~CalibSAABoundary() {}

    // Re-implemented from DataObject
    inline virtual const CLID& clID() const { return classID(); } 
    
    inline static  const CLID& classID() { return CLID_Calib_NAS_SAABoundary; };


    // Re-implemented from CalibBase
    virtual StatusCode   update(CalibBase& other, MsgStream* log);

    const std::vector<double>& getSAAPolygonLatitudes() const {return m_lat;};
    const std::vector<double>& getSAAPolygonLongitudes() const {return m_lon;};
    bool getFirstVertex(std::pair<double,double>& vertex);
    bool getNextVertex(std::pair<double,double>& vertex);
    
    
  protected:

  private:
    std::vector<double> m_lat;
    std::vector<double> m_lon;
    std::vector<double>::const_iterator m_lat_iterator;
    std::vector<double>::const_iterator m_lon_iterator;
  };
}


#endif
