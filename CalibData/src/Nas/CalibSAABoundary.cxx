// $Header$

/** @class CalibSAABoundary
 *    Implementation of a calibration class for SAA boundary descirption
 */

#include "CalibData/Nas/CalibSAABoundary.h"
#include "GaudiKernel/MsgStream.h"

namespace CalibData {
    CalibSAABoundary:: CalibSAABoundary(std::vector<double>& lat, std::vector<double>& lon,
                                        const Gaudi::Time& since, const Gaudi::Time& till,int serNo):
         CalibBase(since, till, serNo), m_lat(lat), m_lon(lon)
       {
    //    m_me = this;
       }  


  StatusCode CalibSAABoundary::update(CalibBase& other, MsgStream* log) {
    // The following dynamic_cast has got to work
    CalibSAABoundary& other1 = dynamic_cast<CalibSAABoundary& >(other);

    CalibBase::update(other1, log);
    m_lat = other1.m_lat;
    m_lon = other1.m_lon;

    return StatusCode::SUCCESS;
  }

  CalibSAABoundary::CalibSAABoundary(const CalibSAABoundary& other) : 
         IValidity(),CalibBase(other), m_lat(other.m_lat), m_lon(other.m_lon) {
  }

  bool CalibSAABoundary::getFirstVertex(std::pair<double,double>& vertex){
    m_lat_iterator=m_lat.begin();
    m_lon_iterator=m_lon.begin();

    if(m_lat_iterator==m_lat.end() || m_lon_iterator==m_lon.end()) return false;   
    vertex.first= *m_lat_iterator;
    vertex.second= *m_lon_iterator;
    return true;
  };
  
  bool CalibSAABoundary::getNextVertex(std::pair<double,double>& vertex){
    ++m_lat_iterator;
    ++m_lon_iterator;
    
    if(m_lat_iterator==m_lat.end() || m_lon_iterator==m_lon.end()) return false;   
    vertex.first= *m_lat_iterator;
    vertex.second= *m_lon_iterator;
    return true;
  };
}
