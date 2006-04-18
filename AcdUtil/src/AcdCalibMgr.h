#ifndef AcdCalibMgr_H
#define AcdCalibMgr_H
// $Header$

// LOCAL

// GLAST
#include "CalibData/RangeBase.h"
#include "CalibData/Acd/AcdCalibBase.h"

// EXTLIB
#include "GaudiKernel/IService.h"

// STD
#include <string>
#include <algorithm>

class AcdCalibSvc;

/** @class AcdCalibMgr
    @author Eric Charles (from Zach Fewtrell's CalibItemMgr)
    \brief abstract class for handling a single calibration data type.

    provides the following services: 
    - TDS calibration data retrieval & indexing
    - validation period checking
    - configurable 'flavor'
    - support for local data store.
*/
class AcdCalibMgr {
public:
  AcdCalibMgr(const std::string &calibTypePath)
    :m_calibTypePath(calibTypePath),
     m_ideal(false),
     owner(0),     
     m_isValid(false),
     m_serNo(-1)
  {}
  
  virtual ~AcdCalibMgr() {};
  
  StatusCode initialize(const std::string &flavor,
                        const AcdCalibSvc &acs);

  /// data should be invalidated at beginning of each event.
  /// just in case there is a change in validity period
  void invalidate() {m_isValid = false;} 

protected:

  /** \brief check calib validity period, (re)build local store if necessary

  needs to be called once per event (_before_ processing calibration data ;).
  Subsequent calls in same event will return immediately.

  */
  virtual StatusCode updateCalib();          

  /// generate full set of local data (if applicable for calib_type)
  virtual StatusCode genLocalStore() {
    // default to a no-op
    return StatusCode::SUCCESS;
  }

  /// TDS path to calib data for my calib_type
  std::string                  m_calibTypePath;

  /// TDS path to calib data for my calib_type and path
  std::string                  m_calibPath;

  /// TDS location for root of my calib_type and path
  CalibData::AcdCalibBase     *m_calibBase;

  /// Should we use the "ideal" calibration instead of the database
  bool                         m_ideal;

  /// ref to owner->CalCalibSvc object
  const AcdCalibSvc           *owner;

  /** retrieve spec'd rangeBase object, update if necessary
      \return NULL if there is no data 
  */
  CalibData::RangeBase *getPmt(idents::AcdId id, unsigned pmt) {
    return m_calibBase->getPmt(id,pmt);
  }
  
 private:

  /// calib flavor
  std::string       m_flavor;     
  /// validity state of AcdCalibMgr data
  bool              m_isValid;    
  /// serial # for current calibration source
  int               m_serNo;      

};
#endif
