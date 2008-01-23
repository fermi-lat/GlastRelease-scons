#ifndef AcdCalibMgr_H
#define AcdCalibMgr_H
// $Header$

// LOCAL
#include "AcdUtil/AcdCalib.h"
#include "AcdCalibSvcBase.h"

// GLAST
#include "CalibData/Acd/AcdCalibObj.h"
#include "CalibData/Acd/AcdCalibBase.h"
#include "CalibSvc/ICalibPathSvc.h"

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
  
  AcdCalibMgr(ICalibPathSvc::CalibItem calibItem)
    :m_calibTypePath(""),
     m_calibItem(calibItem),
     m_owner(0),
     m_ideal(false),
     m_isValid(false),
     m_serNo(-1)
  {}
  
  virtual ~AcdCalibMgr() {};

  /// return the type of calibration handled by this manager
  virtual AcdCalibData::CALTYPE calibType() const = 0;
  
  /// initialize and register this manager
  StatusCode initialize(const std::string &flavor,
			AcdCalibSvcBase &acs);

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

  /// return the calibration
  inline CalibData::AcdCalibBase* calibBase() {
    return m_calibBase;
  }

  /// are we just returning the ideal calibration instead of using the DB
  inline bool ideal() const {
    return m_ideal;
  }

 private:

  /// TDS path to calib data for my calib_type
  std::string                  m_calibTypePath;

  /// TDS path to calib data for my calib_type and path
  std::string                  m_calibPath;
  ICalibPathSvc::CalibItem     m_calibItem;

  /// TDS location for root of my calib_type and path
  CalibData::AcdCalibBase     *m_calibBase;

  AcdCalibSvcBase*          m_owner;

  /// calib flavor
  std::string               m_flavor;     
  /// Should we use the "ideal" calibration instead of the database
  bool                      m_ideal;
  /// validity state of AcdCalibMgr data
  bool                      m_isValid;    
  /// serial # for current calibration source
  int                       m_serNo;      

};

/** @class AcdCalibMgrTmpl
    @author Eric Charles 
    \brief template for handling specific calibration type

*/

template <class T>
class AcdCalibMgrTmpl : public AcdCalibMgr {
  
public:

  /// the Type of object managed by this calibration manager
  typedef typename T::ObjType CalibObjType;
  
public:

  /// Standard c'tor just specifies which calibration type
  AcdCalibMgrTmpl():
    AcdCalibMgr( AcdCalib::calibItem( T::calibType() ) ){;}

  /// D'tor is a no-op
  virtual ~AcdCalibMgrTmpl(){;}

  /// Get a calibration
  StatusCode getCalibration(idents::AcdId id, unsigned pmt, CalibObjType*& calib) {

    if ( ideal() ) {
      // get the ideal value
      CalibData::AcdCalibObj* calibObj = AcdCalib::getIdeal(CalibObjType::calibType(), id, pmt);
      calib = static_cast<CalibObjType*>(calibObj);
      return StatusCode::SUCCESS;
    }  
    
    StatusCode sc = updateCalib();
    if (sc.isFailure()) {
      // null and return failure code
      calib = 0;
      return sc;
    }
    
    T* calibCol = static_cast<T*>(calibBase());
    if ( calibCol == 0 ) {
      return StatusCode::FAILURE;
    }
    
    calib = calibCol->getPmt(id,pmt);
    if ( calib == 0 ) return StatusCode::FAILURE;
    return StatusCode::SUCCESS;	
  }
  
  /// use typename to return the calibration type
  virtual AcdCalibData::CALTYPE calibType() const {
    return T::calibType();
  }
  
};


#endif
