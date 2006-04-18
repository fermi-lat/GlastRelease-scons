// $Header$
/** @file
    @author Eric Charles (From Zach Fewtrell's CalibItemMgr)
 */
// LOCAL
#include "AcdCalibMgr.h"
#include "AcdCalibSvc.h"

// GLAST


// EXTLIB

// STD
#include <sstream>
#include <algorithm>


using namespace std;

StatusCode AcdCalibMgr::initialize(const std::string &flavor, const AcdCalibSvc &acs) {
  
  owner = &acs;
  
  m_flavor = flavor;

  m_calibPath = m_calibTypePath + '/' + flavor;

  static std::string idealName("ideal");
  if ( flavor == idealName ) {
    MsgStream msglog(owner->msgSvc(), owner->name());     
    // else return error (can't find calib)
    msglog << MSG::WARNING << "Overriding calibration to 'ideal' from"
           << m_calibPath << ".  Calib db will not be used" << endreq;
    m_ideal = true;
  }
    
  return StatusCode::SUCCESS;
}

StatusCode AcdCalibMgr::updateCalib() {
  StatusCode sc;

  // if event is already validated return quickly
  if (m_isValid) return StatusCode::SUCCESS;
  
  // ideal is valid by definition
  if (m_ideal) {    
    m_isValid = true;    
    return StatusCode::SUCCESS;
  }

  /////////////////////////////////
  //-- CHECK TDS DATA VALIDITY --//
  /////////////////////////////////

  // Get pointer to Gain tree from TDS
  // usually this f() should return immediately
  // if it fails then we have no valid calib data
  // for the current event.
  DataObject *pObject(0);
  sc = owner->m_dataProviderSvc->retrieveObject(m_calibPath, pObject);
  if (!sc.isFailure())
    m_calibBase = (CalibData::AcdCalibBase *)(pObject);
  else {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    
    // else return error (can't find calib)
    msglog << MSG::ERROR << "can't get " 
           << m_calibPath << " from calib db" << endreq;
    return sc;  
  }

  ///////////////////////////////////////
  //-- CHECK IF TDS DATA HAS CHANGED --//
  ///////////////////////////////////////

  // check serial # to see if we're still valid.
  int curSerNo = m_calibBase->getSerNo();
  if (curSerNo != m_serNo) {
    // create MsgStream only when needed (for performance)
    MsgStream msglog(owner->msgSvc(), owner->name()); 
    msglog << MSG::INFO << "Updating " << m_calibPath << endreq;
    m_serNo = curSerNo;
        
    // generate local data (e.g. splines) if needed
    sc = genLocalStore();
    if (sc.isFailure()) return sc;
  }

  // reset 'in update' flag
  m_isValid = true;

  return StatusCode::SUCCESS;
}
