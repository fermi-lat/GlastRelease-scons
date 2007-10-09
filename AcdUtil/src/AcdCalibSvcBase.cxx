// $Header $
/** @file
    @author Zach Fewtrell
 */
// @file
//
//
// Author: Zachary Fewtrell

// LOCAL
#include "AcdCalibSvcBase.h"

#include "AcdCalibMgr.h"

// GLAST

// EXTLIB
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Incident.h"
#include "GaudiKernel/StatusCode.h"

// STD


AcdCalibSvcBase::AcdCalibSvcBase() 
  :m_mgrs( AcdCalibData::NDESC, std::pair<AcdCalibMgr*, const StringProperty*>(0,0) )
{
}

AcdCalibSvcBase::~AcdCalibSvcBase() {
  for ( std::vector<std::pair< AcdCalibMgr*, const StringProperty*> >::iterator itr = m_mgrs.begin(); itr!= m_mgrs.end(); itr++ ) {
    delete itr->first;
    delete itr->second;
  }
}

StatusCode AcdCalibSvcBase::getCalibrationMgr(AcdCalibData::CALTYPE type, AcdCalibMgr*& calibMgr) {
  calibMgr = m_mgrs[type].first;
  return calibMgr != 0 ? StatusCode::SUCCESS : StatusCode::FAILURE;
}


/// Add a manager (and the associated flavor)
void AcdCalibSvcBase::addMgr(AcdCalibData::CALTYPE type, AcdCalibMgr* mgr, const StringProperty* flavor) {
  m_mgrs[type] = std::make_pair(mgr,flavor);
}


StatusCode AcdCalibSvcBase::prepapreManagers(MsgStream& log, const std::string& defaultFlavor) {

  for ( std::vector<std::pair< AcdCalibMgr*, const StringProperty*> >::iterator itr = m_mgrs.begin(); itr != m_mgrs.end(); itr++ ) {
    AcdCalibMgr* mgr = itr->first;
    if ( mgr == 0 ) continue;    
    const StringProperty* flavorProp = itr->second;
    std::string flavorName = flavorProp->value();
    flavorName = flavorName.length() == 0 ?  defaultFlavor : flavorName;    
    log << MSG::DEBUG << "  " << flavorProp->name() << "\t\t" << flavorName << endreq;
    StatusCode sc = mgr->initialize( flavorName, *this );
    if ( sc.isFailure() ) {
      return sc;
    }
  }
  return StatusCode::SUCCESS;
}

/// Inform that a new incident has occured
void AcdCalibSvcBase::handle ( const Incident& inc ) { 
  if ((inc.type() == "BeginEvent")) {
    for ( std::vector<std::pair< AcdCalibMgr*, const StringProperty*> >::iterator itr = m_mgrs.begin(); itr!= m_mgrs.end(); itr++ ) {
      AcdCalibMgr* mgr = itr->first;
      if ( mgr == 0 ) continue;
      mgr->invalidate();
    }
  }
  return; 
}

