//====================================================================
//	DbEventCnvSvc.cpp
//--------------------------------------------------------------------
//
//	Package    : Db Converters
//
//	Author     : Markus Frank
//  History    :
// +---------+----------------------------------------------+---------
// |    Date |                 Comment                      | Who     
// +---------+----------------------------------------------+---------
// | 21/04/99| Initial version.                             | MF
// +---------+----------------------------------------------+---------
//====================================================================
#define  DBCNV_DbEVENTCNVSVC_CPP
// $Header$

// Include files
#include "Gaudi/Interfaces/ISvcLocator.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/Kernel/SvcFactory.h"
#include "DbEventCnvSvc.h"

static const SvcFactory<DbEventCnvSvc> s_DbEvtCnvSvcFactory;
const ISvcFactory& DbEventCnvSvcFactory = s_DbEvtCnvSvcFactory;

/// Initialize the service.
StatusCode DbEventCnvSvc::initialize()     {
  StatusCode status = DbCnvSvc::initialize();
  if ( status.isSuccess() )   {
    ISvcLocator* svcLocator = serviceLocator();
    //======================================================================
    // Add known converters to the service: 
    //======================================================================
    IDataProviderSvc *pIDP = 0;
    // Set event data service
    status = svcLocator->getService("EventDataSvc", IID_IDataProviderSvc, (IInterface*&)pIDP);
    if ( !status.isSuccess() ) return status;

    // Set proper data store
    status = setStore ( pIDP );
    if ( !status.isSuccess() ) return status;
  }
  return status;
}

/// Update state of the service
StatusCode DbEventCnvSvc::updateServiceState(IOpaqueAddress* /* pAddress */ )    {
  return StatusCode::SUCCESS;
}

/// Standard Constructor
DbEventCnvSvc::DbEventCnvSvc(const std::string& name, ISvcLocator* svc)
 : DbCnvSvc(name, svc)
{
}

/// Standard Destructor
DbEventCnvSvc::~DbEventCnvSvc()   {
}
