//====================================================================
//	DbCnvSvc implementation
//--------------------------------------------------------------------
//
//	Package    : System ( The LHCb Offline System)
//
//  Description: Implementation of the Db data service
//
//	Author     : M.Frank
//  History    :
// +---------+----------------------------------------------+---------
// |    Date |                 Comment                      | Who     
// +---------+----------------------------------------------+---------
// | 29/10/98| Initial version                              | MF
// +---------+----------------------------------------------+---------
//
//====================================================================
#define DBCNV_DBCNVSVC_CPP
// $Header$

// Framework include files
#include "Gaudi/Kernel/DataObject.h"
#include "Gaudi/Interfaces/IFactory.h"
#include "Gaudi/Interfaces/IObjManager.h"
#include "Gaudi/Interfaces/ICnvManager.h"
#include "Gaudi/Interfaces/ISvcLocator.h"
#include "Gaudi/Interfaces/IMessageSvc.h"
#include "Gaudi/Interfaces/IOpaqueAddress.h"
#include "Gaudi/Interfaces/IDataDirectory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/MessageSvc/MsgStream.h"

// Db data service include files
#include "DbCnv/IOODataBase.h"
#include "DbCnv/DbBaseConverter.h"
#include "DbCnv/DbCnvFactory.h"
#include "DbCnv/DbIter.h"
#include "src/DbCnv/DbAddress.h"
#include "DbCnvSvc.h"

class DbDataBaseMgr : public DbOOMs  {
  friend DbCnvSvc;
private:
  IOODataBase*      m_myDataBases[255];
  long m_refCount;
  static DbDataBaseMgr* si_instance;
  DbDataBaseMgr(IMessageSvc* svc)   {
    setInstance(getPointer);
    m_messageSvc = svc;
    m_messageSvc->addRef();
    m_refCount = 0;
    memset(m_myDataBases,0,sizeof(m_myDataBases));
    m_dataBase = m_myDataBases;
  }
  virtual ~DbDataBaseMgr()    {
    m_messageSvc->release();
    setInstance(0);
  }
  static DbDataBaseMgr* instance(IMessageSvc* svc)  {
    if ( 0 == si_instance )   {
      si_instance = new DbDataBaseMgr(svc);
    }
    si_instance->m_refCount++;
    return si_instance;
  }
  long release()   {
    long cnt = --m_refCount;
    if ( m_refCount <= 0 )    {
      delete this;
    }
    return cnt;
  }
public:
  static DbOOMs* getPointer()  {
    return si_instance;
  }
};
DbDataBaseMgr* DbDataBaseMgr::si_instance = 0;

/// Standard constructor
DbCnvSvc::DbCnvSvc(const std::string& nam, ISvcLocator* svc)
: ConversionSvc( nam, svc, DBOOMS_StorageType), 
  m_manager(0),
  m_session(DBOOMS_StorageType), 
  m_federation(DBOOMS_StorageType),
  m_output(DBOOMS_StorageType)
{
  // Set address factory for my type
  setAddressFactory(this);
	declareProperty("FddbName",    m_fddbName          = "Standard.FDDB");
	declareProperty("DbTypeName",  m_implementation    = "OODataBase");

  // I know the following is not terribly nice, but this is the only way to 
  // get hold of the data base type BEFORE the service is attached to the
  // Persistency service
  StatusCode status = setProperties();
  if ( status.isSuccess() )   {
    IOODataBase* db = 0;
    status = requestImplementation(db);
    if ( status.isSuccess() )   {
      status = db->initialize(name());
      if ( status.isSuccess() )   {
        m_type       = db->type();
        m_session    = dbHandle<DbSession>   (repSvcType());
        m_federation = dbHandle<DbFederation>(repSvcType());
        m_output     = dbHandle<DbDataBase>  (repSvcType());
        db->release();
      }
    }
  }
}

/// Standard destructor
DbCnvSvc::~DbCnvSvc()    {
}

/// Query interface
StatusCode DbCnvSvc::queryInterface(const IID& riid, void** ppvInterface)  {
  if ( IID_IDataBaseMgr == riid )  {
    *ppvInterface = (IDataBaseMgr*)this;
  }
  else  {
    // Interface is not directly availible: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}

/// Locate converter factory
const ICnvFactory* DbCnvSvc::findCnvFactory(const CLID& wanted)  {
  ICnvManager* manager = cnvManager();
  if ( 0 != manager )   {
    const ICnvFactory* fac = manager->factory(DBOOMS_StorageType, wanted);
    return fac;
  }
  return 0;
}

/// Create new Converter using factory
IConverter* DbCnvSvc::createConverter(const ICnvFactory& fac)   {
  IConverter* cnv = 0;
  IDbCnvFactory* pFac = 0;
  ICnvFactory* p = const_cast<ICnvFactory*>(&fac);
  StatusCode status = p->queryInterface(IID_IDbCnvFactory, (void**)&pFac);
  if ( status.isSuccess() )   {
    cnv = pFac->instantiate(DbOOMs::instance(repSvcType()), serviceLocator());
    pFac->release();
  }
  return cnv;
}


StatusCode DbCnvSvc::requestImplementation(IOODataBase*& implementation)    {
  IOODataBase* db = 0;
  IObjManager* objManager = 0;
  implementation = 0;
  StatusCode status = serviceLocator()->queryInterface(IID_IObjManager, (void**)&objManager);
  if ( status.isSuccess() )   {
    IObjManager::ObjIterator i, stop;
    for ( i = objManager->objBegin(), stop = objManager->objEnd(); i != stop; i++ )    {
      IFactory* f = const_cast<IFactory*>(*i);
      if ( f->ident() == m_implementation )    {
        IInterface *iface = f->instantiate(serviceLocator());
        status = iface->queryInterface(IID_IOODataBase, (void**)&db);
        if ( !status.isSuccess() )    {
          iface->release();
        }
        else  {
          implementation = db;
        }
        objManager->release(); // Must release interface after use!
        return status;
      }
    }
    objManager->release(); // Must release interface after use!
    return StatusCode::FAILURE;
  }
  return status;
}

/// Initialize the Db data persistency service
StatusCode DbCnvSvc::initialize()    {
  StatusCode status = ConversionSvc::initialize();
  MsgStream log(messageService(), name());
  if ( status.isSuccess() )   {
    status = setProperties();
    if ( status.isSuccess() )   {
      IOODataBase* db = 0;
      m_manager = DbDataBaseMgr::instance(messageService());
      status = requestImplementation(db);
      if ( status.isSuccess() )   {
        status = db->initialize(name());
        if ( status.isSuccess() )   {
          if ( repSvcType() == db->type() )   {
            if ( m_manager->attachDataBase(db) )    {
              if ( m_session.open(DbOOMs::UPDATE) )    {
                if ( m_federation.open(m_session, m_fddbName, DbOOMs::READ) )   {
                  long debugLevel = messageService()->outputLevel();
                  switch ( debugLevel )   {
                  case MSG::NIL:      
                    DbOOMs::DbDbgLevel = DbOOMs::DbFatalError+1;  
                    break;
                  case MSG::VERBOSE:  
                    DbOOMs::DbDbgLevel = DbOOMs::DbNoError;  
                    break;
                  case MSG::DEBUG:   
                    DbOOMs::DbDbgLevel = DbOOMs::DbNoError;  
                    break;
                  case MSG::INFO:     
                    DbOOMs::DbDbgLevel = DbOOMs::DbNoError;  
                    break;
                  case MSG::WARNING:  
                    DbOOMs::DbDbgLevel = DbOOMs::DbWarning;  
                    break;
                  case MSG::ERROR:    
                    DbOOMs::DbDbgLevel = DbOOMs::DbFatalError;  
                    break;
                  case MSG::FATAL:    
                    DbOOMs::DbDbgLevel = DbOOMs::DbFatalError;  
                    break;
                  default:   /* Let's take the default */
                    if ( debugLevel < MSG::NIL  ) debugLevel = MSG::VERBOSE;
                    if ( debugLevel > MSG::FATAL) debugLevel = MSG::FATAL;
                    break;
                  }
                  return StatusCode::SUCCESS;
                }
                else    {
                  log << MSG::ERROR << "Cannot open Federation " << m_fddbName << endreq;
                }
              }
              else    {
                log << MSG::ERROR << "Cannot open DB session " << endreq;
              }
            }
            else  {
              log << MSG::ERROR << "Cannot attach database implementation " << endreq;
            }
          }
          else    {
            log << MSG::ERROR << "Database type does not match sevice type " << endreq;
          }
        }
        else    {
            log << MSG::ERROR << "Database access cannot be initialized " << endreq;
        }
      }
      else  {
        log << MSG::ERROR << "Cannot load implementing factory " << endreq;
      }
      return StatusCode::FAILURE;
    }
  }
  return status;
}

/// Finalize the Db data persistency service
StatusCode DbCnvSvc::finalize()    {
  MsgStream log(messageService(), name());
  StatusCode status = ConversionSvc::finalize();
  m_output.close();
  m_federation.close();
  m_session.close();
  if ( m_manager )    {
    m_manager->shutdown(repSvcType());
    m_manager->release();
  }
  log << MSG::DEBUG << "Number of remaining DbObject instances:" << DbObject::numInstances() << endreq;
  return status;
}

/// Conect output stream (valid until overwritten)
StatusCode DbCnvSvc::connectOutput(const std::string& fname)   {
  if ( !m_output.isValid() || fname != m_output.name() )   {
    if ( m_federation.isValid() )   {
      if ( !(m_federation.openMode()&DbOOMs::CREATE) &&
           !(m_federation.openMode()&DbOOMs::UPDATE) )    {
        m_federation.close();
        if ( !m_federation.open(m_session, m_fddbName, DbOOMs::UPDATE) )   {
          return StatusCode::FAILURE;
        }
      }
      if ( !m_federation.existsDbase(fname) )    {
        m_output.open(m_federation, fname, DbOOMs::CREATE);
      }
      else    {
        m_output.open(m_federation, fname, DbOOMs::UPDATE);
      }
    }
    else    {
      dbHandle<DbDataBase> hdl(repSvcType());
      m_output = hdl;
    }
  }
  return m_output.isValid() ? StatusCode::SUCCESS : StatusCode::FAILURE;
}

/// Request an iterator over a container from the service
dbIter<DbObject>* DbCnvSvc::createIterator(const std::string& dbName, const std::string& cntName)   {
  unsigned char dbType = repSvcType();
  dbHandle<DbDataBase> dbH(dbType);
  if ( m_federation.isValid() )   {
    if ( dbH.open(m_federation, dbName, DbOOMs::READ) )   {
      dbHandle<DbContainer> cntH(dbType);
      if ( cntH.open(dbH, cntName) )   {
        dbIter<DbObject>* iter = new dbIter<DbObject>();
        iter->scan(cntH);
        iter->next();
        return iter;
      }
    }
  }
  return 0;
}

/// Address creation
IOpaqueAddress* DbCnvSvc::instantiate(const GenericLinkBase& link,
                                      const std::string& dbName, 
                                      const std::string& containerName,
                                      const std::string& /* objectName */ ) const    {
  dbHandle<DbDataBase> dbH(repSvcType());
  if ( dbH.open(m_federation, dbName) )   {
    dbHandle<DbContainer> cntH(repSvcType());
    if ( cntH.open(dbH, containerName) )   {
      if ( link.svcType() == repSvcType() )    {
        dbHandle<DbLink> linkH(link);
        DbAddress* pAdd = new DbAddress(cntH, linkH);
        return pAdd;
      }
    }
  }
  return 0;
}
