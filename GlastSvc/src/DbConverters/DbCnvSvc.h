//====================================================================
//	DbCnvSvc definition
//--------------------------------------------------------------------
//
//	Package    : System ( The LHCb Offline System)
//
//  Description: 
//
//	Author     : M.Frank
//====================================================================
#ifndef DBCNV_DBCNVSVC_H
#define DBCNV_DBCNVSVC_H
// $Header$

// Framework include files
#include "Gaudi/Kernel/AddrFactory.h"
#include "Gaudi/PersistencySvc/ConversionSvc.h"

// Db include files
#include "DbCnv/IDataBaseMgr.h"
#include "DbCnv/DbSession.h"
#include "DbCnv/DbDataBase.h"
#include "DbCnv/DbFederation.h"

// Forward declarations
template <class ConcreteService> class SvcFactory;
class DbDataBaseMgr;

/** DbCnvSvc class implementation definition.

    <P> System:  The LHCb Offline System
    <P> Package: DbCnv

     Dependencies:
    <UL>
    <LI> DbCnvSvc interface definition:   DbCnv/IDataBaseMgr.h
    <LI> ConversionSvc definition:          Gaudi/PersistencySvc/ConversionSvc.h
    </UL>

    <P> History    :
    <PRE>
     +---------+----------------------------------------------+---------
     |    Date |                 Comment                      | Who     
     +---------+----------------------------------------------+---------
     | 29/10/98| Initial version                              | MF
     +---------+----------------------------------------------+---------
    </PRE>

   @author Markus Frank
   @version 1.0
*/
class DbCnvSvc : public ConversionSvc, 
                 virtual public IDataBaseMgr,
                 virtual public IAddrFactory
{
  /// Friend: service factory
  friend SvcFactory<DbCnvSvc>;
protected:
  /// This is the one and only one database manager
  DbDataBaseMgr*          m_manager;
  /// Reference to Db main object
  dbHandle<DbSession>     m_session;
  /// Reference to handle of the federation
  dbHandle<DbFederation>  m_federation;
  /// Reference to output file
  dbHandle<DbDataBase>    m_output;
  /// FDBD name
  std::string             m_fddbName;
  /// Name of the OODataBase implementation
  std::string             m_implementation;
  /// Standard constructor
  DbCnvSvc(const std::string& name, ISvcLocator* svc);
  /// Standard destructor
  virtual ~DbCnvSvc();

public:
  /// ConversionSvc overload: Query interface
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);
  /// ConversionSvc overload: initialize Db service
  virtual StatusCode initialize();
  /// ConversionSvc overload: Finalize Db service
  virtual StatusCode finalize();
  /// ConversionSvc overload: Locate converter factory
  virtual const ICnvFactory* findCnvFactory(const CLID& wanted);
  /// ConversionSvc overload: Create new Converter using factory
  virtual IConverter* createConverter(const ICnvFactory& fac);

  /// Give access to the default data base handle for output
  virtual dbHandle<DbDataBase>& defaultDbHdl()   {
    return m_output;
  }
  /** Connect to output file.
      @param fname        String with file name.
      @return             Status code indicating success or failure.
  */
  virtual StatusCode connectOutput(const std::string& fname);

  /** IAddrFactory Interface: instantiate a concrete address
      @param classId      Class id of the object the address refers to
      @param fileName     String with file name.
      @param contName     String with container name.
      @param seqId        Id with sequence number of the object the address refers to
  */
  virtual IOpaqueAddress* instantiate(const CLID& /* classId */ , 
                                      const std::string& /* fileName */, 
                                      const std::string& /* contName */, 
                                      int /* seqId */  ) const    {
    return 0;
  }
  /** Address creation.
    Create an address using the link infotmation together with
    the triple (database name/container name/object name).

    @param refLink        Reference to abstract link information
    @param dbName         Database name
    @param containerName  Object container name
    @param objectName     Object name (if at all relevant)
    @return Opaque address information to retrieve object
  */
  virtual IOpaqueAddress* instantiate(const GenericLinkBase& link,
                                      const std::string& dbName, 
                                      const std::string& containerName,
                                      const std::string& objectName) const;
  /** Request a specific implementation of an OO database
      @param      accessor    Reference to pointer of the requested implementation interface
      @return     Status code indicating success or failure.
  */
  virtual StatusCode requestImplementation(IOODataBase*& accessor);
  /** Request an iterator over a container from the service
      @param      dbName      String containing name of the database
      @param      cntName     String containing the name of the container object
      @return     Reference to the iterator object
  */
  virtual dbIter<DbObject>* createIterator(const std::string& dbName, const std::string& cntName);
};

#endif  // DBCNV_DBCNVSVC_H
