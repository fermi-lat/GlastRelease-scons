//====================================================================
//	DbExportConverter.h
//--------------------------------------------------------------------
//
//	Package    : DbCnv/DbSvc
//
//  Description: Definition of Dll exportable converters
//
//	Author     : M.Frank
//  Created    : 13/1/99
//	Changes    : 
//
//====================================================================
#ifndef DBCNV_DBEXPORTCONVERTER_H
#define DBCNV_DBEXPORTCONVERTER_H 1

// Framework include files
#include "DbCnv/DbContainerConverter.h"
#include "DbCnv/DbUserCnvFactory.h"

extern const CLID& CLID_Any;
extern const CLID& CLID_ObjectVector;
extern const CLID& CLID_ObjectList;

/**

*/
class DbExportConverter : public DbGenericConverter   {
  /// Transient object creator function
  IDbUserCnvFactory::Creator       m_creatorFunc;
public:
  class Object : public DataObject   {
  public:
    static const CLID& classID()    {
      return CLID_Any;
    }
  };
  /// Standard Constructor
  DbExportConverter(IDbUserCnvFactory::Creator fun, const CLID& clid, IOODataBase* db, ISvcLocator* loc);
  /// Standard Destructor
  virtual ~DbExportConverter();
  /// Create transient representation of the requested type
  virtual DataObject* createTransient();
};

/**

*/
class DbExportContCnv : public DbContainerConverter   {
  /// Transient object creator function
  IDbUserCnvFactory::Creator       m_creatorFunc;
public:
  class Vector : public DataObject   {
  public:
    static const CLID& classID()    {
      static CLID clid = CLID_Any|CLID_ObjectVector;
      return clid;
    }
  };
  class List : public DataObject   {
  public:
    static const CLID& classID()    {
      static CLID clid = CLID_Any|CLID_ObjectList;
      return clid;
    }
  };
  /// Standard Constructor
  DbExportContCnv(IDbUserCnvFactory::Creator fun, const CLID& clid, IOODataBase* db, ISvcLocator* loc);
  /// Standard Destructor
  virtual ~DbExportContCnv();
  /// Create transient representation of the requested type
  virtual DataObject* createTransient();
};
#endif // DBCNV_DBEXPORTCONVERTER_H
