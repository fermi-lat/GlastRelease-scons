//====================================================================
//	DbExportConverter.cpp
//--------------------------------------------------------------------
//
//	Package    : Gaudi/System
//
//  Description: Implementation of Dll exportable converters
//
//	Author     : M.Frank
//  Created    : 13/1/99
//	Changes    : 
//
//====================================================================
#define DBCNV_DBEXPORTCONVERTER_CPP

#include "DbCnv/DbObject.h"
#include "src/Converters/DbExportConverter.h"

static const DbUserCnvFactory <DbExportConverter::Object, DbExportConverter> s_EportConverterFactory;
const ICnvFactory& DbOOMsExportConverterFactory = s_EportConverterFactory;

static const DbUserCnvFactory <DbExportContCnv::Vector, DbExportContCnv> s_EportVectorCnvFactory;
const ICnvFactory& DbOOMsExportVectorCnvFactory = s_EportVectorCnvFactory;

static const DbUserCnvFactory <DbExportContCnv::List, DbExportContCnv> s_EportListCnvFactory;
const ICnvFactory& DbOOMsExportListCnvFactory = s_EportListCnvFactory;

/// Standard Constructor
DbExportConverter::DbExportConverter (IDbUserCnvFactory::Creator fun, 
                                      const CLID& clid, 
                                      IOODataBase* db, 
                                      ISvcLocator* loc)
: DbGenericConverter(clid, db, loc), 
  m_creatorFunc(fun)
{
}

/// Standard Destructor
DbExportConverter::~DbExportConverter()    {
}

/// Create transient representation of the requested type
DataObject* DbExportConverter::createTransient ()    {
  DataObject* pObject = 0;
  if ( 0 != m_creatorFunc )   {
    pObject = (*m_creatorFunc)();
  }
  return pObject;
}

/// Standard Constructor
DbExportContCnv::DbExportContCnv (IDbUserCnvFactory::Creator fun, 
                                      const CLID& clid, 
                                      IOODataBase* db, 
                                      ISvcLocator* loc)
: DbContainerConverter(clid, db, loc), 
  m_creatorFunc(fun)
{
}

/// Standard Destructor
DbExportContCnv::~DbExportContCnv()    {
}

/// Create transient representation of the requested type
DataObject* DbExportContCnv::createTransient ()    {
  DataObject* pObject = 0;
  if ( 0 != m_creatorFunc )   {
    pObject = (*m_creatorFunc)();
  }
  return pObject;
}

