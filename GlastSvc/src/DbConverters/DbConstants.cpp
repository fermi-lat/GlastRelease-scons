//====================================================================
//	Constants definition file
//--------------------------------------------------------------------
//
//	Package    : System ( The LHCb Offline System)
//
//  Description: Constants definition
//
//	Author     : M.Frank
//  History    :
// +---------+----------------------------------------------+---------
// |    Date |                 Comment                      | Who     
// +---------+----------------------------------------------+---------
// | 11/07/99| Initial version                              | MF
// +---------+----------------------------------------------+---------
//====================================================================
#define DBCNV_CONSTANTS_CPP
// $Header$

#include "DbCnv/DbOOMs.h"
unsigned char MEMORY_StorageType = 0x5;

const DbAccessMode  DbOOMs::READ          = 1;
const DbAccessMode  DbOOMs::UPDATE        = 2;
const DbAccessMode  DbOOMs::CREATE        = 4;
const DbAccessMode  DbOOMs::NOT_OPEN      = 8;

const DbResult      DbOOMs::DbSuccess     = 1;
const DbResult      DbOOMs::DbError       = 0;
const DbResult      DbOOMs::DbNoError     = DbOOMs::DbSuccess;
const DbResult      DbOOMs::DbWarning     = 2;
const DbResult      DbOOMs::DbUserError   = 4;
const DbResult      DbOOMs::DbSystemError = 6;
const DbResult      DbOOMs::DbFatalError  = 8;
      DbResult      DbOOMs::DbDbgLevel    = DbOOMs::DbFatalError;
