//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Function to load CDM and set bookkeeping values
//
// Environment:
//      Software developed for GLAST.
//
// Author List:
//      Eric Charles
//
// Copyright Information:
//      Copyright (C) 2008      Stanford Linear Accelerator Center
//
//---------------------------------------------------------------------------


// headers
#include "./fsw_iface.h"


// This is a c file, not a c++ file, b/c CDM_pubdefs uses a variable called "new"

// fsw headers
#include "fsw/CDM/CDM_pubdefs.h"

// Just the wrapped functions
const void* fsw_get_data( const char* fileName, 
			  unsigned int option,
			  unsigned int* key,
			  short* schemaID,
			  short* versionID,
			  short* instanceID ) {
  CDM_Database* cdb = CDM_loadDatabase(fileName,option);
  if ( cdb == 0 ) return 0;
  *key = CDM_getKey( cdb );
  *schemaID = CDM_getSchema( cdb );
  *versionID = CDM_getVersion( cdb );
  *instanceID = CDM_getInstance( cdb );
  return CDM_getData(cdb);
}

