//---------------------------------------------------------------------------
// File and Version Information:
//      $Id$
//
// Description:
//      Base class for converting FSW headers to XML
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

// fsw headers
#include "fsw/CDM/CDM_pubdefs.h"

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

