// -*- Mode: c -*-
#ifndef fsw_iface_H
#define fsw_iface_H
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
/** @class fsw_iface
    @author E. Charles

*/


#ifdef __cplusplus
extern "C" {
#endif


  // Load a CDM to memory
  //
  // returns a pointer to the loaded memory
  //
  // fileName is the full file name of the shared object we are loading
  // option is always zero for now
  // key is set the fmx key of the loaded file
  // schemaID is set to the ID of the schema of the loaded file (ie, which CDM);
  // versionID is set to the ID of the schema version of the loaded file
  // instanceID is set to the ID associated with this instance of the schema in the loaded file
  const void* fsw_get_data( const char* fileName, 
			    unsigned int option, 
			    unsigned int* key,
			    short* schemaID,
			    short* versionID,
			    short* instanceID );

#ifdef __cplusplus
}
#endif

#endif    /* CDM_PUBDEFS_H */
