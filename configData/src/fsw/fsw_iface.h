// -*- Mode: c -*-
#ifndef fsw_iface_H
#define fsw_iface_H
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
/** @class fsw_iface
    @author E. Charles

*/


#ifdef __cplusplus
extern "C" {
#endif


  /**
   * @brief Load a CDM to memory and set some values for bookkeeping
   *
   * @param fileName is the full file name of the shared object we are loading
   * @param  option is always zero for now
   * @param key is set the fmx key of the loaded file
   * @param schemaID is set to the ID of the schema of the loaded file (ie, which CDM);
   * @param versionID is set to the ID of the schema version of the loaded file
   * @param instanceID is set to the ID associated with this instance of the schema in the loaded file
   * @return a point to the memory that contains the CDM
   **/
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
