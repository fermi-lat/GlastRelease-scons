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
