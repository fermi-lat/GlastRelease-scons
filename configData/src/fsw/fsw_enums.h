// -*- Mode: c++ -*-
#ifndef fsw_enums_HH
#define fsw_enums_HH
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
/** @namespace configData::fsw_enums
    @author E. Charles

*/


namespace configData {

  namespace fsw_enums {

    // These are the FMX schema IDs
    enum { GFC_SCHEMA_ID = 0xe1,
	   MFC_SCHEMA_ID = 0xe3,
	   HFC_SCHEMA_ID = 0xe4,
	   DFC_SCHEMA_ID = 0xe5 };
  }

}

#endif
