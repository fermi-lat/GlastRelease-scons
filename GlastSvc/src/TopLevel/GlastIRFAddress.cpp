// $Header$
#define GLASTCNV_GLASTADDRESS_CPP 

// Include files
#include "Gaudi/Kernel/AddrFactory.h"
#include "GlastSvc/TopLevel/GlastIRFAddress.h"

// RCS Id for identification of object version
static const char* rcsid = "$Id$";


//------------------------------------------------------------------------------
//
// Implementation of class :  GlastIRFAddress
//
// Author :                   Markus Frank
//
//------------------------------------------------------------------------------
unsigned char GLASTIRF_StorageType = 200;   // don't quite know about this (HSG)

// Address Factory definition
static const AddrFactory<GlastIRFAddress> s_Factory(GLASTIRF_StorageType);
const IAddrFactory& GlastIRFAddressFactory = s_Factory;

/// Standard constructor
GlastIRFAddress::GlastIRFAddress(const CLID& cl, const std::string& fid, int recid, const std::string& containerName)
: GenericAddress(GLASTIRF_StorageType, cl)
{
  setDbName(fid);
  setContainerName(containerName);
  setObjectName("GlastIRFData");
  sicbInfo()->setEntryID(recid);
  setUserParameter(0);
}

/// Standard Constructor for non-initialized addresses
GlastIRFAddress::GlastIRFAddress(const GenericLink& link)
: GenericAddress(link)
{
  setUserParameter(0);
}

/// Standard Constructor for non-initialized addresses
GlastIRFAddress::GlastIRFAddress(const GenericLink* link)
: GenericAddress(*link)
{
  setUserParameter(0);
}

