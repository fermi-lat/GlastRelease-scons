// $Header$
#ifndef GLASTCNV_GLASTADDRESS_H
#define GLASTCNV_GLASTADDRESS_H 1

// Include files
#include "GaudiKernel/GenericAddress.h"
//typedef GenericAddress GlastIRFAddress;

//------------------------------------------------------------------------------
//
// ClassName :    GlastIRFAddress
//  
// Description :  Definition of a Glast relative address
//
// Author :       Markus Frank
//
//------------------------------------------------------------------------------
class GlastIRFAddress : public GenericAddress   {
  /// User parameter used by converter e.g. partition number
  long m_userParameter;
public:
  /// Standard Constructor
    GlastIRFAddress(const CLID& cl, const std::string& fname, int recid, const std::string& containerName);
  /// Constructor for non-initialized addresses
  GlastIRFAddress(const GenericLink& link);
  /// Constructor for non-initialized addresses
  GlastIRFAddress(const GenericLink* link);
  /// Standard destructor
  virtual ~GlastIRFAddress()    {
  }
  /// Set user parameter used by converter e.g. partition number
  void setUserParameter(long par)   {
    m_userParameter = par;
  }
  /// Access user parameter used by converter e.g. partition number
  long userParameter()    const   {
    return m_userParameter;
  }
};
#endif  // GLASTCNV_GLASTADDRESS_H
