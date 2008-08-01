//$Header$
#ifndef IConfigSvc_h
#define IConfigSvc_h 1

#include "GaudiKernel/IInterface.h"

// External constants
#include "GaudiKernel/ClassID.h"
#include "enums/Lsf.h"

// Forward declarations


class TrgConfig;
class FswEfcSampler;


static const InterfaceID IID_IConfigSvc ("IConfigSvc", 1, 0);

/** @class IConfigSvc 

    Abstract interface of a service for access to Configuration information.
    See also data class definitions in configData

    @author Eric Charles, from Martin Kocian's ITrgConfigSvc
*/

class IConfigSvc : virtual public IInterface   {

public:
  
  // Re-implemented from IInterface
  static const InterfaceID& interfaceID() { return IID_IConfigSvc; }
  
  /// get the MOOT Key
  virtual unsigned getMootKey() const = 0;

  /// get the GEM configuration object
  virtual const TrgConfig* getTrgConfig() const = 0;
  
  /// get the information about the prescalers
  virtual const FswEfcSampler* getFSWPrescalerInfo( enums::Lsf::Mode mode, unsigned handlerId, 
						    unsigned int& fmxKey ) const = 0;
  
};

#endif

