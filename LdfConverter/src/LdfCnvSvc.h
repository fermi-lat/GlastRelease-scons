#ifndef LdfCnvSvc_H
#define LdfCnvSvc_H 1

#include <map>

#include "GaudiKernel/ConversionSvc.h"
#include "ILdfCnvSvc.h"

template <class TYPE> class SvcFactory;

class IOpaqueAddress;

/** @class LdfCnvSvc
 * @brief GLAST EBF Conversion Service which coordinates all of our EBF converters.
 *
 * When a component requests an item not yet available on the TDS, 
 * the LdfCnvSvc is called to find the appropriiate converter to gain
 * access to the data and put it on the TDS.
 * Based on SICb service written by Markus Frank.
 *
 * $Header$
 */

class LdfCnvSvc  : virtual public ConversionSvc, virtual public ILdfCnvSvc	{

  friend class SvcFactory<LdfCnvSvc>;

  /// Map with leaf entries containing TDS paths
  typedef std::map<std::string, Leaf*> LeafMap;
  LeafMap m_leaves;

public:

  virtual StatusCode initialize();

  virtual StatusCode finalize();

// HMK these 2 are not in the CalibXmlCnvSvc
  virtual StatusCode updateServiceState(IOpaqueAddress* pAddress);
  /// Associates a path on TDS with a particular converter
  virtual StatusCode declareObject(const ILdfCnvSvc::Leaf& leaf);

  /// Override inherited queryInterface due to enhanced interface
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);

  /** IAddressCreator implementation: Address creation.
    Create an address using the link infotmation together with
    the triple (database name/container name/object name).

    @param refLink        Reference to abstract link information
    @param dbName         Database name
    @param containerName  Object container name
    @param refpAddress    Opaque address information to retrieve object
    @return               StatusCode indicating SUCCESS or failure
  */
  virtual StatusCode createAddress( unsigned char svc_type,
                                    const CLID& clid,
                                    const std::string* par, 
                                    const unsigned long* ip,
                                    IOpaqueAddress*& refpAddress);
protected:

  LdfCnvSvc(const std::string& name, ISvcLocator* svc);

  virtual ~LdfCnvSvc() { };

};


#endif // LdfCnvSvc_H
