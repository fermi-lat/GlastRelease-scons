#ifndef EVENTCNVSVC_H
#define EVENTCNVSVC_H 1

#include <map>

#include "GaudiKernel/ConversionSvc.h"
#include "IEventCnvSvc.h"

template <class TYPE> class SvcFactory;

/** @class EventCnvSvc
 * @brief GLAST Event Conversion Service
 * Based on SICb service written by Markus Frank.
 *
 * $Header$
 */

class EventCnvSvc  : virtual public ConversionSvc, virtual public IEventCnvSvc	{

  friend class SvcFactory<EventCnvSvc>;

  /// Map with leaf entries
  typedef std::map<std::string, Leaf*> LeafMap;
  LeafMap m_leaves;


public:

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  virtual StatusCode updateServiceState(IOpaqueAddress* pAddress);

  /// Override inherited queryInterface due to enhanced interface
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);

  /// Store path on TDS for a particular converter
  virtual StatusCode declareObject(const IEventCnvSvc::Leaf& leaf);

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

  EventCnvSvc(const std::string& name, ISvcLocator* svc);

  virtual ~EventCnvSvc() { };

};


#endif // EVENTCNVSVC_H
