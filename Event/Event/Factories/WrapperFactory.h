// $Id$
// 
//  Original author: Sawyer Gillespie
//                   hgillesp@u.washington.edu
//

#ifndef _H_GlastEvent_WrapperFactory_
#define _H_GlastEvent_WrapperFactory_

// includes
#include "Gaudi/Kernel/ObjectFactory.h"

/** Templated factory to create concrete instances of a given contained object
*/
template <class TYPE>
class ContainedWrapperObjectFactory 
: virtual public IContainedObjectFactory, 
  virtual public ObjectFactoryBase
{
public:
  /// Standard Constructor
  ContainedObjectFactory()   {
    setTypeName(typeid(*this));
  }
  /// Standard destructor
  virtual ~ContainedObjectFactory()    {
  }
  /// Create an instance of a generic ContainedObject object
  virtual ContainedObject* instantiate () const {
    ContainedObject* result = new TYPE();
    return result;
  }
  /// Access to the service type of the converter
  virtual const CLID& clID()   const  {
    return TYPE::classID();
  }
  /// IInterface implementation: queryInterface
  virtual StatusCode queryInterface(const IID& riid, void** ppIF) {
	  if ( IID_IContainedObjectFactory == riid ) {
      *ppIF = (IContainedObjectFactory*) this;
    }
	  else {
      return ObjectFactoryBase::queryInterface(riid, ppIF);
    }
    addRef();
 	  return StatusCode::SUCCESS;
  }
};

/** Templated factory to create concrete instances of a given contained object
*/
template <class TYPE>
class DataObjectFactory 
: virtual public IDataObjectFactory, 
  virtual public ObjectFactoryBase 
{
public:
  /// Standard Constructor
  DataObjectFactory()   {
    setTypeName(typeid(*this));
  }
  /// Standard destructor
  virtual ~DataObjectFactory()    {
  }
  /// Create an instance of a generic ContainedObject object
  virtual DataObject* instantiate () const {
    DataObject* result = new TYPE();
    return result;
  }
  /// Access to the service type of the converter
  virtual const CLID& clID()   const  {
    return TYPE::classID();
  }
  /// IInterface implementation: queryInterface
  virtual StatusCode queryInterface(const IID& riid, void** ppIF) {
	  if ( IID_IDataObjectFactory == riid ) {
      *ppIF = (IDataObjectFactory*) this;
    }
	  else {
      return ObjectFactoryBase::queryInterface(riid, ppIF);
    }
    addRef();
 	  return StatusCode::SUCCESS;
  }
};

//
// ClassName:   ObjectFactory
//
// Description: Templated factory to create concrete instances of a given Gaudi class
//
template <class T>
class ObjectFactory : virtual public ObjectFactoryBase {
public:
  // Default Constructor
  ObjectFactory() {
    setTypeName(typeid(T));
  }
  // Default destructor
  virtual ~ObjectFactory() { 
  }
  // Instantiate an instance of a Gaudi class
	virtual IInterface* instantiate( IInterface *parent ) const {
    T* obj = new T( parent ); 
    return obj;
	}
  /// IInterface implementation: queryInterface
  virtual StatusCode queryInterface(const IID& riid, void** ppIF) {
	  if ( IID_IObjectFactory == riid ) {
      *ppIF = (IFactory*) this;
    }
	  else {
      return ObjectFactoryBase::queryInterface(riid, ppIF);
    }
    addRef();
 	  return StatusCode::SUCCESS;
  }
};

#define _ImplementContainedObjectFactory( Transient )                        \
static ContainedObjectFactory< Transient > s_##Transient##Factory;           \
const IFactory& Transient##Factory = s_##Transient##Factory;

#define _ImplementDataObjectFactory( Transient )                             \
static DataObjectFactory< Transient > s_##Transient##Factory;                \
const IFactory& Transient##Factory = s_##Transient##Factory;

#define _ImplementObjectFactory( Transient )                                 \
static ObjectFactory< Transient > s_##Transient##Factory;                    \
const IFactory& Transient##Factory = s_##Transient##Factory;

#define DLL_DECL_OBJECTFACTORY(x)    extern const IFactory& x##Factory; x##Factory.addRef();

#endif // _H_GlastEvent_WrapperFactory_