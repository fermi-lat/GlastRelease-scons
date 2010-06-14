// $Header$
#define EVENTCNVSVC_CPP


// Include files
#include <iostream>
#include "Gaudi/Kernel/SvcFactory.h"
#include "Gaudi/Kernel/CnvFactory.h"
#include "Gaudi/Kernel/AddrFactory.h"
#include "Gaudi/Kernel/HashTable.h"
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Interfaces/ICnvManager.h"
#include "Gaudi/Interfaces/ISvcLocator.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/RegistryEntry.h"

#include "Gaudi/Kernel/AddrFactory.h"

#include "Address.h"
#include "EventCnvSvc.h"


// RCS Id for identification of object version
static const char* rcsid = "$Id$";


//------------------------------------------------------------------------------
//
// Implementation of class :  EventCnvSvc
//
// modified SicbEventCnvSvc
//
//------------------------------------------------------------------------------


// Instantiation of a static factory class used by clients to create
// instances of this service
static SvcFactory<EventCnvSvc> s_factory;
const ISvcFactory& EventCnvSvcFactory = s_factory;


// External factories

extern const IAddrFactory& DummyAddressFactory;
extern const IID& IID_ISicbEventCnvSvc;
extern const ICnvFactory& EventCnvFactory;


/// Initialize the service.
StatusCode EventCnvSvc::initialize()     {
  StatusCode status = ConversionSvc::initialize();
  if ( status.isSuccess() )   {
    ISvcLocator* svclocator = serviceLocator();
    //======================================================================
    // Add known converters to the service: 
    //======================================================================
    IDataProviderSvc *pIDP = 0;
    // Set event data service
    status = serviceLocator()->getService("EventDataSvc", IID_IDataProviderSvc, (IInterface*&)pIDP);
    if ( status.isSuccess() )   {
      status = setStore ( pIDP );
    }
    else    {
      return status;
    }

    ICnvManager::CnvIterator i, stop;
    for ( i = cnvManager()->cnvBegin(), stop = cnvManager()->cnvEnd(); i != stop; i++ )   {
      if ( repSvcType() == (*i)->repSvcType() )   {
        StatusCode iret = addConverter( *(*i) );
        if ( !iret.isSuccess() )   {
          status = iret;
        }
      }
    }
    
    // Now we have to configure the map of leaves
    for (LeafMap::iterator k = m_leaves.begin(); k != m_leaves.end(); k++ )   {
      std::string path = (*k).first;
      for (LeafMap::iterator j = m_leaves.begin(); j != m_leaves.end(); j++ )   {
        std::string pp = (*j).first.substr(0, (*j).first.rfind("/"));
        if ( path == pp && path != (*j).first )   {
          (*k).second->push_back((*j).second);
        }
      }
    }

    // HMA commented this out....replaced with the above code from SicbEventCnvSvc
    /*
	IConverter* cnv = 0;
    MsgStream log(messageService(), name());
    log << MSG::DEBUG << "initializing: add event factory" << endreq;

	// Event
	cnv = EventCnvFactory.instantiate(svclocator);
    status = addConverter( cnv );
    if ( status.isFailure() ) return status;

    log << MSG::DEBUG << "initializing: event factory added OK" << endreq;
*/
  }

  return status;
}

/// Finalize the service.
StatusCode EventCnvSvc::finalize()     {
  StatusCode status = ConversionSvc::finalize();
  for (LeafMap::iterator k = m_leaves.begin(); k != m_leaves.end(); k++ )   {
    delete (*k).second;
  }
  m_leaves.erase(m_leaves.begin(), m_leaves.end());
  return status;
}

StatusCode EventCnvSvc::declareObject(const IEventCnvSvc::Leaf& leaf)   {
  Leaf* ll = new Leaf(leaf);
  std::pair<LeafMap::iterator, bool> p = m_leaves.insert(LeafMap::value_type( leaf.path, ll) );
  if( p.second )    {
    return StatusCode::SUCCESS;
  }
  delete ll;
  return StatusCode::FAILURE;
}

/// Update state of the service
StatusCode EventCnvSvc::updateServiceState(IOpaqueAddress* pAddress)    {

    // not sure about the use of recid or bank...

  MsgStream log(messageService(), name());
  log << MSG::DEBUG << "inside EventCnvSvc::updateServiceState()" << endreq;

  // convert to an address of an Event
	Address* addr = dynamic_cast<Address*>(pAddress); 
    if ( 0 != addr ) {
      const std::string& new_fid = addr->dbName();
      int new_recid = addr->sicbInfo()->entryID(); // HMA???

      IDataDirectory* dir = addr->directory();
      if ( 0 != dir )   {
        RegistryEntry* ent = dynamic_cast<RegistryEntry*>(dir);
        if ( 0 != ent && 0 == ent->object() )   {
          std::string path = dir->fullpath();
          LeafMap::iterator itm = m_leaves.find(path);
          if ( itm != m_leaves.end() )    {
            IEventCnvSvc::Leaf* leaf = (*itm).second;
            if ( 0 != leaf )    {
              //addr->setUserParameter(leaf->userParameter);
              for ( Leaf::iterator il = leaf->begin(); il != leaf->end(); il++ )   {
                RegistryEntry* newEntry = new RegistryEntry((*il)->path, 0);
                newEntry->setAddress(new Address((*il)->clid, new_fid, new_recid, (*il)->bank));
                ent->add( newEntry );
              }
            }
          }
        }
      }
        return StatusCode::SUCCESS;
    }
    return StatusCode::FAILURE;
}

/// Query interface
StatusCode EventCnvSvc::queryInterface(const IID& riid, void** ppvInterface)  {
  if ( IID_ISicbEventCnvSvc == riid )  {
    *ppvInterface = (IEventCnvSvc*)this;
  }
  else  {
    // Interface is not directly availible: try out a base class
    return ConversionSvc::queryInterface(riid, ppvInterface);
  }
  addRef();
  return StatusCode::SUCCESS;
}


/// Standard Constructor
EventCnvSvc::EventCnvSvc(const std::string& name, ISvcLocator* svc)
  : ConversionSvc(name, svc, SICB_StorageType)              
{
  // Set address factory for my type
  setAddressFactory(&DummyAddressFactory);
}


/// Standard Destructor
EventCnvSvc::~EventCnvSvc()                         
{
}
