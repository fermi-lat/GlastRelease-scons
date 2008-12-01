// File and Version Information:
//      $Header$
//
// Description:
//      EventCnvSvc is the GLAST converter service.
//
// Author(s):

#define EVENTCNVSVC_CPP

#include <iostream>
#include "GaudiKernel/ConversionSvc.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/ICnvManager.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/RegistryEntry.h"

#include "Address.h"

#include "GlastSvc/EventSelector/IGlastCnv.h"

#include <map>

static const InterfaceID IID_EventCnvSvc("EventCnvSvc", 1, 0);

/** @class EventCnvSvc
 * @brief GLAST Event Conversion Service which coordinates all of our converters.
 *
 * When a component requests an item not yet available on the TDS, 
 * the EventCnvSvc is called to find the appropriiate converter to gain
 * access to the data and put it on the TDS.
 * Based on SICb service written by Markus Frank.
 *
 * $Header$
 */

class EventCnvSvc  : virtual public ConversionSvc	
{

  friend class SvcFactory<EventCnvSvc>;

public:

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  virtual StatusCode updateServiceState(IOpaqueAddress* pAddress);

  /// Override inherited queryInterface due to enhanced interface
  virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);

  /** IAddressCreator implementation: Address creation.
    Create an address using the link infotmation together with
    the triple (database name/container name/object name).

    @param refLink        Reference to abstract link information
    @param dbName         Database name
    @param containerName  Object container name
    @param refpAddress    Opaque address information to retrieve object
    @return               StatusCode indicating SUCCESS or failure
  */
  virtual StatusCode createAddress( long  svc_type,
                                    const CLID& clid,
                                    const std::string* par, 
                                    const unsigned long* ip,
                                    IOpaqueAddress*& refpAddress);
 // virtual StatusCode createAddress( unsigned char svc_type,
 //                                   const CLID& clid,
 //                                   const std::string* par, 
 //                                   const unsigned long* ip,
 //                                   IOpaqueAddress*& refpAddress);
protected:

  EventCnvSvc(const std::string& name, ISvcLocator* svc);

  virtual ~EventCnvSvc() { };

private:
    typedef std::map<std::string, IGlastCnv*>                PathToCnvMap;
    typedef std::map<std::string, std::vector<std::string> > SubPathMap;

    PathToCnvMap m_pathToCnvMap;
    SubPathMap   m_subPathMap;
};

#include "GlastSvc/EventSelector/IGlastCnv.h"

//static const InterfaceID IID_IBaseCnv(902, 1 , 0); 
// RCS Id for identification of object version
static const char* rcsid = "$Id$";


// Instantiation of a static factory class used by clients to create
// instances of this service
static const SvcFactory<EventCnvSvc> s_EventCnvSvcFactory;
const ISvcFactory& EventCnvSvcFactory = s_EventCnvSvcFactory;

EventCnvSvc::EventCnvSvc(const std::string& name, ISvcLocator* svc)
: ConversionSvc(name, svc, SICB_StorageType)              
{
    m_pathToCnvMap.clear();
    m_subPathMap.clear();
}

StatusCode EventCnvSvc::initialize()     {
    // Purpose and Method:  Setup GLAST's Event Converter Service.
    //   Associate EventCnvSvc with the EventDataSvc
    //   Associate the list of known converters with this EventCnvSvc.
    MsgStream log(msgSvc(), name());

    StatusCode status = ConversionSvc::initialize();

    if ( status.isSuccess() )   
    {
        ISvcLocator* svclocator = serviceLocator();
        IDataProviderSvc *pIDP = 0;
        // Set event data service
        status = service("EventDataSvc", pIDP, true);
        if ( status.isSuccess() )   
        {
            status = setDataProvider ( pIDP );
        }
        else    
        {
            log << MSG::ERROR << "Failed to set data provider service" << endreq;
            return status;
        }
        
        // Add known converters to the service: 
        ICnvManager::CnvIterator i, stop;
        for ( i = cnvManager()->cnvBegin(), stop = cnvManager()->cnvEnd(); i != stop; i++ )   
        {
            // Make sure the converters are of our "type" (SICB_StorageType)
            if ( repSvcType() == (*i)->repSvcType() )   
            {
                // Add the converter
                StatusCode iret = addConverter( (*i)->objType() );  

                // Was there a problem?
                if ( iret.isFailure() )   
                {
                    log << MSG::ERROR << "Unable to add converter! " << endreq;
                    return iret;
                }

                // Caste back to our conversion type
                IGlastCnv* glastCnv = dynamic_cast<IGlastCnv*>(converter((*i)->objType()));

                // Store the pointer to the converter in our map
                m_pathToCnvMap[glastCnv->getPath()] = glastCnv;
            }
        }

        // Loop back through the map to build the list of immediate daughter paths
        for(PathToCnvMap::iterator mapIter = m_pathToCnvMap.begin(); mapIter != m_pathToCnvMap.end(); mapIter++)
        {
            bool               orphan = true;
            const std::string& path   = mapIter->first;

            // Another loop through to find the matches
            for(PathToCnvMap::iterator inIter = m_pathToCnvMap.begin(); inIter != m_pathToCnvMap.end(); inIter++)
            {
                // Skip if our self
                if (path == inIter->first) continue;

                // Find the parent path of this path
                std::string prntPath = inIter->first.substr(0, inIter->first.rfind("/"));

                // If a match then add this as a daughter of the path
                if (path == prntPath)
                {
                    m_subPathMap[path].push_back(inIter->first);
                    orphan = false;
                }
            } 

            // An orphan?
            if (orphan)
            {
                int j = 0;
            }
        }
    }
    
    return status;
}

StatusCode EventCnvSvc::finalize()     {
    StatusCode status = ConversionSvc::finalize();

    return status;
}

/// Update state of the service
StatusCode EventCnvSvc::updateServiceState(IOpaqueAddress* pAddress)    
{    
    // not sure about the use of recid or bank...
    
    MsgStream log(msgSvc(), name());
    StatusCode status = INVALID_ADDRESS;
    IRegistry* ent    = pAddress->registry();

    // Need a registry entry to proceed
    if ( 0 != ent )   
    {
        SmartIF<IDataManagerSvc> iaddrReg(IID_IDataManagerSvc, dataProvider());

        status = StatusCode::SUCCESS;

        // Path we are updating from
        std::string path = ent->identifier();

        // Find the list of possible daughter paths
        SubPathMap::iterator pathIter = m_subPathMap.find(path);

        // Loop over the list of daughter paths
        for(std::vector<std::string>::iterator subIter = pathIter->second.begin(); subIter != pathIter->second.end(); subIter++)
        {
            // Now look up the pointer to the converter
            PathToCnvMap::iterator cnvIter  = m_pathToCnvMap.find(*subIter);
            IGlastCnv*             glastCnv = cnvIter->second;

            if (glastCnv->getPath() == path) continue;

            if (glastCnv)
            {
                IOpaqueAddress* newAddr = 0;
                unsigned long ipars[2] = {0, 0}; //{(*il)->userParameter, new_rid};
                const std::string spars[2] = {"", ""}; //{par[0], (*il)->bank};
                StatusCode ir =
            
                addressCreator()->createAddress(SICB_StorageType, 
                                                glastCnv->objType(), 
                                                spars, 
                                                ipars,
                                                newAddr);
                if ( ir.isSuccess() )   
                {
                    ir = iaddrReg->registerAddress(glastCnv->getPath(), newAddr);
                    if ( !ir.isSuccess() )    
                    {
                        newAddr->release();
                        status = ir;
                    }
                }
            }
        }
    }
    return status;
}

StatusCode EventCnvSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
    if ( IID_EventCnvSvc == riid )  {
        *ppvInterface = (EventCnvSvc*)this;
    }
    else  {
        // Interface is not directly availible: try out a base class
        return ConversionSvc::queryInterface(riid, ppvInterface);
    }
    addRef();
    return StatusCode::SUCCESS;
}

StatusCode EventCnvSvc::createAddress( long svc_type,
                                       const CLID& clid,
                                       const std::string* par, 
                                       const unsigned long* ip,
                                       IOpaqueAddress*& refpAddress)
{
    if ( svc_type == repSvcType() )   
    {
        Address* pAdd = new Address(svc_type, clid, par[0]);
        if ( pAdd )
        {
            refpAddress = pAdd;
            return StatusCode::SUCCESS;
        }
        // Call to GenericAddress routine
        pAdd->release();
    }
    return StatusCode::FAILURE;
}


