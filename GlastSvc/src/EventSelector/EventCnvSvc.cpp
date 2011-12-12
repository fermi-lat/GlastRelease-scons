// File and Version Information:
//      $Header$
//
// Description:
//      EventCnvSvc is the GLAST converter service.
//
// Author(s):

#define EVENTCNVSVC_CPP

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/MCEvent.h"
#include "GaudiUtils/IIODataManager.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/IPersistencySvc.h"
#include "GaudiKernel/GenericAddress.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/strcasecmp.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IRegistry.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IAddressCreator.h"
//#include "GaudiKernel/IConverter.h"

//#include "EventCnvSvc.h"

#include <iostream>
#include "GaudiKernel/ConversionSvc.h"
//#include "GaudiKernel/SvcFactory.h"
//#include "GaudiKernel/CnvFactory.h"
//#include "GaudiKernel/MsgStream.h"
//#include "GaudiKernel/SmartIF.h"
//#include "GaudiKernel/ICnvManager.h" HMK No longer availablein Gaudi v21r7
//#include "GaudiKernel/ISvcLocator.h"
//#include "GaudiKernel/IDataProviderSvc.h"
//#include "GaudiKernel/IDataManagerSvc.h"
//#include "GaudiKernel/RegistryEntry.h"

//#include "Address.h"

//#include "GlastSvc/EventSelector/IGlastCnv.h"

//#include <map>


template <class TYPE> class SvcFactory;


//static const InterfaceID IID_EventCnvSvc("EventCnvSvc", 1, 0);

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

class EventCnvSvc  : public ConversionSvc	
{

  friend class SvcFactory<EventCnvSvc>;

public:

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  virtual StatusCode updateServiceState(IOpaqueAddress* pAddress);

  /// Override inherited queryInterface due to enhanced interface
  virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);

    /// Concrete class type
   virtual const CLID& objType() const;


   /** Object creation callback
      *  @param[in]   pAddr     Reference to opaque object address
      *  @param[out]  refpObj   Reference to store pointer to created object
      *
      *  @return Status code indicating success or failure
      */
  //virtual StatusCode createObj(IOpaqueAddress* pAddr, DataObject*& refpObj);


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

    /// Retrieve converter from list
// virtual IConverter* converter(const CLID& clid);

protected:

  EventCnvSvc(const std::string& name, ISvcLocator* svc);

  virtual ~EventCnvSvc() { };

  /// Reference to data manager interface
 IDataManagerSvc*       m_dataMgr;
  /// Reference to file manager service
 // Gaudi::IIODataManager* m_ioMgr;


private:
    //typedef std::map<std::string, IGlastCnv*>                PathToCnvMap;
    //typedef std::map<std::string, std::vector<std::string> > SubPathMap;

    //PathToCnvMap m_pathToCnvMap;
    //SubPathMap   m_subPathMap;
    IAddressCreator *m_addressCreator;
};

//#include "GlastSvc/EventSelector/IGlastCnv.h"

//static const InterfaceID IID_IBaseCnv(902, 1 , 0); 
// RCS Id for identification of object version
static const char* rcsid = "$Id$";


// Instantiation of a static factory class used by clients to create
// instances of this service
//static const SvcFactory<EventCnvSvc> s_EventCnvSvcFactory;
//const ISvcFactory& EventCnvSvcFactory = s_EventCnvSvcFactory;

DECLARE_SERVICE_FACTORY( EventCnvSvc );

EventCnvSvc::EventCnvSvc(const std::string& name, ISvcLocator* svc)
: ConversionSvc(name, svc, SICB_StorageType),m_dataMgr(0)
{
  //  m_pathToCnvMap.clear();
  //  m_subPathMap.clear();
}

StatusCode EventCnvSvc::initialize()     {
    // Purpose and Method:  Setup GLAST's Event Converter Service.
    //   Associate EventCnvSvc with the EventDataSvc
    //   Associate the list of known converters with this EventCnvSvc.
    MsgStream log(msgSvc(), name());

    StatusCode status = ConversionSvc::initialize();
    if ( !status.isSuccess() ) return status;


  IPersistencySvc *pSvc = 0;
  status = service("EventPersistencySvc",pSvc,true);
  if ( !status.isSuccess() )  {
    log << MSG::ERROR << "Unable to localize EventPersistencySvc." << endmsg;
    return status;
  }

  status = pSvc->addCnvService( this );
  if ( !status.isSuccess() )  {
    log << MSG::ERROR << "Unable to add conversion service" << endmsg;
    return status;
  }

//  m_rootCLID = CLID_Event;

  m_addressCreator = 0;
  //pSvc->addrCreator();
   status = service("EventPersistencySvc",m_addressCreator);
  if ( !status.isSuccess() )  {
    log << MSG::ERROR << "Unable to localize IAddressCreator." << endmsg;
    return status;
  }


  status = service ( "EventDataSvc" , m_dataMgr , true ) ;

  // get the IDataManagerSvc interface from the EventPersistencySvc
  //status = dataProvider()->queryInterface(IDataManagerSvc::interfaceID(),
  //                                    (void**)&m_dataMgr);
  if ( !status.isSuccess() )  {
    log << MSG::ERROR << "Conversion service " << name() 
        << "not registered to EventPersistencySvc." << endmsg;
    return status;
  }

  //m_rootCLID = m_dataMgr->rootCLID();

  //IOpaqueAddress* newAddr = 0;
  //unsigned long ipars[2] = {0, 0}; //{(*il)->userParameter, new_rid};
  //const std::string spars[2] = {"", ""}; //{par[0], (*il)->bank};
  //StatusCode ir = m_addressCreator->createAddress(SICB_StorageType,
  //    CLID_Event, //glastCnv->objType(), 
  //    spars, 
  //    ipars,
  //    newAddr);
  //if ( ir.isSuccess() )   
  //{
  //    ir = m_dataMgr->registerAddress("/Event", newAddr);
  //    if ( !ir.isSuccess() )    
  //    {
  //        newAddr->release();
  //        status = ir;
   //   }
  //}
  //newAddr = 0;
  //ir = m_addressCreator->createAddress(SICB_StorageType,
  //    CLID_McEvent, //glastCnv->objType(), 
  //    spars, 
  //    ipars,
  //    newAddr);
  //if ( ir.isSuccess() )   
  //{
  //    ir = m_dataMgr->registerAddress("/Event/MC", newAddr);
  //    if ( !ir.isSuccess() )    
  //    {
  //        newAddr->release();
  //        status = ir;
   //   }
  //}



  // Retrieve conversion service handling event iteration
//  status = Service::service(m_ioMgrName, m_ioMgr);
//  if( !status.isSuc/cess() ) {
//    log << MSG::ERROR 
//        << "Unable to localize interface IID_IIODataManager from service:" 
//        << m_ioMgrName << endmsg;
//    return status;
 // }

    /*
    if ( status.isSuccess() )   
    {
        // HMK Unused ISvcLocator* svclocator = serviceLocator();
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
                // HMK Unused int j = 0;
            }
        }
    }
    */
    
    return status;
}

StatusCode EventCnvSvc::finalize()     {
  if ( m_dataMgr ) m_dataMgr->release();
  m_dataMgr = 0;
//  if ( m_ioMgr ) m_ioMgr->release();
//  m_ioMgr = 0;
      StatusCode status = ConversionSvc::finalize();

    return status;
}

/// Update state of the service

StatusCode EventCnvSvc::updateServiceState(IOpaqueAddress* pAddress)    
{    
    // not sure about the use of recid or bank...
    
    MsgStream log(msgSvc(), name());
    static bool flag = false;
    StatusCode status = INVALID_ADDRESS;
    IRegistry* ent    = pAddress->registry();

//    if (flag) return StatusCode::SUCCESS;
//    flag = true;

    // Need a registry entry to proceed
    if ( 0 != ent )   
    {
//        SmartIF<IDataManagerSvc> iaddrReg(IID_IDataManagerSvc, dataProvider());
//        SmartIF<IDataManagerSvc> iaddrReg(IDataManagerSvc::interfaceID(), dataProvider());

        status = StatusCode::SUCCESS;

        // Path we are updating from
        std::string path = ent->identifier();

        if (path != "/Event") return StatusCode::SUCCESS;

        // Find the list of possible daughter paths
        //SubPathMap::iterator pathIter = m_subPathMap.find(path);

        // Check to be sure something was found
        //if (pathIter != m_subPathMap.end())
       // {
            // Loop over the list of daughter paths
         //   for(std::vector<std::string>::iterator subIter = pathIter->second.begin(); subIter != pathIter->second.end(); subIter++)
         //   {
                // Now look up the pointer to the converter
           //     PathToCnvMap::iterator cnvIter  = m_pathToCnvMap.find(*subIter);
            //    IGlastCnv*             glastCnv = cnvIter->second;

            //    if (glastCnv->getPath() == path) continue;

              //  if (glastCnv)
               // {
                    IOpaqueAddress* newAddr = 0;
                    unsigned long ipars[2] = {0, 0}; //{(*il)->userParameter, new_rid};
                    const std::string spars[2] = {"", ""}; //{par[0], (*il)->bank};
                    StatusCode ir =
            
                    //addressCreator()->createAddress(SICB_StorageType, 
                    m_addressCreator->createAddress(SICB_StorageType,
                                                    CLID_McEvent, //glastCnv->objType(), 
                                                    spars, 
                                                    ipars,
                                                    newAddr);
                    if ( ir.isSuccess() )   
                    {
                        //ir = iaddrReg->registerAddress(glastCnv->getPath(), newAddr);
                        ir = m_dataMgr->registerAddress("/Event/MC", newAddr);
                        if ( !ir.isSuccess() )    
                        {
                            newAddr->release();
                            status = ir;
                        }
                    }
                //}
            //}
       // }
    }
    return status;
}


/// Concrete class type  Do we need this?  is that the right type?
const CLID& EventCnvSvc::objType() const  {
  return DataObject::classID();     
}


StatusCode EventCnvSvc::queryInterface(const InterfaceID& riid, void** ppvInterface)  {
//    if ( IID_EventCnvSvc == riid )  {
//    if ( EventCnvSvc::interfaceID() == riid )  {
//    if (IID_EventCnvSvc.versionMatch(riid) ) {
 //       *ppvInterface = (EventCnvSvc*)this;
  //      addRef();
  //      return StatusCode::SUCCESS;
  //  }
    // Interface is not directly availible: try out a base class
    return ConversionSvc::queryInterface(riid,ppvInterface);
}




/// Object creation callback
/*StatusCode EventCnvSvc::createObj(IOpaqueAddress* pA, DataObject*& refpObj) {
  if ( pA )  {
    //bool default_reader = m_bankLocation == RawEventLocation::Default;
    // default_reader = 
    //if ( default_reader )   {
      if ( pA->clID() == CLID_DataObject )  {
        refpObj = new DataObject();
        return StatusCode::SUCCESS;
      }  else if (pA->clID() == CLID_Event )  {
         refpObj = new Event::EventHeader();
         return StatusCode::SUCCESS;
      } else if ( pA->clID() == CLID_McEvent) {
          refpObj = new Event::MCEvent();
          return StatusCode::SUCCESS;
      }

     // RawDataAddress* rawAdd = dynamic_cast<RawDataAddress*>(pA);
    //  rawAdd->setSvcType(ROOT_StorageType);
   // return dataProvider()->retrieveObject("/Event",refpObj);
   // }
  }
  return StatusCode::FAILURE;
}
*/

StatusCode EventCnvSvc::createAddress( long svc_type,
                                       const CLID& clid,
                                       const std::string* par, 
                                       const unsigned long* /*ip*/,
                                       IOpaqueAddress*& refpAddress)
{
    MsgStream log(msgSvc(), name() );

    if ( svc_type == repSvcType() )   
    {
        // Now create the address
        refpAddress = new GenericAddress( svc_type, clid);
        return StatusCode::SUCCESS;

    }
    log << MSG::ERROR << "Unable to createAddress" << endmsg;
    return StatusCode::FAILURE;
}


//IConverter* EventCnvSvc::converter(const CLID& clid) {
//  IConverter* cnv = 0;
//  cnv = ConversionSvc::converter(clid);
//  if ( cnv ) {
//    return cnv;
//  }
//  else {
//    return ConversionSvc::converter(CLID_Any);
 // }  
//}

