// File and Version Information:
//      $Header$
//
// Description:
//      LdfCnvSvc is the GLAST converter service.
//
// Author(s):

#define LdfCnvSvc_CXX

#include <iostream>
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/CnvFactory.h"
//#include "GaudiKernel/HashTable.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ICnvManager.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/RegistryEntry.h"
#include "GaudiKernel/GenericAddress.h"

#include "LdfBaseCnv.h"
//#include "Address.h"
#include "LdfCnvSvc.h"

//static const InterfaceID IID_IBaseEbfCnv(902, 1 , 0); 
static const InterfaceID IID_ILdfBaseCnv("ILdfBaseCnv", 1 , 0); 

// RCS Id for identification of object version
static const char* rcsid = "$Id$";


// Instantiation of a static factory class used by clients to create
// instances of this service
static const SvcFactory<LdfCnvSvc> s_LdfCnvSvcFactory;
const ISvcFactory& LdfCnvSvcFactory = s_LdfCnvSvcFactory;

LdfCnvSvc::LdfCnvSvc(const std::string& name, ISvcLocator* svc)
: ConversionSvc(name, svc, TEST_StorageType)              
{
}

StatusCode LdfCnvSvc::initialize()     {
    // Purpose and Method:  Setup GLAST's Event Converter Service.
    //   Associate LdfCnvSvc with the EventDataSvc
    //   Associate the list of known converters with this LdfCnvSvc.

    StatusCode status = ConversionSvc::initialize();
    if ( status.isSuccess() )   {
        //        ISvcLocator* svclocator = serviceLocator();  [unused for now]
        IDataProviderSvc *pIDP = 0;
        // Set event data service
        status = service("EventDataSvc", pIDP, true);
        if ( status.isSuccess() )   {
            status = setDataProvider ( pIDP );
        }
        else    {
            return status;
        }
        
        // Add known converters to the service: 
        ICnvManager::CnvIterator i, stop;
        for ( i = cnvManager()->cnvBegin(), stop = cnvManager()->cnvEnd(); i != stop; i++ )   {
            if ( repSvcType() == (*i)->repSvcType() )   {
                StatusCode iret = addConverter( (*i)->objType() );        
                if ( !iret.isSuccess() )   {
                    status = iret;
                }
            }
        }
        
        // Now we have to configure the map of leaves
        // Which should contain the association of converters with 
        // paths on the TDS
        for (LeafMap::iterator k = m_leaves.begin(); k != m_leaves.end(); k++ )   {
            std::string path = (*k).first;
            for (LeafMap::iterator j = m_leaves.begin(); j != m_leaves.end(); j++ )   {
                std::string pp = (*j).first.substr(0, (*j).first.rfind("/"));
                if ( path == pp && path != (*j).first )   {
                    (*k).second->push_back((*j).second);
                }
            }
        }
        
    }
    
    return status;
}

StatusCode LdfCnvSvc::finalize()     {
    StatusCode status = ConversionSvc::finalize();
    for (LeafMap::iterator k = m_leaves.begin(); k != m_leaves.end(); k++ )   {
        delete (*k).second;
    }
    m_leaves.erase(m_leaves.begin(), m_leaves.end());
    return status;
}

StatusCode LdfCnvSvc::declareObject(const ILdfCnvSvc::Leaf& leaf)   {
    // Purpose and Method:  Callback from each of the individual converters that allows
    //  association of TDS path and converter.
    Leaf* ll = new Leaf(leaf);
    std::pair<LeafMap::iterator, bool> p = m_leaves.insert(LeafMap::value_type( leaf.path, ll) );
    if( p.second )    {
        return StatusCode::SUCCESS;
    }
    delete ll;
    return StatusCode::FAILURE;
}

/// Update state of the service
StatusCode LdfCnvSvc::updateServiceState(IOpaqueAddress* pAddress)    {
    
    // not sure about the use of recid or bank...
    
    MsgStream log(msgSvc(), name());
    StatusCode status = INVALID_ADDRESS;
    IRegistry* ent = pAddress->registry();
    if ( 0 != ent )   {
        SmartIF<IDataManagerSvc> iaddrReg(IID_IDataManagerSvc, dataProvider());
 //       if ( 0 != iaddrReg )   {
            status = StatusCode::SUCCESS;
            std::string path = ent->identifier();
            LeafMap::iterator itm = m_leaves.find(path);
            //  int size = m_leaves.size();   [unused for now]
            if ( itm != m_leaves.end() )    {
                ILdfCnvSvc::Leaf* leaf = (*itm).second;
                if ( 0 != leaf )    {
                    //ipar[0] = leaf->userParameter;
                    for ( Leaf::iterator il = leaf->begin(); il != leaf->end(); il++ )   {
                        IOpaqueAddress* newAddr = 0;
                        unsigned long ipars[2] = {0, 0}; //{(*il)->userParameter, new_rid};
                        const std::string spars[2] = {"", ""}; //{par[0], (*il)->bank};
                        StatusCode ir =
                            addressCreator()->createAddress(TEST_StorageType, 
                            (*il)->clid, 
                            spars, 
                            ipars,
                            newAddr);
                        if ( ir.isSuccess() )   {
                            ir = iaddrReg->registerAddress((*il)->path, newAddr);
                            if ( !ir.isSuccess() )    {
                                newAddr->release();
                                status = ir;
                            }
                        }
                    }
                }
            }
//     }
//    else  {
//        status = NO_INTERFACE;
//    }
    }
    return status;
}

StatusCode LdfCnvSvc::queryInterface(const IID& riid, void** ppvInterface)  {
    if ( IID_ILdfBaseCnv == riid )  {
        *ppvInterface = (ILdfCnvSvc*)this;
    }
    else  {
        // Interface is not directly availible: try out a base class
        return ConversionSvc::queryInterface(riid, ppvInterface);
    }
    addRef();
    return StatusCode::SUCCESS;
}

StatusCode LdfCnvSvc::createAddress(unsigned char svc_type,
                                      const CLID& clid,
                                      const std::string* , 
                                      const unsigned long* ,
                                      IOpaqueAddress*& refpAddress) 
{
    MsgStream log( msgSvc(), name() );
    
    if (svc_type != repSvcType() ) {
        log << MSG::ERROR << "bad storage type" << svc_type << endreq;
        return StatusCode::FAILURE;
    }
/*
        Address* pAdd = new Address(svc_type, clid, par[0]);
        if ( pAdd )
        {
            refpAddress = pAdd;
            return StatusCode::SUCCESS;
        }
        // Call to GenericAddress routine
        pAdd->release();
        */

  refpAddress = new GenericAddress(TEST_StorageType,
                                   clid,
                                   "",  
                                   "",
                                   0);


    return StatusCode::SUCCESS;
}


