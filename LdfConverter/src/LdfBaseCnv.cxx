// File and Version Information:
//      $Header$
//
// Description:
//      LdfBaseCnv is the base class defining all GLAST converters.
//
// Author(s):

#define _LdfBaseCnv_CXX

#include "GaudiKernel/IService.h"
#include "GaudiKernel/ISvcLocator.h"

//#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/MsgStream.h"

#include "LdfBaseCnv.h"

static const InterfaceID IID_ILdfBaseCnv("ILdfBaseCnv", 1 , 0); 

LdfBaseCnv::LdfBaseCnv(const CLID& clid, ISvcLocator* svc)
: Converter(TEST_StorageType, clid, svc), m_CnvSvc(0)   {
    
    // Declarations are commented out since variables are currently unused
    //  StatusCode sc;
    //  MsgStream log(msgSvc(), "LdfBaseCnv");
    
    // provide access to the Glast Detector Service, so that we may call the 
    //  GlastDetSvc::accept method from within our converters
    //IService *isvc=0;
    //sc = serviceLocator()->getService ("GlastDetSvc", isvc, true);
    //if (sc.isSuccess()) {
    //    sc = isvc->queryInterface(IID_IGlastDetSvc, (void**)&m_detSvc);
    //}
    //if(sc.isFailure()){
    //    log << MSG::ERROR << "Unable start Glast detector service within LdfBaseCnv" << endreq;
    //} 
}


StatusCode LdfBaseCnv::createRep(DataObject* , 
                              IOpaqueAddress*& )   {
    // Purpose and Method: Convert the transient object to the requested 
    //     representation.  It is expected that derived classes will override
    //     this method.
    return StatusCode::FAILURE;
}

StatusCode LdfBaseCnv::fillRepRefs(IOpaqueAddress* ,
                                DataObject* )    {
    // Purpose and Method:  Resolve the references of the converted object.
    //     It is expected that derived classes will override this method.
    return StatusCode::FAILURE;
}

StatusCode LdfBaseCnv::updateRep(IOpaqueAddress* , 
                              DataObject* )   {
    // Purpose and Method:  Update the converted representation of a transient 
    //     object.  It is expected that derived classes will override this.
    return StatusCode::FAILURE;
}

StatusCode LdfBaseCnv::updateRepRefs(IOpaqueAddress* , 
                                  DataObject* ) {
    // Purpose and Method:  Update the references of an already converted object.
    //   It is expected that derived classes will override this method.
    return StatusCode::FAILURE;
}


StatusCode LdfBaseCnv::initialize()   {
    // Purpose and Method:  Perform standard converter initialization.
    //   Access the EventCnvSvc to create an association between converters 
    //   and paths within the TDS, using the vector of leaves and the
    //   declareObject methods available in each specific converter.
    StatusCode status = Converter::initialize();
    if ( status.isSuccess() )   {
        IService* isvc = 0;
        status = serviceLocator()->service("EventCnvSvc", isvc, true);
        if ( status.isSuccess() )   {
            status = isvc->queryInterface(IID_ILdfBaseCnv, (void**)&m_CnvSvc);
            if ( status.isSuccess() )   {
                for ( std::vector<ILdfCnvSvc::Leaf>::iterator i = m_leaves.begin(); i != m_leaves.end(); i++ )    {
                    m_CnvSvc->declareObject(*i);
                }
            }
        }
    }
    return status;
}

StatusCode LdfBaseCnv::finalize()   {
    if ( m_CnvSvc )     {
        m_CnvSvc->release();
        m_CnvSvc = 0;
    }
    return Converter::finalize();
}

void LdfBaseCnv::declareObject(const std::string& path, const CLID& cl, 
                            const std::string& bank, long par)  {
    // Purpose and Method:  Save the path on the TDS, in the m_leaves vector, 
    //   corresponding to the DataObject that the converter handles.
    m_leaves.push_back(ILdfCnvSvc::Leaf(path, cl, bank, par));
}


