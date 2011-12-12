// File and Version Information:
//      $Header$
//
// Description:
//      MetaEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define MetaEventCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "ldfReader/data/LatData.h"
#include "Event/TopLevel/Definitions.h"
#include "lsfData/LsfMetaEvent.h"
#include "LdfEvent/LsfMetaEvent.h"

class  MetaEventCnv : public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<MetaEventCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    MetaEventCnv(ISvcLocator* svc);

    virtual ~MetaEventCnv() {};

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_MetaEvent;}
    static const unsigned char storageType() {return TEST_StorageType;}

    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_MetaEvent;}

    /// Retrieve the class type of the data store the converter uses.
    // MSF: Masked to generate compiler error due to interface change
    virtual long repSvcType() const {return Converter::i_repSvcType();}

    /// Create the transient representation of an object.
    virtual StatusCode createObj(IOpaqueAddress* pAddress,DataObject*& refpObject);

    /// Methods to set and return the path in TDS for output of this converter
  //  virtual void setPath(const std::string& path) {m_path = path;}
    virtual const std::string& getPath() const    {return m_path;}

private:

    std::string m_path;

};

// Instantiation of a static factory class used by clients to create
// instances of this service
//static CnvFactory<MCEventCnv> s_factory;
//const ICnvFactory& MCEventCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY ( MetaEventCnv );


MetaEventCnv::MetaEventCnv(ISvcLocator* svc) : Converter(TEST_StorageType, CLID_MetaEvent, svc)
{
    // Here we associate this converter with the /Event/Gem path on the TDS.
    m_path = "/Event/MetaEvent";
}

StatusCode MetaEventCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode MetaEventCnv::finalize() 
{
    return Converter::finalize();
}


StatusCode MetaEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    const lsfData::MetaEvent& lsfMetaEvent = ldfReader::LatData::instance()->getMetaEvent();
    LsfEvent::MetaEvent *me = new LsfEvent::MetaEvent(lsfMetaEvent);
    refpObject = me;

    return StatusCode::SUCCESS;
}

//StatusCode MetaEventCnv::updateObj(int* , LsfEvent::MetaEvent*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "MetaEventCnv");
//    log << MSG::DEBUG << "MetaEventCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}//


