// File and Version Information:
//      $Header$
//
// Description:
//      AncillaryEventCnv is the concrete converter for the digi event header 
//      on the TDS /Event/AncillaryEvent
//
// Author(s):

#define AncillaryEventCnv_CXX 

//#include "GaudiKernel/RegistryEntry.h"
//#include "AncillaryEventCnv.h"
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "LdfBaseCnv.h"

#include "AdfEvent/AdfEvent.h"

class  AncillaryEventCnv : public LdfBaseCnv
//public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<AncillaryEventCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    AncillaryEventCnv(ISvcLocator* svc);

    virtual ~AncillaryEventCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_AncillaryEvent;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_AncillaryEvent;}

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
//static CnvFactory<AncillaryEventCnv> s_factory;
//const ICnvFactory& AncillaryEventCnvFactory = s_factory;

DECLARE_CONVERTER_FACTORY ( AncillaryEventCnv );

AncillaryEventCnv::AncillaryEventCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
//Converter(TEST_StorageType, CLID_AncillaryEvent, svc)
//: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/AncillaryEvent", objType(), "PASS");
    m_path="/Event/AncillaryEvent";
}
/*
StatusCode AncillaryEventCnv::initialize() { 
    return Converter::initialize();
}

StatusCode AncillaryEventCnv::finalize() { 
    return Converter::finalize();
}
*/

StatusCode AncillaryEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    DataObject *digi = new DataObject();
    refpObject = digi;
    return StatusCode::SUCCESS;
};

//StatusCode AncillaryEventCnv::updateObj(int* , DataObject* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "AncillaryEventCnv");
//    log << MSG::DEBUG << "AncillaryEventCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;//
//}


