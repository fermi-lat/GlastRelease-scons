// File and Version Information:
//      $Header$
//
// Description:
//      LdfDigiCnv is the concrete converter for the digi event header on the TDS /Event/Digi
//
// Author(s):

#define LdfDigiCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "Event/TopLevel/DigiEvent.h"
//#include "Event/TopLevel/DigiEvent.h"

class  LdfDigiCnv : public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfDigiCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfDigiCnv(ISvcLocator* svc);

    virtual ~LdfDigiCnv() {  };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_DigiEvent;}
    static const unsigned char storageType() {return TEST_StorageType;}

    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_DigiEvent;}

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
DECLARE_CONVERTER_FACTORY ( LdfDigiCnv );

LdfDigiCnv::LdfDigiCnv(ISvcLocator* svc) : Converter(TEST_StorageType, CLID_DigiEvent, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    m_path = "/Event/Digi";
}

StatusCode LdfDigiCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfDigiCnv::finalize() 
{
    return Converter::finalize();
}

StatusCode LdfDigiCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    Event::DigiEvent *digi = new Event::DigiEvent();
    // Set fromMc to false
    digi->initialize(false);
    refpObject = digi;
    return StatusCode::SUCCESS;
};

//StatusCode LdfDigiCnv::updateObj(int* , Event::DigiEvent* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.
//
//    MsgStream log(msgSvc(), "LdfDigiCnv");
//    log << MSG::DEBUG << "LdfDigiCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;//
//}

