// File and Version Information:
//      $Header$
//
// Description:
//      LdfTimeCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfTimeCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "ldfReader/data/LatData.h"
#include "LdfEvent/LdfTime.h"

#include "LdfBaseCnv.h"

class  LdfTimeCnv : public LdfBaseCnv
// public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfTimeCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfTimeCnv(ISvcLocator* svc);

    virtual ~LdfTimeCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_LdfTime;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*

    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_LdfTime;}

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
DECLARE_CONVERTER_FACTORY ( LdfTimeCnv );

LdfTimeCnv::LdfTimeCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
//Converter(TEST_StorageType, CLID_LdfTime, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    m_path = "/Event/Time";
    declareObject("/Event/Time", objType(), "PASS");
}
StatusCode LdfTimeCnv::initialize() 
{
    StatusCode status = Converter::initialize();

/*
StatusCode LdfTimeCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfTimeCnv::finalize() 
{
    return Converter::finalize();
}
*/

StatusCode LdfTimeCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store Ebf time on the TDS.
    LdfEvent::LdfTime *time = new LdfEvent::LdfTime();
    refpObject = time;
    // Retrieve the LAT data for this event 
    unsigned int summary = ldfReader::LatData::instance()->summaryData().summary();
    time->initialize(ldfReader::LatData::instance()->summaryData().timeSec(), 
        ldfReader::LatData::instance()->summaryData().timeNanoSec(),
        ldfReader::LatData::instance()->summaryData().upperPpcTimeBaseWord(),
        ldfReader::LatData::instance()->summaryData().lowerPpcTimeBaseWord());
    return StatusCode::SUCCESS;
}

//StatusCode LdfTimeCnv::updateObj(int* , LdfEvent::LdfTime*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "LdfTimeCnv");
//    log << MSG::DEBUG << "LdfTimeCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;//
//}


