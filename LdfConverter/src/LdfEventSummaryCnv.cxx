// File and Version Information:
//      $Header$
//
// Description:
//      LdfEventSummaryCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfEventSummaryCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "Event/TopLevel/Event.h"
#include "ldfReader/data/LatData.h"
#include "LdfEvent/EventSummaryData.h"
#include "LdfBaseCnv.h"

class  LdfEventSummaryCnv : public LdfBaseCnv //public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfEventSummaryCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfEventSummaryCnv(ISvcLocator* svc);

    virtual ~LdfEventSummaryCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_LdfEventSummaryData;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_LdfEventSummaryData;}

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
DECLARE_CONVERTER_FACTORY ( LdfEventSummaryCnv );


LdfEventSummaryCnv::LdfEventSummaryCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
//Converter(TEST_StorageType, CLID_LdfEventSummaryData, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    m_path = "/Event/EventSummary";
    declareObject("/Event/EventSummary", objType(), "PASS");
}

/*
StatusCode LdfEventSummaryCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfEventSummaryCnv::finalize() 
{
    return Converter::finalize();
}
*/

StatusCode LdfEventSummaryCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    // Retrieve the LAT data for this event 
    MsgStream log(msgSvc(), "LdfEventSummaryCnv");
    unsigned int ebfSummary = ldfReader::LatData::instance()->summaryData().summary();
    LdfEvent::EventSummaryData *summary = new LdfEvent::EventSummaryData(ebfSummary);
    summary->initEventFlags(ldfReader::LatData::instance()->getEventFlags());

    //summary->initOswEvtSequence(ldfReader::LatData::instance()->getOsw().evtSequence());

    summary->initEventSizeInBytes(ldfReader::LatData::instance()->eventSizeInBytes());

    // Retrieve the contribution lengths and store them in the EventSummaryData
    unsigned int gemLen = ldfReader::LatData::instance()->getGem().lenInBytes();
    unsigned int oswLen = ldfReader::LatData::instance()->getOsw().lenInBytes();
    unsigned int aemLen = ldfReader::LatData::instance()->getAem().lenInBytes();
    unsigned int temLen[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned int diagLen[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned int errLen[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned int i;
    for (i = 0; i < 16; i++) {
        if(ldfReader::LatData::instance()->getTower(i)) {
            if (!ldfReader::LatData::instance()->getTower(i)->getTem().exist()) continue;
            temLen[i] = ldfReader::LatData::instance()->getTower(i)->getTem().lenInBytes();
            diagLen[i] = ldfReader::LatData::instance()->getTower(i)->getTem().getDiagnostic()->lenInBytes();
            errLen[i] = ldfReader::LatData::instance()->getTower(i)->getTem().getErr()->lenInBytes();
        }
    }

    summary->initContribLen(temLen, gemLen, oswLen, errLen, diagLen, aemLen);

    refpObject = summary;

    return StatusCode::SUCCESS;
}

//StatusCode LdfEventSummaryCnv::updateObj(int* , LdfEvent::EventSummaryData*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "LdfEventSummaryCnv");
//    log << MSG::DEBUG << "LdfEventSummaryCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


