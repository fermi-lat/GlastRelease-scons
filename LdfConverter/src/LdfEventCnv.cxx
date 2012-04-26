// File and Version Information:
//      $Header$
//
// Description:
//      LdfEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfEventCnv_CXX 

//#include "GaudiKernel/RegistryEntry.h"
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/IOpaqueAddress.h"
//#include "LdfEventCnv.h"
#include "GaudiKernel/MsgStream.h"
#include "Event/TopLevel/Event.h"
#include "LdfBaseCnv.h"

#include "ldfReader/data/LatData.h"
#include "facilities/Timestamp.h"
#include "astro/JulianDate.h"
#include "Event/Utilities/TimeStamp.h"
#include <math.h>

class  LdfEventCnv : public LdfBaseCnv // public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfEventCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfEventCnv(ISvcLocator* svc);

    virtual ~LdfEventCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_Event;}
    static const unsigned char storageType() {return TEST_StorageType;}

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_Event;}

   /// Retrieve the class type of the data store the converter uses.
    // MSF: Masked to generate compiler error due to interface change
//    virtual long repSvcType() const {return Converter::i_repSvcType();}

    /// Initialize the converter
 //   virtual StatusCode initialize();

    /// Initialize the converter
 //   virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
   // virtual const CLID& objType() const {return CLID_Event;}

    /// Retrieve the class type of the data store the converter uses.
    // MSF: Masked to generate compiler error due to interface change
//    virtual long repSvcType() const {return Converter::i_repSvcType();}

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
//static CnvFactory<LdfEventCnv> s_factory;
//const ICnvFactory& LdfEventCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY ( LdfEventCnv );

LdfEventCnv::LdfEventCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
//: Converter(TEST_StorageType, CLID_Event, svc)
//: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event", objType(), "PASS");
    m_path = "/Event";
}

//StatusCode LdfEventCnv::initialize() {
    //StatusCode status = Converter::initialize();
//    StatusCode status = LdfBaseCnv::inialize();
//    return status;
//}

//StatusCode LdfEventCnv::finalize() {
//    StatusCode status = Converter::finalize();
//    return status;
//}

StatusCode LdfEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an empty EventHeader on
    //   the TDS.
    MsgStream log(msgSvc(), "LdfEventCnv");
    log << MSG::DEBUG << "LdfEventCnv::createObj called" << endreq;
    Event::EventHeader *header = new Event::EventHeader();
    refpObject = header;
    // Retrieve the LAT data for this event 
    unsigned int summary = ldfReader::LatData::instance()->summaryData().summary();
    // According to LAT Inter-module Communications Reference Manual
    // the event number is a 17 bit value, the 2 least significant are the
    // tag and the 15 most significant are in the eventNumber field of the summary
    header->setEvent(ldfReader::LatData::instance()->eventId());
    header->setRun(ldfReader::LatData::instance()->runId());
	
    // Also set the time in the Event::EventHeader
//    unsigned upperPpcTime = ldfReader::LatData::instance()->summaryData().upperPpcTimeBaseWord();
//    unsigned lowerPpcTime = ldfReader::LatData::instance()->summaryData().lowerPpcTimeBaseWord();
//    const double sixteenMHz = 16000000.;
    // To handle th 64-bit value - we separate the computation
    // The upper 32 bits would have to be shifted by 31 (or multiplied)
    // by 2^32 - we divide this by 16000000 to get the upperMultiplier
//    const double upperMultiplier = 268.435456; 
//    double lower = lowerPpcTime / sixteenMHz;
//    double upper = upperPpcTime * upperMultiplier;
//    double ppcSeconds =  upper + lower;
//    float ppcWholeSeconds = floor(ppcSeconds);
//    double frac = ppcSeconds - ppcWholeSeconds;
//    int ppcNanoSec = frac / 0.000000001;
//    log << MSG::DEBUG << "ppcSeconds: " << ppcSeconds << " frace = " << frac << endreq;
 //   log << MSG::DEBUG << "ppcNanoSec: " << ppcNanoSec << endreq;

    // To eliminate duplicate times, we make use of the PPC time as suggested 
    // by Jim Panetta.  Here we ignore the fractional seconds from the "real"
    // absolute time and add in the fractional seconds as computed by the PPC
    // word.  This should produce a unique time stamp for each event.
//    log << MSG::DEBUG << "seconds: " << ldfReader::LatData::instance()->summaryData().timeSec() << endreq;
//    log << MSG::DEBUG << "nanosec: " << ldfReader::LatData::instance()->summaryData().timeNanoSec() << endreq;
//    facilities::Timestamp facTimeStamp(
//             ldfReader::LatData::instance()->summaryData().timeSec(), 
//             ppcNanoSec);
//    double julTimeStamp = facTimeStamp.getJulian();
//    astro::JulianDate julDate(julTimeStamp);
    // Find number of seconds since missionStart
//    double tdsTimeInSeconds = julDate.seconds() - astro::JulianDate::missionStart().seconds();
//    TimeStamp tdsTimeStamp(tdsTimeInSeconds);
    TimeStamp tdsTimeStamp(ldfReader::LatData::instance()->getTimeInSecTds());
    header->setTime(tdsTimeStamp);

    // Set the LiveTime on the TDS to -1 
    // Need the livetime DB to set livetime since mission start, so we'll just 
    // set it to some dummy value for now
    header->setLivetime(-1.);


    // initialize the trigger to zero to promote processing if TriggerAlg is
    // run
    header->setTrigger(-1);

    return StatusCode::SUCCESS;
}

//StatusCode LdfEventCnv::updateObj(int* , Event::EventHeader*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "LdfEventCnv");
//    log << MSG::DEBUG << "LdfEventCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


