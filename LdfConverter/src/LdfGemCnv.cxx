// File and Version Information:
//      $Header$
//
// Description:
//      LdfGemCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfGemCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "ldfReader/data/LatData.h"
#include "LdfEvent/Gem.h"
#include "LdfBaseCnv.h"

class  LdfGemCnv : public LdfBaseCnv
//public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfGemCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfGemCnv(ISvcLocator* svc);

    virtual ~LdfGemCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_LdfGem;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_LdfGem;}

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
DECLARE_CONVERTER_FACTORY ( LdfGemCnv );

LdfGemCnv::LdfGemCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
///Converter(TEST_StorageType, CLID_LdfGem, svc)
{
    // Here we associate this converter with the /Event/Gem path on the TDS.
    m_path  = "/Event/Gem";
    declareObject("/Event/Gem", objType(), "PASS");
}

/*
StatusCode LdfGemCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfGemCnv::finalize() 
{
    return Converter::finalize();
}

*/

StatusCode LdfGemCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    LdfEvent::Gem *gem = new LdfEvent::Gem();
    refpObject = gem;
    // Retrieve the GEM data for this event 
    ldfReader::GemData ldfGem = ldfReader::LatData::instance()->getGem();
    ldfReader::GemDataTileList ldfTileList = ldfGem.tileList();
    LdfEvent::GemTileList tdsTileList(ldfTileList.XZM(), ldfTileList.XZP(),
            ldfTileList.YZM(), ldfTileList.YZP(), ldfTileList.XY(), 
            ldfTileList.RBN(), ldfTileList.NA());

    gem->initTrigger(ldfGem.tkrVector(), ldfGem.roiVector(), 
        ldfGem.calLEvector(),
        ldfGem.calHEvector(), ldfGem.cnoVector(), ldfGem.conditionSummary(), 
        ldfGem.missed(), tdsTileList);

    LdfEvent::GemOnePpsTime ppsTime(ldfGem.onePpsTime().timebase(),
        ldfGem.onePpsTime().seconds());

    LdfEvent::GemDataCondArrivalTime condArr;
    condArr.init(ldfGem.condArrTime().condArr());
    gem->initSummary(ldfGem.liveTime(), ldfGem.prescaled(),
            ldfGem.discarded(), condArr,
            ldfGem.triggerTime(), ppsTime, ldfGem.deltaEventTime(), ldfGem.deltaWindowOpenTime());
   

    return StatusCode::SUCCESS;
}

//StatusCode LdfGemCnv::updateObj(int* , LdfEvent::Gem*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.
//
//    MsgStream log(msgSvc(), "LdfGemCnv");
//    log << MSG::DEBUG << "LdfGemCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


