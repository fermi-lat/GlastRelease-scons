// File and Version Information:
//      $Header$
//
// Description:
//      LdfErrorCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfErrorCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "ldfReader/data/LatData.h"

#include "LdfEvent/ErrorData.h"

class  LdfErrorCnv : public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfErrorCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfErrorCnv(ISvcLocator* svc);

    virtual ~LdfErrorCnv()  { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_LdfErrorData;}
    static const unsigned char storageType() {return TEST_StorageType;}

    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_LdfErrorData;}

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
DECLARE_CONVERTER_FACTORY ( LdfErrorCnv );

LdfErrorCnv::LdfErrorCnv(ISvcLocator* svc) : Converter(TEST_StorageType, CLID_LdfErrorData, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    m_path = "/Event/Error";
}

StatusCode LdfErrorCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfErrorCnv::finalize() 
{
    return Converter::finalize();
}


StatusCode LdfErrorCnv::createObj(IOpaqueAddress* , 
                                    DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF diagnostic data on the TDS

    LdfEvent::ErrorData *errCol = new LdfEvent::ErrorData();
    refpObject = errCol; 
    
    // Retrieve the LAT data for this event and print out its contents
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    unsigned int i;
    for (i = 0; i < 16; i++) {
        if(myLatData->getTower(i)) {
            if (!myLatData->getTower(i)->getTem().exist()) continue;

            const ldfReader::ErrData* errLdf = myLatData->getTower(i)->getTem().getErr();

            if (!errLdf->exist()) continue;
            unsigned short towerId = myLatData->getTower(i)->getTowerId();
            LdfEvent::TowerErrorData err(towerId, errLdf->cal(), errLdf->tkr(), errLdf->phs(), errLdf->tmo());
            errCol->addTowerError(err);
        }
    }
    return StatusCode::SUCCESS;
}

//StatusCode LdfErrorCnv::updateObj(int* , LdfEvent::ErrorData* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "LdfErrorCnv");
//    log << MSG::DEBUG << "LdfErrorCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


