// File and Version Information:
//      $Header$
//
// Description:
//      CcsdsCnv is the concrete converter for the CCSDS items on TDS 
//
// Author(s):

#define CcsdsCnv_CXX 

//#include "GaudiKernel/RegistryEntry.h"
//#include "CcsdsCnv.h"
//#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "LdfEvent/LsfCcsds.h"
#include "LdfBaseCnv.h"


#include "ldfReader/data/LatData.h"
#include "lsfData/LsfCcsds.h"

class  CcsdsCnv : public LdfBaseCnv //public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<CcsdsCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    CcsdsCnv(ISvcLocator* svc);

    virtual ~CcsdsCnv(){};

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_LsfCcsds;}
    static const unsigned char storageType() {return TEST_StorageType;}

    /// Initialize the converter
    //virtual StatusCode initialize();

    /// Initialize the converter
    //virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_LsfCcsds;}

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
//static CnvFactory<CcsdsCnv> s_factory;
//const ICnvFactory& CcsdsCnvFactory = s_factory;
DECLARE_CONVERTER_FACTORY ( CcsdsCnv );

CcsdsCnv::CcsdsCnv(ISvcLocator* svc) 
: LdfBaseCnv(classID(), svc)
//: Converter(TEST_StorageType, CLID_LsfCcsds, svc)
{
    // Here we associate this converter with the /Event/Ccsds path on the TDS.
    declareObject("/Event/Ccsds", objType(), "PASS");
    m_path = "/Event/Ccsds";
}

/*
StatusCode CcsdsCnv::initialize() {
    return Converter::initialize();
}

StatusCode CcsdsCnv::finalize() {
    return Converter::finalize();
}
*/

StatusCode CcsdsCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    const lsfData::LsfCcsds& lsfCcsds = ldfReader::LatData::instance()->getCcsds();
    LsfEvent::LsfCcsds *me = new LsfEvent::LsfCcsds(lsfCcsds);
    refpObject = me;

    return StatusCode::SUCCESS;
}

//StatusCode CcsdsCnv::updateObj(int* , LsfEvent::LsfCcsds*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "CcsdsCnv");
//    log << MSG::DEBUG << "CcsdsCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


