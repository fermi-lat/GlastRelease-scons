// File and Version Information:
//      $Header$
//
// Description:
//      CcsdsCnv is the concrete converter for the CCSDS items on TDS 
//
// Author(s):

#define EbfCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "ldfReader/data/LatData.h"
#include "lsfData/Ebf.h"
#include "EbfWriter/Ebf.h"

#include "LdfBaseCnv.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
//static CnvFactory<EbfCnv> s_factory;
//const ICnvFactory& EbfCnvFactory = s_factory;
class  EbfCnv : public LdfBaseCnv
//public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<EbfCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    EbfCnv(ISvcLocator* svc);

    virtual ~EbfCnv() {};

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_Ebf;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_Ebf;}

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
DECLARE_CONVERTER_FACTORY ( EbfCnv );

EbfCnv::EbfCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)//Converter (TEST_StorageType, CLID_Ebf, svc)
{
    // Here we associate this converter with the /Event/Filter/Ebf path on the TDS.
    m_path = "/Event/Filter/Ebf";
    declareObject("/Event/Filter/Ebf", objType(), "PASS");
    return;
}

/*
StatusCode EbfCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode EbfCnv::finalize() 
{
    return Converter::finalize();
}
*/

StatusCode EbfCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF Gem on the TDS.
    MsgStream log(msgSvc(), "EbfCnv");
    log << MSG::DEBUG << "EbfCnv::createObj" << endreq;
    const lsfData::Ebf& lsfEbf = ldfReader::LatData::instance()->getEbf();
    EbfWriterTds::Ebf *me = new EbfWriterTds::Ebf();
    unsigned int len;
    char *data = lsfEbf.get(len);
    log << MSG::DEBUG << "EBF len " << len  << " first byte "
        << std::hex << (unsigned int)(data[0]) << std::dec << endreq;
    me->set(data,len);
    refpObject = me;

    return StatusCode::SUCCESS;
}

//StatusCode EbfCnv::updateObj(int* , EbfWriterTds::Ebf*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "EbfCnv");
//    log << MSG::DEBUG << "EbfCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}//

