// File and Version Information:
// $Header$
//
// Description:
// Concrete converter for the McEvent header on the TDS /Event/MC

#define CNV_MCEventCnv_CPP 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

#include "GlastSvc/EventSelector/IGlastCnv.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/MCEvent.h"

class  MCEventCnv : public Converter 
{

friend class CnvFactory<MCEventCnv>;


protected:
   /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    MCEventCnv(ISvcLocator* svc);

    virtual ~MCEventCnv();
public:

    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_McEvent;}
    static const unsigned char storageType() {return SICB_StorageType;}

    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_McEvent;}

    /// Retrieve the class type of the data store the converter uses.
    // MSF: Masked to generate compiler error due to interface change
    virtual long repSvcType() const {return Converter::i_repSvcType();}

    /// Create the transient representation of an object.
    virtual StatusCode createObj(IOpaqueAddress* pAddress,DataObject*& refpObject);

    /// Methods to set and return the path in TDS for output of this converter
    //virtual void setPath(const std::string& path) {m_path = path;}
    virtual const std::string& getPath() const    {return m_path;}

private:

    std::string m_path;

};

// Instantiation of a static factory class used by clients to create
// instances of this service
//static CnvFactory<MCEventCnv> s_factory;
DECLARE_CONVERTER_FACTORY ( MCEventCnv );


 MCEventCnv::MCEventCnv( ISvcLocator* svc) : Converter (SICB_StorageType, CLID_McEvent, svc) 
{
    m_path="/Event/MC";

    return;
}

MCEventCnv::~MCEventCnv() 
{
    return;
}

StatusCode MCEventCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode MCEventCnv::finalize() 
{
    return Converter::finalize();
}

// (To TDS) Conversion stuff
StatusCode MCEventCnv::createObj(IOpaqueAddress* /*addr*/, DataObject*& refpObject) 
{
    StatusCode ret = StatusCode::SUCCESS;

    refpObject = new Event::MCEvent();

    return ret;
}
