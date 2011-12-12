// File and Version Information:
// $Header$
//
// Description:
// Concrete converter for the McEvent header on the TDS /Event/MC

#define CNV_MCEVENTCNV_CPP 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"

// Event for creating the McEvent stuff
#include "Event/TopLevel/MCEvent.h"

#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/ObjectVector.h"

// RCS Id for identification of object version
static const char* rcsid = "$Id$";

class  McEventCnv : public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<McEventCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    McEventCnv(ISvcLocator* svc);

    virtual ~McEventCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_McEvent;}
    static const unsigned char storageType() {return TEST_StorageType;}

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
  //  virtual void setPath(const std::string& path) {m_path = path;}
    virtual const std::string& getPath() const    {return m_path;}

private:

    std::string m_path;

};


McEventCnv::McEventCnv(ISvcLocator* svc) : Converter(TEST_StorageType, CLID_McEvent, svc)
{
  m_path = "/Event/MC";
}
StatusCode McEventCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode McEventCnv::finalize() 
{
    return Converter::finalize();
}


StatusCode McEventCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)
{
    // Purpose and Method:  Converter just creates an empty instance of MCEvent
    //   for the TDS.  The data members will be initialized by other components.

    refpObject = new Event::MCEvent();
    StatusCode sc=StatusCode::SUCCESS;
    return sc;
}


//const CLID& McEventCnv::classID(){ return Event::MCEvent::classID();}

