// File and Version Information:
// $Header$
//
// Description:
// Concrete converter for the McEvent header on the TDS /Event/MC

#define CNV_FILTERCNV_CPP 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "GaudiKernel/DataObject.h"

#include "LdfBaseCnv.h"
#include "FilterCnv.h"

//#include "GaudiKernel/ObjectVector.h"

// RCS Id for identification of object version
static const char* rcsid = "$Id$";

// Instantiation of a static factory class used by clients to create
// instances of this service
//static CnvFactory<FilterCnv> s_factory;
//const ICnvFactory& FilterCnvFactory = s_factory;
class  FilterCnv : public LdfBaseCnv
//public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<FilterCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    FilterCnv(ISvcLocator* svc);

    virtual ~FilterCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    //static const CLID&         classID()     {return CLID_DataObject;}
    static const CLID&         classID()     {return CLID_Filter;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_Filter;}
    //virtual const CLID& objType() const {return CLID_DataObject;}

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


DECLARE_CONVERTER_FACTORY ( FilterCnv );

FilterCnv::FilterCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
//Converter(TEST_StorageType, CLID_DataObject, svc)
{
  m_path = "/Event/Filter";
  declareObject("/Event/Filter", objType(), "PASS");
}

/*
StatusCode FilterCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode FilterCnv::finalize() 
{
    return Converter::finalize();
}
*/

StatusCode FilterCnv::createObj(IOpaqueAddress* pAddress, DataObject*& refpObject)
{
    // Purpose and Method:  Converter just creates an empty instance of MCEvent
    //   for the TDS.  The data members will be initialized by other components.

    refpObject = new DataObject();
    StatusCode sc=StatusCode::SUCCESS;
    return sc;
}





