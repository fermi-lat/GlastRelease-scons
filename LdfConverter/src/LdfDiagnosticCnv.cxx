// File and Version Information:
//      $Header$
//
// Description:
//      LdfDiagnosticCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define LdfDiagnosticCnv_CXX 

#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "ldfReader/data/LatData.h"

#include "LdfEvent/DiagnosticData.h"
#include "LdfBaseCnv.h"

class LdfDiagnosticCnv : public LdfBaseCnv
//public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<LdfDiagnosticCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    LdfDiagnosticCnv(ISvcLocator* svc);

    virtual ~LdfDiagnosticCnv() { };

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_LdfDiagnosticData;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_LdfDiagnosticData;}

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
DECLARE_CONVERTER_FACTORY ( LdfDiagnosticCnv );

LdfDiagnosticCnv::LdfDiagnosticCnv(ISvcLocator* svc) 
: LdfBaseCnv(classID(), svc)
//: Converter(TEST_StorageType, CLID_LdfDiagnosticData, svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    m_path = "/Event/Diagnostic";
    declareObject("/Event/Diagnostic", objType(), "PASS");
}

/*
StatusCode LdfDiagnosticCnv::initialize() 
{
    StatusCode status = Converter::initialize();

    return status;
}

StatusCode LdfDiagnosticCnv::finalize() 
{
    return Converter::finalize();
}
*/

StatusCode LdfDiagnosticCnv::createObj(IOpaqueAddress* , 
                                    DataObject*& refpObject) {
    // Purpose and Method:  This converter will store LDF diagnostic data on the TDS
    LdfEvent::DiagnosticData *diag = new LdfEvent::DiagnosticData();
    refpObject = diag;
    
    // Retrieve the LAT data for this event and print out its contents
    ldfReader::LatData* myLatData = ldfReader::LatData::instance();
    unsigned int i;
    for (i = 0; i < 16; i++) {
        if(myLatData->getTower(i)) {
            if (!myLatData->getTower(i)->getTem().exist()) continue;

            const ldfReader::DiagnosticData* diagnostic = myLatData->getTower(i)->getTem().getDiagnostic();

            if (!diagnostic->exist()) continue;

	    // Added 03.18.2005 by awb to get all Tkr and Cal diagnostics contributions. 
            // Had to hardcode '8' - this should be fixed later.  
            unsigned int numCal = diagnostic->getNumCalDiagnostic();
	    unsigned int ind;
	    for (ind = 0; ind < numCal; ind++) {
	      ldfReader::CalDiagnosticData ldfCalDiag = diagnostic->getCalDiagnosticByIndex(ind);
	      //LdfEvent::CalDiagnosticData cal(diagnostic->getCalDiagnosticByIndex(ind).dataWord());
	      LdfEvent::CalDiagnosticData cal(ldfCalDiag.dataWord(), ldfCalDiag.tower(), ldfCalDiag.layer());
	      diag->addCalDiagnostic(cal);
	    }
            unsigned int numTkr = diagnostic->getNumTkrDiagnostic();
	    for (ind = 0; ind < numTkr; ind++) {
	      ldfReader::TkrDiagnosticData ldfTkrDiag = diagnostic->getTkrDiagnosticByIndex(ind);
	      LdfEvent::TkrDiagnosticData tkr(ldfTkrDiag.dataWord(), ldfTkrDiag.tower(), ldfTkrDiag.gtcc());
	      diag->addTkrDiagnostic(tkr);
	    }

	    // Taken out 03.18.2005 by awb:
            //ldfReader::CalDiagnosticData ldfCalDiag = diagnostic->getCalDiagnostic();
            //LdfEvent::CalDiagnosticData cal(ldfCalDiag.dataWord(), ldfCalDiag.tower(), ldfCalDiag.layer());
            //diag->addCalDiagnostic(cal);

            //ldfReader::TkrDiagnosticData ldfTkrDiag = diagnostic->getTkrDiagnostic();
            // LdfEvent::TkrDiagnosticData tkr(ldfTkrDiag.dataWord(), ldfTkrDiag.tower(), ldfTkrDiag.gtcc());
            // diag->addTkrDiagnostic(tkr);
        }
    }
    return StatusCode::SUCCESS;
}

//StatusCode LdfDiagnosticCnv::updateObj(int* , LdfEvent::DiagnosticData* ) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "LdfDiagnosticCnv");
//    log << MSG::DEBUG << "LdfDiagnosticCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;//
//}


