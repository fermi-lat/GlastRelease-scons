// File and Version Information:
//      $Header$
//
// Description:
//      AdfEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define AdfEventCnv_CXX 

//#include "GaudiKernel/RegistryEntry.h"
//#include "AdfEventCnv.h"
#include "GaudiKernel/Converter.h"
#include "GaudiKernel/CnvFactory.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IOpaqueAddress.h"
#include "AdfEvent/AdfEvent.h"

//#include "GaudiKernel/MsgStream.h"

#include "AdfEvent/AncillaryWord.h"
#include "ldfReader/data/LatData.h"

#include "LdfBaseCnv.h"

class  AdfEventCnv :public LdfBaseCnv
// public Converter //virtual public IGlastCnv, public Converter 
{

  friend class CnvFactory<AdfEventCnv>;


protected:

    /**
        Constructor for this converter
        @param svc a ISvcLocator interface to find services
        @param clid the type of object the converter is able to convert
    */
    AdfEventCnv(ISvcLocator* svc);

    virtual ~AdfEventCnv();

public:
    /// Query interfaces of Interface
    //virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    static const CLID&         classID()     {return CLID_AncillaryEvent;}
    static const unsigned char storageType() {return TEST_StorageType;}

/*
    /// Initialize the converter
    virtual StatusCode initialize();

    /// Initialize the converter
    virtual StatusCode finalize();
*/

    /// Retrieve the class type of objects the converter produces. 
    virtual const CLID& objType() const {return CLID_AncillaryEvent;}

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
DECLARE_CONVERTER_FACTORY ( AdfEventCnv );


AdfEventCnv::AdfEventCnv(ISvcLocator* svc) : LdfBaseCnv(classID(), svc)
//Converter(TEST_StorageType, CLID_AncillaryEvent, svc)
//: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/AncillaryEvent/AdfEvent", objType(), "PASS");
    m_path = "/Event/AncillaryEvent/AdfEvent";
    return;
}

AdfEventCnv::~AdfEventCnv() { }

/*
StatusCode AdfEventCnv::initialize() { return Converter::initialize(); }
StatusCode AdfEventCnv::finalize() { return Converter::finalize(); }
*/

StatusCode AdfEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an AdfEvent on
    //   the TDS.
    MsgStream log(msgSvc(), "AdfEventCnv");
    // Retrieve the LAT data for this event 
    const ldfReader::AdfData& adfLdf = ldfReader::LatData::instance()->getAdf();
    if (!adfLdf.exist()){
        refpObject = 0;
        return StatusCode::FAILURE;
    }

    AncillaryData::AdfEvent *adf = new AncillaryData::AdfEvent();
    refpObject  = adf;
    const unsigned char* buffer = adfLdf.buffer();
    unsigned int *header = &((unsigned int*)buffer)[0];
    AncillaryData::EventSummaryData summaryData;
    summaryData.setData(header);
    adf->setEventSummaryData(summaryData);

  unsigned remainingBytes = adf->getEventSummaryData().getEventLength(); 
  remainingBytes-=4*AncillaryData::ANCILLARY_WORD_LENGTH;
  if(remainingBytes<=0) return StatusCode::SUCCESS;
  unsigned int index = 4;
  while(remainingBytes>0)
    {
      unsigned int rawword = header[index++]; //readRawWord();
      remainingBytes-=AncillaryData::ANCILLARY_WORD_LENGTH;
      AncillaryData::AncillaryWord OneWord;
      OneWord.setData(rawword);
      unsigned int Header = OneWord.getHeader();
      AncillaryData::QdcHeaderWord qdcHW;
      AncillaryData::FadcHeaderWord fadcHW;
      AncillaryData::ScalerHeaderWord scalerHW;
      
      switch (Header)
	{
	case AncillaryData::ANCILLARY_QDC_HID:
	  qdcHW.setData(rawword);

	  for(int i=0; i < qdcHW.getQdcFifo(); i++)
	    {
	      AncillaryData::QdcDataWord qdcDW;
	      qdcDW.setData(header[index++]);
	      qdcDW.checkHeader();
	      adf->appendQdcDataWord(qdcDW);
	    }
	  remainingBytes-=qdcHW.getQdcFifo()*AncillaryData::ANCILLARY_WORD_LENGTH;
	  break;
	case AncillaryData::ANCILLARY_FADC_HID:
	  fadcHW.setData(rawword);
	  for(int i=0; i < fadcHW.getFadcFifo(); i++)
	    {
	      AncillaryData::FadcDataWord fadcDW;
	      fadcDW.setData(header[index++]);
	      fadcDW.checkHeader();
	      adf->appendFadcDataWord(fadcDW);
	    }
	  remainingBytes-=fadcHW.getFadcFifo()*AncillaryData::ANCILLARY_WORD_LENGTH;
	  break;
	case AncillaryData::ANCILLARY_SCALER_ID:
	  scalerHW.setData(rawword);
	  for(int i=0; i < scalerHW.getScalerCounters(); i++)
	    {
	      AncillaryData::ScalerDataWord scalerDW;
	      scalerDW.setData(header[index++]);
	      adf->appendScalerDataWord(scalerDW);
	    }
	  remainingBytes-= scalerHW.getScalerCounters()*AncillaryData::ANCILLARY_WORD_LENGTH;
	  break;
	default:
	  break;
	}
    }


    return StatusCode::SUCCESS;
}

//StatusCode AdfEventCnv::updateObj(int* , AncillaryData::AdfEvent*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

//    MsgStream log(msgSvc(), "AdfEventCnv");
//    log << MSG::DEBUG << "AdfEventCnv::updateObj" << endreq;
//    return StatusCode::SUCCESS;
//}


