// File and Version Information:
//      $Header$
//
// Description:
//      AdfEventCnv is the concrete converter for the event header on the TDS /Event
//
// Author(s):

#define AdfEventCnv_CXX 

#include "GaudiKernel/RegistryEntry.h"
#include "AdfEventCnv.h"
#include "GaudiKernel/MsgStream.h"

#include "AdfEvent/AncillaryWord.h"
#include "ldfReader/data/LatData.h"

// Instantiation of a static factory class used by clients to create
// instances of this service
static CnvFactory<AdfEventCnv> s_factory;
const ICnvFactory& AdfEventCnvFactory = s_factory;

AdfEventCnv::AdfEventCnv(ISvcLocator* svc)
: LdfBaseCnv(classID(), svc)
{
    // Here we associate this converter with the /Event path on the TDS.
    declareObject("/Event/AncillaryEvent/AdfEvent", objType(), "PASS");
}


StatusCode AdfEventCnv::createObj(IOpaqueAddress* , 
                               DataObject*& refpObject) {
    // Purpose and Method:  This converter will create an AdfEvent on
    //   the TDS.
    MsgStream log(msgSvc(), "AdfEventCnv");
    AncillaryData::AdfEvent *adf = new AncillaryData::AdfEvent();
    refpObject  = adf;
    // Retrieve the LAT data for this event 
    const ldfReader::AdfData& adfLdf = ldfReader::LatData::instance()->getAdf();
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

StatusCode AdfEventCnv::updateObj(int* , AncillaryData::AdfEvent*) {
    // Purpose and Method:  This method does nothing other than announce it has
    //   been called.

    MsgStream log(msgSvc(), "AdfEventCnv");
    log << MSG::DEBUG << "AdfEventCnv::updateObj" << endreq;
    return StatusCode::SUCCESS;
}

