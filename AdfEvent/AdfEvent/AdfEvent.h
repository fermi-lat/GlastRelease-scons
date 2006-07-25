#ifndef ADFEVENT_HH
#define ADFEVENT_HH

#include <vector>
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/IInterface.h"

#include "AdfEvent/EventSummaryData.h"
#include "AdfEvent/QdcDataWord.h"
#include "AdfEvent/FadcDataWord.h"
#include "AdfEvent/ScalerDataWord.h"
#include "AdfEvent/QdcHeaderWord.h"
#include "AdfEvent/FadcHeaderWord.h"
#include "AdfEvent/ScalerHeaderWord.h"

static const CLID& CLID_AncillaryEvent = InterfaceID("AncillaryEvent", 1, 0);
static const CLID& CLID_AncillaryDataFormat = InterfaceID("AncillaryDataFormat", 1, 0);

namespace AncillaryData
{
  class AdfEvent : public DataObject 
    {
    public:
      AdfEvent(){;}
      static const CLID& classID()       { return CLID_AncillaryDataFormat; };
      // Event Summary:
      void setEventSummaryData(EventSummaryData eventSummaryData) {m_eventSummaryData=eventSummaryData;}
      EventSummaryData getEventSummaryData() const {return m_eventSummaryData;}
      
      void appendQdcDataWord(QdcDataWord word){QdcData.push_back(word);}
      void appendFadcDataWord(FadcDataWord word){FadcData.push_back(word);}
      void appendScalerDataWord(ScalerDataWord word){ScalerData.push_back(word);}
      
      std::vector<QdcDataWord> getQdcDataWord(){return QdcData;}
      std::vector<FadcDataWord> getFadcDataWord(){return FadcData;}
      std::vector<ScalerDataWord> getScalerDataWord(){return ScalerData;}
      void print()
	{
	  m_eventSummaryData.printEventHeader();
	  std::cout<<" QDC:    Size: "<<QdcData.size()<<std::endl;
	  std::cout<<" FADC:   Size: "<<FadcData.size()<<std::endl;
	  std::cout<<" SCALER: Size: "<<ScalerData.size()<<std::endl;
	}
      void dump()
	{
	  
	  for (std::vector<QdcDataWord>::iterator qdcIterator=QdcData.begin(); qdcIterator!=QdcData.end(); qdcIterator++)
	    {
	      std::cout<<(*qdcIterator).getQdcChannel()<<" "<<(*qdcIterator).getQdcValue()<<std::endl;
	    }
	  
	  for (std::vector<FadcDataWord>::iterator fadcIterator=FadcData.begin(); fadcIterator!=FadcData.end(); fadcIterator++)
	    {
	      std::cout<<(*fadcIterator).getFadcModule();
	      std::cout<<" "<<(*fadcIterator).getFadcChannel();
	      std::cout<<" "<<(*fadcIterator).getFadcValue()<<std::endl;
	    }
	  
	  for (std::vector<ScalerDataWord>::iterator scalerIterator=ScalerData.begin(); scalerIterator!=ScalerData.end(); scalerIterator++)
	    {
	      std::cout<<(*scalerIterator).getScalerValue()<<std::endl;
	    }
	  //	  std::vector
	}
      // 
    private:
      EventSummaryData m_eventSummaryData;      
      std::vector<QdcDataWord> QdcData;
      std::vector<FadcDataWord> FadcData;
      std::vector<ScalerDataWord> ScalerData;
    };
}
#endif
