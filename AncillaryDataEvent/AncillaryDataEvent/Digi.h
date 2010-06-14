#include "TaggerHit.h"
#include <vector>
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"


#include "AdfEvent/AdfEvent.h"
#include "AncillaryDataEvent/TaggerHit.h"
#include "AncillaryDataEvent/QdcHit.h"

static const CLID& CLID_AncillaryDataDigiEvent = InterfaceID("AncillaryDataDigiEvent", 1, 0);

namespace AncillaryData
{
  class Digi : public DataObject 
    {
    public:
      Digi(){;}
      Digi(AncillaryData::AdfEvent *adfEvent);
      
      static const CLID& classID()       { return CLID_AncillaryDataDigiEvent; };
      void appendTaggerHit(TaggerHit h){TaggerHitColl.push_back(h);}
      void appendQdcHit(QdcHit h){QdcHitColl.push_back(h);}
      
      void print();
      
    private:
      std::vector<TaggerHit> TaggerHitColl;
      std::vector<QdcHit> QdcHitColl;
      
      unsigned m_eventNumber;
      unsigned m_spillNumber;
    };
};
