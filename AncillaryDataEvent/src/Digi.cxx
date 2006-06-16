#include "AncillaryDataEvent/Digi.h"
#include "AdfEvent/AdfEvent.h"
namespace AncillaryData
{
  Digi::Digi(AncillaryData::AdfEvent *adfEvent)
  {
    std::vector<QdcDataWord>  qdcDataColl   = adfEvent->getQdcDataWord();
    std::vector<FadcDataWord>  fadcDataColl = adfEvent->getFadcDataWord();
    std::vector<QdcDataWord>::iterator qdcIterator;
    std::vector<FadcDataWord>::iterator fadcIterator;
    
    for (qdcIterator = qdcDataColl.begin(); qdcIterator!=qdcDataColl.end(); ++qdcIterator)
      {
	appendQdcHit(QdcHit(*qdcIterator));
      }
    for (fadcIterator = fadcDataColl.begin(); fadcIterator!=fadcDataColl.end(); ++fadcIterator)
      {
	appendTaggerHit(TaggerHit(*fadcIterator));
      }
  }
  
  void Digi::print()
  {
    for(std::vector<TaggerHit>::iterator pos=TaggerHitColl.begin();pos!=TaggerHitColl.end(); ++pos)
      (*pos).print();
    
    for(std::vector<QdcHit>::iterator pos=QdcHitColl.begin();pos!=QdcHitColl.end(); ++pos)
      (*pos).print();
  }
}
