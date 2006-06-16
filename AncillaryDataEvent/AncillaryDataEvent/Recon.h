#include "TaggerCluster.h"
#include "TaggerHit.h"
#include <vector>
#include "GaudiKernel/DataObject.h"

static const CLID& CLID_AncillaryDataReconEvent = InterfaceID("AncillaryDataReconEvent", 1, 0);

namespace AncillaryData
{
  class Recon : public DataObject 
    {
    public:
      Recon();
      static const CLID& classID()       { return CLID_AncillaryDataReconEvent; };
    private:
      std::vector<TaggerCluster> TaggerClusterColl;
    } 
};
