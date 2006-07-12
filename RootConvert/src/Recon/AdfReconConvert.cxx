#include "AncillaryDataEvent/Recon.h"
#include "AncillaryDataEvent/TaggerHit.h"
#include "AncillaryDataEvent/TaggerCluster.h"

#include "reconRootData/AdfRecon.h"

#include "TCollection.h"

#include <RootConvert/Utilities/Toolkit.h>
#include "RootConvert/Digi/AdfDigiConvert.h"


namespace RootPersistence {


/*
  void convert( const AncillaryData::TaggerHit* tdsObj,
                commonRootData::TaggerHit& rootObj) {
      rootObj.initialize(tdsObj->getModuleId(), tdsObj->getLayerId(), tdsObj->getStripId(),
            tdsObj->getPulseHeight(), tdsObj->getPedestalSubtract());
  };

  void convert( const commonRootData::TaggerHit* rootObj,
                AncillaryData::TaggerHit& tdsObj) {
       tdsObj.setLayerId(rootObj->getLayerId());
       tdsObj.setModuleId(rootObj->getModuleId());
       tdsObj.setPulseHeight(rootObj->getPulseHeight());
       if (rootObj->isPedestalSubtracted()) tdsObj.setPedestalSubtract();
       tdsObj.setStripId(rootObj->getStripId());
  };
*/

  void convert( const AncillaryData::TaggerCluster& tdsObj,
                reconRootData::TaggerCluster*& rootObj) {

    rootObj->initialize(tdsObj.getPosition(), tdsObj.getPulseHeight(),
                      tdsObj.getEta());

    AncillaryData::TaggerHit tagHitTds = tdsObj.getHighestHit();
    commonRootData::TaggerHit tagHitRoot;
    convert(&tagHitTds, tagHitRoot);

    rootObj->setHighestHit(tagHitRoot);
    
    const std::vector<AncillaryData::TaggerHit>& taggerHitColTds = tdsObj.getHits();
    std::vector<AncillaryData::TaggerHit>::const_iterator tagIt;
    for (tagIt = taggerHitColTds.begin(); tagIt != taggerHitColTds.end(); tagIt++) {
    rootObj->addTaggerHit(tagIt->getModuleId(), tagIt->getLayerId(), 
          tagIt->getStripId(), 
	  tagIt->getPulseHeight(), tagIt->getPedestalSubtract());
    }
  };

  void convert( const reconRootData::TaggerCluster* rootObj,
                AncillaryData::TaggerCluster& tdsObj) {

    const commonRootData::TaggerHit& tagHitRoot = rootObj->getHighestHit();
    AncillaryData::TaggerHit tagHitTds;
    convert(&tagHitRoot, tagHitTds);

    const TClonesArray* taggerHitColRoot = rootObj->getTaggerHitCol();
    TIter tagIt(taggerHitColRoot);
    commonRootData::TaggerHit* tagRoot;
    while ((tagRoot = (commonRootData::TaggerHit*)tagIt.Next())) {
	  AncillaryData::TaggerHit tagTds;
	  convert(tagRoot, tagTds);
	  tdsObj.append(tagTds);
    }

    tdsObj.initialize(rootObj->getPosition(), rootObj->getPulseHeight(),
                      rootObj->getEta(), tagHitTds);
  };

  void convert( const AncillaryData::Recon& tdsObj, reconRootData::AdfRecon& rootObj) {
    rootObj.setEventNumber(tdsObj.getEventNumber());
    rootObj.setSpillNumber(tdsObj.getSpillNumber());

    const std::vector<AncillaryData::TaggerCluster>& taggerClusterColTds = tdsObj.getTaggerClusters();
    std::vector<AncillaryData::TaggerCluster>::const_iterator tagIt;
    for (tagIt = taggerClusterColTds.begin(); tagIt != taggerClusterColTds.end(); tagIt++) {
        reconRootData::TaggerCluster *clusterRoot = new reconRootData::TaggerCluster();
        convert(*tagIt, clusterRoot);
        rootObj.addTaggerCluster(clusterRoot); 
    }

    const std::vector<AncillaryData::QdcHit>& qdcHitColTds = tdsObj.getQdcHitCol();
    std::vector<AncillaryData::QdcHit>::const_iterator qdcIt;
    for(qdcIt = qdcHitColTds.begin(); qdcIt != qdcHitColTds.end(); qdcIt++) {
         rootObj.addQdcHit(qdcIt->getQdcChannel(), qdcIt->getPulseHeight(), 
                           qdcIt->getQdcModule(),
                           qdcIt->getPedestalSubtract());
    }

  };

  void convert( const reconRootData::AdfRecon& rootObj, AncillaryData::Recon& tdsObj) {

	  tdsObj.setEventNumber(rootObj.getEventNumber());
	  tdsObj.setSpillNumber(rootObj.getSpillNumber());

	  const TObjArray* taggerClusterColRoot = rootObj.getTaggerClusterCol();
	  TIter tagIt(taggerClusterColRoot);
	  reconRootData::TaggerCluster *tagRoot;
	  while ((tagRoot = (reconRootData::TaggerCluster*)tagIt.Next())) {
		  AncillaryData::TaggerCluster tagTds;
		  convert(tagRoot, tagTds);
		  tdsObj.appendTaggerCluster(tagTds);
	  }


          const TClonesArray* qdcHitColRoot = rootObj.getQdcHitCol();
          TIter qdcIt(qdcHitColRoot);
          commonRootData::QdcHit *qdcRoot;
          while ((qdcRoot = (commonRootData::QdcHit*)qdcIt.Next())) {
                  AncillaryData::QdcHit qdcTds;
                  convert(qdcRoot, qdcTds);
                  tdsObj.appendQdcHit(qdcTds);
          }
  };

 
}
