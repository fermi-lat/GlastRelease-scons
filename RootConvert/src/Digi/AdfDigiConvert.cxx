#include "AncillaryDataEvent/Digi.h"
#include "AncillaryDataEvent/TaggerHit.h"
#include "AncillaryDataEvent/QdcHit.h"

#include "digiRootData/AdfDigi.h"

#include "TCollection.h"

#include <RootConvert/Utilities/Toolkit.h>


namespace RootPersistence {

  void convert( const AncillaryData::ScalerHit* tdsObj,
                commonRootData::ScalerHit& rootObj) {
       rootObj.initialize(tdsObj->getScalerChannel(), tdsObj->getScalerValue());
  };

  void convert( const commonRootData::ScalerHit* rootObj,
                AncillaryData::ScalerHit& tdsObj) {
    tdsObj.setScalerChannel(rootObj->getChannel());
    tdsObj.setScalerValue(rootObj->getValue());
  };

  void convert( const AncillaryData::QdcHit* tdsObj,
                commonRootData::QdcHit& rootObj) {
       rootObj.initialize(tdsObj->getQdcChannel(), tdsObj->getPulseHeight(), tdsObj->getQdcModule(), tdsObj->getSigma(), tdsObj->getPedestalSubtract());
  };

  void convert( const commonRootData::QdcHit* rootObj,
                AncillaryData::QdcHit& tdsObj) {
    if (rootObj->isPedestalSubtracted()) tdsObj.setPedestalSubtract();
    tdsObj.setPulseHeight(rootObj->getPulseHeight());
    tdsObj.setQdcChannel(rootObj->getChannel());
    tdsObj.setQdcModule(rootObj->getModule());
    tdsObj.setSigma(rootObj->getSigma());
  };

  void convert( const AncillaryData::TaggerHit* tdsObj,
                commonRootData::TaggerHit& rootObj) {
      rootObj.initialize(tdsObj->getModuleId(), tdsObj->getLayerId(), 
                        tdsObj->getStripId(), tdsObj->getPulseHeight(), 
                        tdsObj->getSigma(), tdsObj->getPedestalSubtract());
  };

  void convert( const commonRootData::TaggerHit* rootObj,
                AncillaryData::TaggerHit& tdsObj) {
       tdsObj.setLayerId(rootObj->getLayerId());
           tdsObj.setModuleId(rootObj->getModuleId());
           tdsObj.setPulseHeight(rootObj->getPulseHeight());
           if (rootObj->isPedestalSubtracted()) tdsObj.setPedestalSubtract();
           tdsObj.setStripId(rootObj->getStripId());
           tdsObj.setSigma(rootObj->getSigma());
  };

  void convert( const AncillaryData::Digi& tdsObj, AdfDigi& rootObj) {
      rootObj.setEventNumber(tdsObj.getEventNumber());
          rootObj.setSpillNumber(tdsObj.getSpillNumber());

          const std::vector<AncillaryData::TaggerHit>& taggerHitColTds = tdsObj.getTaggerHitCol();
          std::vector<AncillaryData::TaggerHit>::const_iterator tagIt;
          for (tagIt = taggerHitColTds.begin(); tagIt != taggerHitColTds.end(); tagIt++) {
                  rootObj.addTaggerHit(tagIt->getModuleId(), 
                          tagIt->getLayerId(), tagIt->getStripId(), 
                          tagIt->getPulseHeight(), tagIt->getSigma(),
                          tagIt->getPedestalSubtract());
          }

          const std::vector<AncillaryData::ScalerHit>& scalerHitColTds = tdsObj.getScalerHitCol();
          std::vector<AncillaryData::ScalerHit>::const_iterator scalerIt;
          for(scalerIt = scalerHitColTds.begin(); scalerIt != scalerHitColTds.end(); scalerIt++) {
                  rootObj.addScalerHit(scalerIt->getScalerChannel(), scalerIt->getScalerValue());
          }

          const std::vector<AncillaryData::QdcHit>& qdcHitColTds = tdsObj.getQdcHitCol();
          std::vector<AncillaryData::QdcHit>::const_iterator qdcIt;
          for(qdcIt = qdcHitColTds.begin(); qdcIt != qdcHitColTds.end(); qdcIt++) {
                  rootObj.addQdcHit(qdcIt->getQdcChannel(), qdcIt->getPulseHeight(), qdcIt->getQdcModule(), qdcIt->getSigma(), qdcIt->getPedestalSubtract());
          }
  };

  void convert( const AdfDigi& rootObj, AncillaryData::Digi& tdsObj) {

          tdsObj.setEventNumber(rootObj.getEventNumber());
          tdsObj.setSpillNumber(rootObj.getSpillNumber());

          const TClonesArray* taggerHitColRoot = rootObj.getTaggerHitCol();
          TIter tagIt(taggerHitColRoot);
          commonRootData::TaggerHit *tagRoot;
          while ((tagRoot = (commonRootData::TaggerHit*)tagIt.Next())) {
                  AncillaryData::TaggerHit tagTds;
                  convert(tagRoot, tagTds);
                  tdsObj.appendTaggerHit(tagTds);
          }

          const TClonesArray* scalerHitColRoot = rootObj.getScalerHitCol();
          TIter scalerIt(scalerHitColRoot);
          commonRootData::ScalerHit *scalerRoot;
          while ((scalerRoot = (commonRootData::ScalerHit*)scalerIt.Next())) {
                  AncillaryData::ScalerHit scalerTds;
                  convert(scalerRoot, scalerTds);
                  tdsObj.appendScalerHit(scalerTds);
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
