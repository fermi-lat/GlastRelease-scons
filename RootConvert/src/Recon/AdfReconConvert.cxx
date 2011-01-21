#include "AncillaryDataEvent/Recon.h"
#include "AncillaryDataEvent/TaggerHit.h"
#include "AncillaryDataEvent/TaggerCluster.h"

#include "reconRootData/AdfRecon.h"

#include "TCollection.h"

#include <RootConvert/Utilities/Toolkit.h>
#include "RootConvert/Digi/AdfDigiConvert.h"


namespace RootPersistence {


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
          tagIt->getPulseHeight(), tagIt->getSigma(), 
          tagIt->getPedestalSubtract());
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
    rootObj.initEventNumber(tdsObj.getEventNumber());
    rootObj.initSpillNumber(tdsObj.getSpillNumber());
    rootObj.initNumHighestCluster(tdsObj.getNumberOfHigestClusters());
    unsigned int m;
    for(m=0; m<AncillaryData::N_MODULES; m++)
      {
        rootObj.initNumCluster(0,m,tdsObj.getNumberOfClusters(0,m));
        rootObj.initNumCluster(1,m,tdsObj.getNumberOfClusters(1,m));
     }


    Double_t x[AncillaryData::N_MODULES], y[AncillaryData::N_MODULES];
    Double_t z[AncillaryData::N_MODULES];
    unsigned int i;
    for (i=0;i<AncillaryData::N_MODULES;i++) {
        x[i] = tdsObj.getX(i);
        y[i] = tdsObj.getY(i);
        z[i] = tdsObj.getZ(i);
    }
    rootObj.initXyz(x,y,z,AncillaryData::N_MODULES);
    rootObj.init(tdsObj.getPX(), tdsObj.getPY(), tdsObj.getPZ());
    rootObj.initPhi(tdsObj.getPhiIn(), tdsObj.getPhiOut());
    rootObj.initThetaPhi(tdsObj.getTheta(), tdsObj.getDeltaPhi());
    rootObj.initEnergy(tdsObj.getReconstructedEnergy(), 
                       tdsObj.getCorrectedEnergy());
    rootObj.initError(tdsObj.getErrorErec(), tdsObj.getErrorEcorr());

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
                           qdcIt->getQdcModule(), qdcIt->getSigma(),
                           qdcIt->getPedestalSubtract());
    }

    const std::vector<AncillaryData::ScalerHit>& scalerHitColTds = tdsObj.getScalerHitCol();
         std::vector<AncillaryData::ScalerHit>::const_iterator scalerIt;
         for(scalerIt = scalerHitColTds.begin(); scalerIt != scalerHitColTds.end(); scalerIt++) {
              rootObj.addScalerHit(scalerIt->getScalerChannel(), scalerIt->getScalerValue());
          }


  };

  void convert( const reconRootData::AdfRecon& rootObj, AncillaryData::Recon& tdsObj) {

    tdsObj.setEventNumber(rootObj.getEventNumber());
    tdsObj.setSpillNumber(rootObj.getSpillNumber());
    tdsObj.setNumberOfHighestClusters(rootObj.getNumHighestCluster());
    unsigned int m;
    for(m=0; m<AncillaryData::N_MODULES; m++)
      {
        tdsObj.setNumberOfClusters(0,m,rootObj.getNumCluster(0,m));
        tdsObj.setNumberOfClusters(1,m,rootObj.getNumCluster(1,m));
      }

    tdsObj.setIntersection(rootObj.getPx(), rootObj.getPy(), rootObj.getPz());
    tdsObj.setPhi(rootObj.getPhiIn(), rootObj.getPhiOut());
    tdsObj.setThetaPhi(rootObj.getTheta(), rootObj.getDeltaPhi());
    tdsObj.setEnergy(rootObj.getReconstructedEnergy(), rootObj.getCorrectedEnergy());
    tdsObj.setEnergyError(rootObj.getErrorReconstructedEnergy(),
                          rootObj.getErrorCorrectedEnergy());
    double x[AncillaryData::N_MODULES], y[AncillaryData::N_MODULES];
    double z[AncillaryData::N_MODULES];
    unsigned int i;
    for (i=0;i<AncillaryData::N_MODULES;i++) {
        x[i] = rootObj.getX(i);
        y[i] = rootObj.getY(i);
        z[i] = rootObj.getZ(i);
    }
    tdsObj.setPos(x,y,z);

    const TObjArray* taggerClusterColRoot = rootObj.getTaggerClusterCol();
    TIter tagIt(taggerClusterColRoot);
    reconRootData::TaggerCluster *tagRoot;
    while ((tagRoot = (reconRootData::TaggerCluster*)tagIt.Next())) {
        AncillaryData::TaggerCluster tagTds;
        convert(tagRoot, tagTds);
        tdsObj.appendTaggerCluster(tagTds);
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
