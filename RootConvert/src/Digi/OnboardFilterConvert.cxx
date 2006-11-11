#include "RootConvert/Digi/OnboardFilterConvert.h"

#include <RootConvert/Utilities/Toolkit.h>


namespace RootPersistence {

  void convert(const TFC_projection &tdsObj, TfcProjection &rootObj) {
      rootObj.initialize(tdsObj.intercept, tdsObj.slope, tdsObj.acdTopMask,
              tdsObj.acdXMask, tdsObj.acdYMask, tdsObj.layers, 
              tdsObj.hits, tdsObj.skirtMask, tdsObj.min, tdsObj.max,
              tdsObj.nhits);
  }

  void convert(const TfcProjection &rootObj, TFC_projection &tdsObj) {
      tdsObj.intercept = rootObj.getIntercept();
      tdsObj.slope = rootObj.getSlope(); 
      tdsObj.acdTopMask = rootObj.getAcdTopMask();
      tdsObj.acdXMask = rootObj.getAcdXMask();
      tdsObj.acdYMask = rootObj.getAcdYMask();
      tdsObj.skirtMask = rootObj.getSkirtMask();
      tdsObj.min = rootObj.getMin();
      tdsObj.max = rootObj.getMax();
      tdsObj.nhits = rootObj.getNhits();
      tdsObj.layers = rootObj.getLayers(); 
      const Short_t *hitsRoot = rootObj.getHits();
      unsigned int i;
      for (i=0;i<18;i++) 
          tdsObj.hits[i] = hitsRoot[i]; 
  }

  void convert(const TFC_projectionDir &tdsObj, TfcProjectionDir& rootObj) {
      rootObj.initialize(tdsObj.idx, tdsObj.xCnt, tdsObj.yCnt);
  }

  void convert(const TfcProjectionDir &rootObj, TFC_projectionDir& tdsObj) {
      tdsObj.idx = rootObj.getIndex();
      tdsObj.xCnt = rootObj.getXCount();
      tdsObj.yCnt = rootObj.getYCount();
  }

  void convert( const TFC_projections* tdsObj, TfcProjectionCol& rootObj) {

      int i;
      TfcProjectionDir dirRoot[16];
      for (i=0;  i<16; i++) 
          convert(tdsObj->dir[i], dirRoot[i]);

      TfcProjection *projectColRoot = 0;
      projectColRoot = new TfcProjection[tdsObj->curCnt];

      for (i=0; i<tdsObj->curCnt; i++) 
          convert(tdsObj->prjs[i], projectColRoot[i]);

      rootObj.initialize(tdsObj->maxCnt, tdsObj->curCnt, tdsObj->twrMsk, 
          dirRoot, projectColRoot);

      delete [] projectColRoot;
  }

  void convert ( const TfcProjectionCol& rootObj, TFC_projections*& tdsObj ) {
    tdsObj->maxCnt = rootObj.getMaxCount();
    tdsObj->curCnt = rootObj.getCurrentCount();
    tdsObj->twrMsk = rootObj.getTwrMsk();
    const TfcProjectionDir *dirRoot = rootObj.getProjectionDirCol();
    UInt_t i;
    for (i=0;i<16;i++) 
        convert(dirRoot[i], tdsObj->dir[i]);
 
    const TfcProjection* projectionArrRoot = rootObj.getProjectionCol();
    // Check to be sure the number of elements in the array doesn't exceed 
    // 1000, that's the setting for this fixed array on the TDS
    unsigned int numProj = rootObj.getCurrentCount();
    if (rootObj.getCurrentCount() >= 1000) 
        numProj=1000;
    for (i=0; i<numProj; i++) 
       convert(projectionArrRoot[i], tdsObj->prjs[i]); 

  }

  void convert( const std::vector<OnboardFilterTds::track>& tdsObj, 
                TObjArray& rootObj) {
      std::vector<OnboardFilterTds::track>::const_iterator trackIt;
      for(trackIt = tdsObj.begin(); trackIt != tdsObj.end(); trackIt++) {
          ObfTrack *rootTrack = new ObfTrack();
          rootTrack->initialize(trackIt->numhits, trackIt->phi_rad,
                     trackIt->theta_rad, trackIt->length, trackIt->lowCoord, 
                     trackIt->highCoord, trackIt->exLowCoord, 
                     trackIt->exHighCoord);
          rootObj.Add(rootTrack);
      }
  }

  void convert ( const TObjArray* rootObj, 
                 std::vector<OnboardFilterTds::track>& tdsObj) {


      TIter rootIt(rootObj);
      ObfTrack* rootTrack;
      while((rootTrack = (ObfTrack*)rootIt.Next())) {
          OnboardFilterTds::track curTrack;
          curTrack.numhits = rootTrack->getNumHits();
          curTrack.phi_rad = rootTrack->getPhi();
          curTrack.theta_rad = rootTrack->getTheta();
          curTrack.lowCoord = rootTrack->getLowCoord();
          curTrack.highCoord = rootTrack->getHighCoord();
          curTrack.exLowCoord = rootTrack->getExLowCoord();
          curTrack.exHighCoord = rootTrack->getExHighCoord();
          curTrack.length = rootTrack->getLength();
          tdsObj.push_back(curTrack);
      }
  }


  void convert( OnboardFilterTds::FilterStatus& tdsObj, 
                FilterStatus& rootObj) {

    Int_t ebfSize = 0;
    rootObj.init(tdsObj.get(),tdsObj.getStageEnergy(),tdsObj.getTcids(),
                 ebfSize);
    rootObj.initGem(tdsObj.getGemThrTkr(), tdsObj.getGemCalHiLo(),
            tdsObj.getGemCondsumCno(),
            tdsObj.getGemAcd_vetoes_XZ(),tdsObj.getGemAcd_vetoes_YZ(),
            tdsObj.getGemAcd_vetoes_XY(), tdsObj.getGemAcd_vetoes_RU(),
            tdsObj.getGemLivetime(), tdsObj.getGemTrgtime(), 
            tdsObj.getGemPpstime(), tdsObj.getGemDiscarded(),
            tdsObj.getGemPrescaled(), tdsObj.getGemSent());

    int xz, yz, xy;
    tdsObj.getAcdMap(xz, yz, xy);
    const Int_t *acdStatus = tdsObj.getAcdStatus();

    rootObj.initAcd(xz, yz, xy, tdsObj.getVetoWord(),acdStatus);

    rootObj.initTkr(tdsObj.getXcapture(), tdsObj.getYcapture(),
                    tdsObj.getXY00(), tdsObj.getXY11(), tdsObj.getXY22(),
                    tdsObj.getXY33(), tdsObj.getTmsk());

    rootObj.initLayers(tdsObj.getLayers());

    Int_t numLogsHit = 0; //  can't find this right now
    rootObj.initCal(numLogsHit, tdsObj.getLayerEnergy());

    rootObj.initSeparation(tdsObj.getSeparation());
    
    int xHits, yHits;
    double slopeXZ,slopeYZ,intXZ,intYZ;
    tdsObj.getBestTrack(xHits,yHits,slopeXZ,slopeYZ,intXZ,intYZ);
    rootObj.initBestTrack(xHits,yHits,slopeXZ,slopeYZ,intXZ,intYZ);

    TObjArray trackRootCol;
    if (tdsObj.tracksExist()) {
        const std::vector<OnboardFilterTds::track> trackTdsCol = tdsObj.getTracks();
        convert(trackTdsCol, trackRootCol);
    }
    rootObj.initTracks(trackRootCol);

    TfcProjectionCol projRootCol;
    const TFC_projections* projTds = tdsObj.getProjections();
    convert(projTds, projRootCol);
    rootObj.initProjectionCol(projRootCol);

  
  };

  void convert( const FilterStatus& rootObj, 
                OnboardFilterTds::FilterStatus& tdsObj) {

      tdsObj.set(rootObj.getStatus());
      tdsObj.setStageEnergy(rootObj.getStageEnergy());
      tdsObj.setTcids(rootObj.getTcids());

      tdsObj.setGemThrTkr(rootObj.getGemThrTkr());
      tdsObj.setGemCalHiLo(rootObj.getGemCalHiLo());
      tdsObj.setGemCondsumCno(rootObj.getGemCondsumCno());
      tdsObj.setGemAcd_vetoes_XZ(rootObj.getGemAcd_vetoes_XZ());
      tdsObj.setGemAcd_vetoes_YZ(rootObj.getGemAcd_vetoes_YZ());
      tdsObj.setGemAcd_vetoes_XY(rootObj.getGemAcd_vetoes_XY());
      tdsObj.setGemAcd_vetoes_RU(rootObj.getGemAcd_vetoes_RU());
      tdsObj.setGemLivetime(rootObj.getGemLivetime());
      tdsObj.setGemTrgtime(rootObj.getGemTrgtime());
      tdsObj.setGemPpstime(rootObj.getGemPpstime());
      tdsObj.setGemDiscarded(rootObj.getGemDiscarded());
      tdsObj.setGemPrescaled(rootObj.getGemPrescaled());
      tdsObj.setGemSent(rootObj.getGemSent());
      tdsObj.setBestTrack(rootObj.getXhitsBestTrack(),rootObj.getYhitsBestTrack(),
          rootObj.getSlopeXzBestTrack(), rootObj.getSlopeYzBestTrack(),rootObj.getIntXzBestTrack(),
          rootObj.getIntYzBestTrack());

      Int_t xz, yz, xy;
      rootObj.getAcdMap(xz,yz,xy);
      tdsObj.setAcdMap(xz,yz,xy);
      
      Int_t iAcd;
      const Int_t* statusRoot = rootObj.getAcdStatus();
      for (iAcd=0; iAcd < 16; iAcd++) {
          tdsObj.setAcdStatus(iAcd, statusRoot[iAcd]);
      }

      tdsObj.setVetoWord(rootObj.getVetoWord());
    
      tdsObj.setLayers(rootObj.getLayerCol());

      tdsObj.setSeparation(rootObj.getSeparation());

      //Ask Richard to provide suitable set method using floats rather than int
      //tdsObj->setLayerEnergy(rootObj.getLayerEnergy());

      tdsObj.setCapture(rootObj.getXcapture(), rootObj.getYcapture());

      tdsObj.setXY(rootObj.getXY00(), rootObj.getXY11(), rootObj.getXY22(),
                   rootObj.getXY33());

      tdsObj.setTmsk(rootObj.getTmsk());



  };


}
