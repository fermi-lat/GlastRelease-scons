
#include <RootConvert/Recon/TkrTruncationInfoConvert.h>
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

#include <Event/Recon/TkrRecon/TkrTruncatedPlane.h>
#include <TObject.h>
#include <TCollection.h>

#include <iostream>

namespace RootPersistence {

void convert(const std::pair<Event::SortId, Event::TkrTruncatedPlane> & tdsTruncationPair , 
             TkrTruncationData & rootTruncationData)
{
  Event::SortId id =  tdsTruncationPair.first;
  Event::TkrTruncatedPlane plane = tdsTruncationPair.second;
  rootTruncationData.setStatus(plane.getStatus()); 
  //hack so I can get what I want to do done! 
  const Event::intVector& tdsCntVec = plane.getStripCount();
  std::vector<Int_t> rootVec; //intVector rootVec;
  for(Event::intVector::const_iterator tdsVecIter = tdsCntVec.begin(); tdsVecIter != tdsCntVec.end(); tdsVecIter++)
  {
      rootVec.push_back(*tdsVecIter);
  }
  rootTruncationData.setStripCount(rootVec);
  const Event::intVector& tdsNumVec = plane.getStripNumber();
  rootVec.clear();
  for(Event::intVector::const_iterator tdsVecIter = tdsNumVec.begin(); tdsVecIter != tdsNumVec.end(); tdsVecIter++)
  {
      rootVec.push_back(*tdsVecIter);
  }
  rootTruncationData.setStripNumber(rootVec);
  rootTruncationData.setLocalX(plane.getLocalX());
  rootTruncationData.setPlaneZ(plane.getPlaneZ());   
  rootTruncationData.setTower(id.getTower());
  rootTruncationData.setTray(id.getTray());
  rootTruncationData.setFace(id.getFace());
  rootTruncationData.setView(id.getView());  
}


  void convert( const TCollection &  truncationDataColRoot, 
                Event::TkrTruncationInfo & tdsTruncationInfo )
  {
    TIter truncationDataIter(&truncationDataColRoot);
    TkrTruncationData*  truncationDataRoot = 0;
    int nRCTrunc = 0;
    int nCCTrunc = 0;
    while ((truncationDataRoot = 
            (TkrTruncationData*)truncationDataIter.Next())!=0) {
      int status = truncationDataRoot->getStatus();
      if( (status & Event::TkrTruncatedPlane::RCSET)!=0 ) nRCTrunc++;
      if( (status & Event::TkrTruncatedPlane::CCSET)!=0 ) nCCTrunc++;

      //hack
      const std::vector<Int_t>& rootCntVec = truncationDataRoot->getStripCount();
      Event::intVector tdsCntVec;

      for(std::vector<Int_t>::const_iterator rootVecIter = rootCntVec.begin(); rootVecIter != rootCntVec.end(); rootVecIter++)
      {
          tdsCntVec.push_back(*rootVecIter);
      }
      
      const std::vector<Int_t>& rootNumVec = truncationDataRoot->getStripNumber();
      Event::intVector tdsNumVec;

      for(std::vector<Int_t>::const_iterator rootVecIter = rootNumVec.begin(); rootVecIter != rootNumVec.end(); rootVecIter++)
      {
          tdsCntVec.push_back(*rootVecIter);
      }
      
      Event::TkrTruncatedPlane truncatedPlaneTds 
        ( status,
          tdsCntVec,
          tdsNumVec,
          truncationDataRoot->getLocalX(),
          truncationDataRoot->getPlaneZ());
      
      tdsTruncationInfo.addItem(
                                truncationDataRoot->getTower(), 
                                truncationDataRoot->getTray(), 
                                truncationDataRoot->getFace(), 
                                truncationDataRoot->getView(), 
                                truncatedPlaneTds);
    }
    tdsTruncationInfo.setCCTrunc(nCCTrunc);
    tdsTruncationInfo.setRCTrunc(nRCTrunc);
    // mark it done!
   tdsTruncationInfo.setDone();    
  }
}
