
#include <RootConvert/Recon/TkrTruncationInfoConvert.h>
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

#include <Event/Recon/TkrRecon/TkrTruncatedPlane.h>
#include <TObject.h>

#include <iostream>

namespace RootPersistence {

void convert(const std::pair<Event::SortId, Event::TkrTruncatedPlane> & tdsTruncationPair , 
             TkrTruncationData & rootTruncationData)
{
  Event::SortId id =  tdsTruncationPair.first;
  Event::TkrTruncatedPlane plane = tdsTruncationPair.second;
  rootTruncationData.setStatus(plane.getStatus()); 
  rootTruncationData.setStripCount(plane.getStripCount());
  rootTruncationData.setStripNumber(plane.getStripNumber());
  rootTruncationData.setLocalX(plane.getLocalX());
  rootTruncationData.setPlaneZ(plane.getPlaneZ());   
  rootTruncationData.setTower(id.getTower());
  rootTruncationData.setTray(id.getTray());
  rootTruncationData.setFace(id.getFace());
  rootTruncationData.setView(id.getView());  
}


  void convert( const TObjArray &  truncationDataColRoot, 
                Event::TkrTruncationInfo & tdsTruncationInfo )
  {
    TIter truncationDataIter(&truncationDataColRoot);
    TkrTruncationData*  truncationDataRoot = 0;
    static int nRCTrunc = 0;
    static int nCCTrunc = 0;
    while ((truncationDataRoot = 
            (TkrTruncationData*)truncationDataIter.Next())!=0) {
      int status = truncationDataRoot->getStatus();
      if( (status & Event::TkrTruncatedPlane::RCSET)!=0 ) nRCTrunc++;
      if( (status & Event::TkrTruncatedPlane::CCSET)!=0 ) nCCTrunc++;
      
      Event::TkrTruncatedPlane truncatedPlaneTds 
        ( status,
          truncationDataRoot->getStripCount(),
          truncationDataRoot->getStripNumber(),
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
