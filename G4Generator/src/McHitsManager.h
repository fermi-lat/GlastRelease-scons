#ifndef MCHITSMANAGER_H
#define MCHITSMANAGER_H

#include <algorithm>

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/MsgStream.h"

// GlastEvent for creating the McEvent stuff
#include "GlastEvent/MonteCarlo/McPositionHit.h"
//#include "GaudiKernel/IDataProviderSvc.h"
#include <vector>

/** 
 *  @class McHitsManager
 *
 *  @brief 
 *
 *  @author 
 */

class CompareHits {
  public:
  bool operator()(mc::McPositionHit *left, mc::McPositionHit *right)
    {return left->volumeID() < right->volumeID();}

    };

class McHitsManager {
 public:

  static McHitsManager* getPointer();

  void addHit(mc::McPositionHit *hit);

  McPositionHitVector* getVector(){
    return m_posHit;
  }

  void init(){m_posHit = new McPositionHitVector;}
  void sort(){std::sort(m_posHit->begin(),m_posHit->end(), CompareHits());}

 private:
  McHitsManager(){}; 

  static McHitsManager* pointer;
  
  McPositionHitVector *m_posHit;  
};

#endif //McHitsManager_H
