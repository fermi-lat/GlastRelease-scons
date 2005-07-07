
#include <RootConvert/MonteCarlo/McPositionHitConvert.h>
#include <RootConvert/Recon/CalClusterConvert.h>
#include <RootConvert/Recon/CalEventEnergyConvert.h>
#include <commonRootData/RootDataUtil.h>
#include "TRandom.h"
#include <iostream>

template <class TdsClass, class RootClass>
bool check( Float_t randNum ) {
    TdsClass tdsObj ;
    RootClass rootObj1, rootObj2 ;
    UInt_t i ;
    for ( i=0 ; i<10 ; ++i ) {

        rootObj1.Fake(i,randNum) ;
        RootPersistence::convert(rootObj1,tdsObj) ;
        RootPersistence::convert(tdsObj,rootObj2) ;
        Bool_t theSame = rootObj1.CompareInRange(rootObj2) ;
        if (!theSame) {
            std::cout
              <<"RootConvert test FAILED for "
              <<rootObj1.ClassName()
              <<std::endl ;
            return false ;
        }
    }
   std::cout
     <<"RootConvert test succeeded for "
     <<rootObj1.ClassName()
     <<std::endl ;
   return true ;
}

int main( int /* argc */, char ** /* argv */ ) {

    bool result = true ;
    TRandom randGen ;
    Float_t randNum = randGen.Rndm() ;
    result = result && check<Event::McPositionHit,McPositionHit>(randNum) ;
    result = result && rootdatautil::CompareInRange(ROOT_NUMCALLAYERS,NUMCALLAYERS,"NUMCALLAYERS") ;
    result = result && check<Event::CalCluster,CalCluster>(randNum) ;
    result = result && check<Event::CalEventEnergy,CalEventEnergy>(randNum) ;
    
    if (result) {
      std::cout<<"RootConvert test suceeded"<<std::endl ;
      return 0 ;
    } else {
      std::cout<<"RootConvert test ERROR"<<std::endl ;
      return 1 ;
    }
}


