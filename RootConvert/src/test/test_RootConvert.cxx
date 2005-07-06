
#include "RootConvert/MonteCarlo/McPositionHitConvert.h"
#include "RootConvert/Recon/CalClusterConvert.h"
#include "TRandom.h"
#include <iostream>

template <class TdsClass, class RootClass>
bool check( Float_t randNum ) {
    TdsClass tdsObj ;
    RootClass rootObj1, rootObj2 ;
    UInt_t i ;
    for ( i=0 ; i<10 ; ++i ) {

        RootPersistence::random(rootObj1,i,randNum) ;
        RootPersistence::convert(rootObj1,tdsObj) ;
        RootPersistence::convert(tdsObj,rootObj2) ;
        if (!RootPersistence::compare(rootObj1,rootObj2)) {
            std::cout
              <<"RootConvert test FAILED for "
              <<rootObj1.ClassName()
              <<std::endl ;
            return false ;
        }
    }
   return true ;
}

int main( int /* argc */, char ** /* argv */ ) {

    bool result = true ;
    TRandom randGen ;
    Float_t randNum = randGen.Rndm() ;
    result = result && check<Event::McPositionHit,McPositionHit>(randNum) ;
    result = result && check<Event::CalCluster,CalCluster>(randNum) ;
    
    if (result) {
      std::cout<<"RootConvert test suceeded"<<std::endl ;
      return 0 ;
    } else {
      std::cout<<"RootConvert test FAILED"<<std::endl ;
      return 1 ;
    }
}


