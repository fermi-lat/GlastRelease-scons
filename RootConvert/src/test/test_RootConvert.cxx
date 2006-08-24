
#include <RootConvert/MonteCarlo/McPositionHitConvert.h>
#include <RootConvert/Recon/CalClusterConvert.h>
#include <RootConvert/Recon/CalEventEnergyConvert.h>
#include <RootConvert/Recon/CalXtalRecDataConvert.h>
#include <RootConvert/Recon/AcdReconConvert.h>
#include <RootConvert/Digi/OnboardFilterConvert.h>
#include <commonRootData/RootDataUtil.h>
#include "TRandom.h"
#include <iostream>

template <class TdsClass, class RootClass>
bool checkTopObject( Float_t randNum ) {
    TdsClass tdsObj ;
    RootClass rootObj1, rootObj2 ;
    Int_t ievent ;
    for ( ievent=0 ; ievent<2 ; ++ievent ) {
    
            rootObj1.Fake(ievent,randNum) ;
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

template <class TdsClass, class RootClass>
bool checkElements( Float_t randNum ) {
    TdsClass tdsObj ;
    RootClass rootObj1, rootObj2 ;
    Int_t ievent ;
    UInt_t rank ;
    for ( ievent=0 ; ievent<2 ; ++ievent ) {
        for ( rank=0 ; rank<10 ; ++rank ) {
    
            rootObj1.Fake(ievent,rank,randNum) ;
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
    result = result && checkElements<Event::McPositionHit,McPositionHit>(randNum) ;
    result = result && rootdatautil::CompareInRange(ROOT_NUMCALLAYERS,NUMCALLAYERS,"NUMCALLAYERS") ;
    result = result && checkElements<Event::CalCluster,CalCluster>(randNum) ;
    result = result && checkElements<Event::CalEventEnergy,CalEventEnergy>(randNum) ;
    result = result && checkElements<Event::CalXtalRecData,CalXtalRecData>(randNum) ;
    result = result && checkTopObject<Event::AcdRecon,AcdRecon>(randNum) ;
    
    result = result && checkTopObject<OnboardFilterTds::FilterStatus, FilterStatus>(randNum);

    if (result) {
      std::cout<<"RootConvert test succeeded"<<std::endl ;
      return 0 ;
    } else {
      std::cout<<"RootConvert test ERROR"<<std::endl ;
      return 1 ;
    }
}


