
#include <reconRootData/GcrTrack.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(GcrTrack)

//-----------------------------------------------------------------------------------------------------------------
void GcrTrack::Print(Option_t *option) const
{
  std::cout << "---> writeOutPrint GcrTrack" << std::endl;
 
}


//===================================================
// dummy data, just for tests
//===================================================

void GcrTrack::Fake( Int_t ievent, Float_t randNum) {
    setDirection(TVector3(1.0,1.5, 1.9)) ;
    setDirError(TVector3(2.0, 2.5, 2.9)) ;
    setCalEntryPoint(TVector3(3.0, 4.0, 5.0));
    setCalExitPoint(TVector3(6.0, 7.0, 8.0));
 }

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t GcrTrack::CompareInRange( const GcrTrack & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(Direction) && result ;
    result = COMPARE_IN_RANGE(DirError) && result ;
    result = COMPARE_IN_RANGE(CalEntryPoint) && result ;
    result = COMPARE_IN_RANGE(CalExitPoint) && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<ClassName()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

