
#include <reconRootData/GcrXtal.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(GcrXtal)

//-----------------------------------------------------------------------------------------------------------------
void GcrXtal::print() const
{
  std::cout << "---> writeOutPrint GcrXtal" << std::endl;
  std::cout << "---> xtalID =" << getXtalId() << std::endl;
  std::cout << "---> pathLength =" << getPathLength() << std::endl;
 
}

void GcrXtal::Print(Option_t* /*option*/) const
{
  std::cout << "---> writeOutPrint GcrXtal" << std::endl;
  std::cout << "---> xtalID =" << getXtalId() << std::endl;
  std::cout << "---> pathLength =" << getPathLength() << std::endl;

}

//===================================================
// dummy data, just for tests
//===================================================

void GcrXtal::Fake( Int_t /*ievent*/, Int_t ixtal, Float_t /*randNum*/) 
{
    setXtalId(CalXtalId(0,1,2)) ;
    setPathLength(1.*ixtal) ;
    setClosestFaceDist(2.*ixtal) ;
    setCrossedFaces(3*ixtal) ;
    setEntryPoint(TVector3(4.,5.,6.)) ;
    setExitPoint(TVector3(7.,5.,6.)) ;
}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t GcrXtal::CompareInRange( const GcrXtal & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(XtalId) && result ;
    result = COMPARE_IN_RANGE(PathLength) && result ;
    result = COMPARE_IN_RANGE(ClosestFaceDist) && result ;
    result = COMPARE_IN_RANGE(CrossedFaces) && result ;
    result = COMPARE_IN_RANGE(EntryPoint) && result ;
    result = COMPARE_IN_RANGE(ExitPoint) && result ;

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

