
#include <gcrSelectRootData/GcrSelectedXtal.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(GcrSelectedXtal)

//-----------------------------------------------------------------------------------------------------------------
void GcrSelectedXtal::print() const
{
  std::cout << "---> writeOutPrint GcrSelectedXtal" << std::endl;
  std::cout << "---> xtalID =" << getXtalId() << std::endl;
  std::cout << "---> pathLength =" << getPathLength() << std::endl;
 
}

void GcrSelectedXtal::Print(Option_t *opt) const 
{
  std::cout << "---> writeOutPrint GcrSelectedXtal" << std::endl;
  std::cout << "---> xtalID =" << getXtalId() << std::endl;
  std::cout << "---> pathLength =" << getPathLength() << std::endl;
}
 


//===================================================
// dummy data, just for tests
//===================================================

void GcrSelectedXtal::Fake( Int_t ievent, UInt_t ixtal, Float_t randNum) 
{
    setXtalId(CalXtalId(ixtal));
    setRawEnergy(ievent*randNum);
    setPathLength(ievent*randNum*2.);
    setCorrEnergy(ievent*randNum*3.);
    setSelectGrade(ievent);
    setClosestFaceDist(ievent*randNum*4.);
    setCrossedFaces(ievent*5);
    setEntryPoint(TVector3(ievent*randNum*ixtal, ievent*randNum*ixtal*6.,
                           ievent*randNum*ixtal*7.));
    setExitPoint(TVector3(ievent*randNum*ixtal*8.,ievent*randNum*9.,ievent*randNum*10.));
    

}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t GcrSelectedXtal::CompareInRange( const GcrSelectedXtal & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(RawEnergy) && result ;
    result = COMPARE_IN_RANGE(PathLength) && result ;
    result = COMPARE_IN_RANGE(CorrEnergy) && result ;
    result = COMPARE_IN_RANGE(SelectGrade) && result ;
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

