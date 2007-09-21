/*
#include "reconRootData/TkrTruncationData.h"
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(TkrTruncationData)
*/
//void TkrTruncationData::Clear(Option_t* /* option */) 
/*
{
  m_status=0;
  m_stripCount.clear();
  m_stripNumber.clear();
  m_localX.clear();
  m_planeZ=0.;
  m_tower=0;
  m_tray=0;
  m_face=0;
  m_view=0;
}

void TkrTruncationData::Print( Option_t * option) const
{
  TObject::Print(option);
  using namespace std;
  
  cout << "status : " << m_status << endl;
  //cout << "stripCount : " << m_stripCount << endl;
  //cout << "stripNumber : " << m_stripNumber << endl;
  //cout << "localX : " << m_localX << endl;
  cout << "tower, tray, face, view " << m_tower  << ", "
       << m_tray << ", " << m_face << ", " << m_view 
       << endl;
}

void TkrTruncationData::Fake(Int_t ievent, UInt_t rank, Float_t randNum ) {
  m_status=ievent;
  m_planeZ=randNum;
  m_tower = 13;
  m_tray=10;
  m_face =0;
  m_view=1;
}

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t TkrTruncationData::CompareInRange( const TkrTruncationData & ref, const std::string & name ) const {

    bool result = true ;

    // the test comes before "&& result" because we want to enforce all tests
    result = COMPARE_IN_RANGE(Status) && result ;
    result = COMPARE_IN_RANGE(StripCount) && result ;
    result = COMPARE_IN_RANGE(StripNumber) && result ;
    //result = COMPARE_IN_RANGE(LocalX) && result ;
    result = COMPARE_IN_RANGE(PlaneZ) && result ;
    result = COMPARE_IN_RANGE(Tower) && result ;
    result = COMPARE_IN_RANGE(Tray) && result ;
    result = COMPARE_IN_RANGE(Face) && result ;
    result = COMPARE_IN_RANGE(View) && result ;

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

*/