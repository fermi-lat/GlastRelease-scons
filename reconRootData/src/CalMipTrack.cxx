
#include <reconRootData/CalMipTrack.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(CalMipTrack)

//-----------------------------------------------------------------------------------------------------------------
void CalMipTrack::print() const
{
  std::cout << "---> writeOutPrint calMipTrack" << std::endl;
  std::cout << "---> point =" << getPoint().x() << " " << getPoint().y() << " " << getPoint().z() << std::endl;
  std::cout << "---> dir   =" << getDir().x()   << " " << getDir().y()   << " " << getDir().z()   << std::endl;
  std::cout << "---> d2C / d2Edge / calEdge / arcLen = " << m_d2C << " " << m_d2Edge << " " << m_calEdge << " "<< m_arcLen << std::endl;
  std::cout << "---> ecor / ecorRms / chi2 / erm = " << m_ecor << " " << m_ecorRms << " " << m_chi2 << " " << m_erm << std::endl;
 
}


//===================================================
// dummy data, just for tests
//===================================================

void CalMipTrack::Fake( Int_t /* ievent */, UInt_t /* rank*/ , Float_t /* randNum */ ) {
    setPoint(TVector3(1.0,1.5, 1.9)) ;
    setDir(TVector3(2.0, 2.5, 2.9)) ;
    setD2C(1.) ;
    setD2Edge(2.) ;
    setCalEdge(3) ;
    setArcLen(4.) ;
    setEcor(5.) ;
    setEcorRms(6.) ;
    setChi2(7.) ;
    setErm(8.) ;
 }

#define COMPARE_IN_RANGE(att) rootdatautil::CompareInRange(get ## att(),ref.get ## att(),#att)

Bool_t CalMipTrack::CompareInRange( const CalMipTrack & ref, const std::string & name ) const {

    bool result = true ;

    result = COMPARE_IN_RANGE(Point) && result ;
    result = COMPARE_IN_RANGE(Dir) && result ;
    result = COMPARE_IN_RANGE(D2C) && result ;
    result = COMPARE_IN_RANGE(D2Edge) && result ;
    result = COMPARE_IN_RANGE(CalEdge) && result ;
    result = COMPARE_IN_RANGE(ArcLen) && result ;
    result = COMPARE_IN_RANGE(Ecor) && result ;
    result = COMPARE_IN_RANGE(EcorRms) && result ;
    result = COMPARE_IN_RANGE(Chi2) && result ;
    result = COMPARE_IN_RANGE(Erm) && result ;

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

