#include "reconRootData/AcdTkrLocalCoordsV2.h"
#include "Riostream.h"

#include "TMatrixD.h"

ClassImp(AcdTkrLocalCoordsV2)

AcdTkrLocalCoordsV2::AcdTkrLocalCoordsV2() 
  :m_volume(-1),
   m_arcLengthPlane(0.),
   m_cosTheta(0.),
   m_global(0.,0.,0.),
   m_localCovProj(2),
   m_localCovProp(2){
  m_local[0] = 0.;
  m_local[1] = 0.;  
  m_active[0] = m_active[1] = 0.;
}


AcdTkrLocalCoordsV2::AcdTkrLocalCoordsV2(Int_t volume, Float_t arcLength, Float_t cosTheta, 
                                     const TVector3& global, 
                                     const Float_t local[2],const Float_t active[2],
                                     const TMatrixDSym& localCovProj, const TMatrixDSym& localCovProp)
  :m_volume(volume),
   m_arcLengthPlane(arcLength),
   m_cosTheta(cosTheta),
   m_global(global),
   m_localCovProj(localCovProj),
   m_localCovProp(localCovProp){
  m_local[0] = local[0];
  m_local[1] = local[1];  
  m_active[0] = active[0];
  m_active[1] = active[1];
}

AcdTkrLocalCoordsV2::AcdTkrLocalCoordsV2(const AcdTkrLocalCoordsV2& other)
  :m_volume(other.m_volume),
   m_arcLengthPlane(other.m_arcLengthPlane),
   m_cosTheta(other.m_cosTheta),
   m_global(other.m_global),
   m_localCovProj(other.m_localCovProj),
   m_localCovProp(other.m_localCovProp){
  m_local[0] = other.m_local[0];
  m_local[1] = other.m_local[1];  
  m_active[0] = other.m_active[0];
  m_active[1] = other.m_active[1];  
}

void AcdTkrLocalCoordsV2::clear()
{
  m_volume = -1;
  m_arcLengthPlane = 0.;
  m_cosTheta = 0.;
  m_global.SetXYZ(0.,0.,0.);
  m_local[0] = 0.;
  m_local[1] = 0.;
  m_active[0] = m_active[1] = 0.;
  m_localCovProj.Zero();
  m_localCovProp.Zero();
}

void AcdTkrLocalCoordsV2::print(Option_t * /*option*/ ) const 
{
  using namespace std;
  cout << "AcdTkrLocalCoordsV2 vol: " << m_volume
       << ".  s " << m_arcLengthPlane
       << "  angle: " << m_cosTheta 
       << "   Global: <" << m_global.X() << ','  << m_global.Y() << ',' << m_global.Z()
       << ">.  Local:  {" << m_local[0] << ' ' << m_local[1] 
       << "}.  Active:  {" << m_active[0] << ' ' << m_active[1] 
       << "}  CovProj: " << std::endl;
  m_localCovProj.Print();
  cout << "}  CovProp: " << std::endl;
  m_localCovProp.Print();
}

AcdTkrLocalCoordsV2& AcdTkrLocalCoordsV2::operator=(const AcdTkrLocalCoordsV2& params)
{
  copy(params);
  return *this;
}

void AcdTkrLocalCoordsV2::set(Int_t volume, Float_t arcLength, Float_t cosTheta, 
                              const TVector3& global, 
                              const Float_t local[2], const Float_t active[2],
                              const TMatrixDSym& localCovProj, const TMatrixDSym& localCovProp)
{
  m_volume = volume;
  m_arcLengthPlane = arcLength;
  m_cosTheta = cosTheta;
  m_global = global;
  m_local[0] = local[0];
  m_local[1] = local[1];
  m_active[0] = active[0];
  m_active[1] = active[1];
  m_localCovProj = localCovProj;
  m_localCovProp = localCovProp;
}  

void AcdTkrLocalCoordsV2::copy(const AcdTkrLocalCoordsV2& other)
{
  set(other.m_volume,other.m_arcLengthPlane,other.m_cosTheta,
      other.m_global,
      other.m_local,other.m_active,
      other.m_localCovProj,other.m_localCovProp);
}  


//======================================================
// For Unit Tests
//======================================================

void AcdTkrLocalCoordsV2::Fake( Int_t /* ievent */, UInt_t /*rank*/, Float_t /*randNum*/ ) {

   // NOT TESTED YET

 }

Bool_t AcdTkrLocalCoordsV2::CompareInRange( const AcdTkrLocalCoordsV2 & /*ref*/, const std::string & name ) const {

    bool result = true ;

    // NOT TESTED YET
    //result = ... && result ;

    if (!result) {
        if ( name == "" ) {
            std::cout<<"Comparison ERROR for "<<Class_Name()<<std::endl ;
        }
        else {
            std::cout<<"Comparison ERROR for "<<name<<std::endl ;
        }
    }
    return result ;

}

