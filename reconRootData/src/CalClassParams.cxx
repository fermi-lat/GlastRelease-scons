

//  Implementation file of Root clone for CalClassParams
//  
// Authors:
//
//    Luca Baldini
//    Johan Bregeon
//

#include <reconRootData/CalClassParams.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalClassParams)
  
Double_t CalClassParams::getProb(const std::string &className) const
{
  if (hasClass(className)){
    return m_probMap.find(className)->second;
  }
  else {
    return -1.;
  }
}

Bool_t CalClassParams::hasClass(const std::string &className) const
{
  return m_probMap.count(className);
}

void CalClassParams::Clear( Option_t * )
{
  m_probMap.clear();
  init();
}

void CalClassParams::init()
{
  setProb("gam", -1.);
}

void CalClassParams::Print( Option_t * ) const
{
  std::cout << "CalClassParams::Print() not implemented, yet." << std::endl;
}

// dummy data, just for tests
void CalClassParams::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */ )
 {
   init();
 }

Bool_t CalClassParams::CompareInRange( const CalClassParams & cp, const std::string & name ) const
{  
  Bool_t result = true ;
  // Not implemented, yet.
  std::cout << "CalClassParams::CompareInRange() not implemented, yet." << std::endl;
  return result;
}

