//									  
// The AcdId class contains a tile ID number.	
// AcdId is used in the AcdTile class.
//									  

// Apr 2000 Daniel Flath - Minor changes to function names, etc.
// Jan 1999 Daniel Flath - ROOT HTML Documentation added
// Dec 1999 Daniel Flath - Rewrite for GLAST
// Version 1.1 25 Oct 1999 R.Dubois Clone from LCD towerID

#include <commonRootData/idents/AcdGapId.h>
#include <commonRootData/RootDataUtil.h>
#include <iostream>

ClassImp(AcdGapId)

void AcdGapId::Clear(Option_t *) {
  m_val = 0;
}

void AcdGapId::Print(Option_t * /* option */) const {
  using namespace std;
  cout.precision(2);
  cout << "Id: " << m_val << endl;
}

Bool_t AcdGapId::CompareInRange( const AcdGapId & gap, const std::string & name ) const {

  bool result = true ;

  result = m_val == gap.asShort();
  if (!result) {
    std::cout<<"Comparison ERROR for "<<name<<std::endl ;
  }
  return result ;
}

void AcdGapId::Fake( Int_t /* ievent */, UInt_t /* rank */, Float_t /* randNum */) {
  return;
}
