#include "reconRootData/AcdHit.h"
#include "Riostream.h"

ClassImp(AcdHit)

AcdHit::AcdHit() {
  Clear();
}

void AcdHit::Print(Option_t* option) const {
  using namespace std;
  cout << "AcdHit.  id: " << m_acdId.getId() << "  phas: {" << m_pha[0] << ',' << m_pha[1] 
       << "}  mips: [" << m_mipsPmt[0]  << ',' << m_mipsPmt[1] 
       << "]  flags: <" << m_flags[0] << ',' << m_flags[1] << ">" << endl;
}
