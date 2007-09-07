#include "reconRootData/TkrTruncationData.h"
#include "Riostream.h"

ClassImp(TkrTruncationData)

void TkrTruncationData::Clear(Option_t* /* option */) 
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
