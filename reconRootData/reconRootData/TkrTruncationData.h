#ifndef ROOT_TKRTRUNCATIONDATA_H
#define ROOT_TKRTRUNCATIONDATA_H

#include "TObject.h"

#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif
typedef std::vector<UInt_t>   intVector;
typedef std::vector<Float_t> floatVector;

/** @class TkrTruncationData
 * This ROOT class is the flattened version of the TDS 
 *  TkrTruncationPlane class, with the sortId info 
 *  from the TkrTruncationInfo added
 * $Header$
 */
class TkrTruncationData : public TObject
{
 public:
  TkrTruncationData() {Clear();} 
  virtual ~TkrTruncationData() {Clear();}

  const Int_t getStatus() const { return m_status; }
  const intVector&   getStripCount() const   { return m_stripCount; }
  const intVector&   getStripNumber() const { return m_stripNumber; }
  const floatVector& getLocalX() const      { return m_localX; }
  const float        getPlaneZ() const      { return m_planeZ; }
  const int getTower() const { return m_tower;}
  const int getTray()  const { return m_tray;}
  const int getFace()  const { return m_face;}
  const int getView()  const { return m_view; }

  void setStatus(int status) { m_status = status; }
  void setStripCount(const intVector& stripCount)  { m_stripCount = stripCount; }
  void setStripNumber(const intVector& stripNumber) { m_stripNumber = stripNumber; }
  void setLocalX(const floatVector& localX)      { m_localX = localX; }
  void setPlaneZ(float planeZ)      { m_planeZ = planeZ; }
  void setTower(int tower)  { m_tower = tower;}
  void setTray(int tray)  { m_tray = tray;}
  void setFace(int face) { m_face = face;}
  void setView(int view) { m_view = view; }

  void Clear( Option_t * option ="" ) ;
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; 
  Bool_t CompareInRange( const TkrTruncationData & ref, const std::string & name = "" ) const ;
  void Print( Option_t * option ="" ) const ;

 private:
  ///From TkrTruncationPlane
  int m_status;
  intVector m_stripCount;
  intVector m_stripNumber;
  floatVector m_localX;
  float m_planeZ;
  ///From TkrTruncationInfo
  int m_tower;
  int m_tray;
  int m_face;
  int m_view;

  ClassDef(TkrTruncationData,1) 
};
#endif
