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
// rather not pollute the namespace with these typedefs
//typedef std::vector<UInt_t>   intVector;
//typedef std::vector<Float_t> floatVector;


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
  const vector<UInt_t>&   getStripCount() const   { return m_stripCount; }
  const vector<UInt_t>&   getStripNumber() const { return m_stripNumber; }
  const vector<Float_t>& getLocalX() const      { return m_localX; }
  const Float_t        getPlaneZ() const      { return m_planeZ; }
  const Int_t getTower() const { return m_tower;}
  const Int_t getTray()  const { return m_tray;}
  const Int_t getFace()  const { return m_face;}
  const Int_t getView()  const { return m_view; }

  void setStatus(int status) { m_status = status; }
  void setStripCount(const vector<UInt_t>& stripCount)  { m_stripCount = stripCount; }
  void setStripNumber(const vector<UInt_t>& stripNumber) { m_stripNumber = stripNumber; }
  void setLocalX(const vector<Float_t>& localX)      { m_localX = localX; }
  void setPlaneZ(float planeZ)      { m_planeZ = planeZ; }
  void setTower(int tower)  { m_tower = tower;}
  void setTray(int tray)  { m_tray = tray;}
  void setFace(int face) { m_face = face;}
  void setView(int view) { m_view = view; }

  void Clear( Option_t * option ="" ) ;
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; 
  Bool_t CompareInRange( const TkrTruncationData & ref, const std::string & name = "" ) const ;
  void Print( Option_t * option ="" ) const ;

  /// Override the "new" operator in order to use the ReconObjectManager pool management
  void* operator new(size_t size);

  void* operator new(size_t size, void* vp);

  /// Override the "delete" operator in order to use the ReconObjectManager pool management
  void  operator delete(void* p);


 private:
  ///From TkrTruncationPlane
  Int_t m_status;
  ///From TkrTruncationInfo
  Int_t m_tower;
  Int_t m_tray;
  Int_t m_face;
  Int_t m_view;
  Float_t m_planeZ;
  vector<UInt_t> m_stripCount;
  vector<UInt_t> m_stripNumber;
  vector<Float_t> m_localX;

  ClassDef(TkrTruncationData,1) 
};

#endif
