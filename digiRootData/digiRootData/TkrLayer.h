
#ifndef TkrLayer_H
#define TkrLayer_H

#include "TObject.h"
#include "TObjArray.h"
#include "StripId.h"

/*! \class TkrLayer
\brief
 Stores all hit Silicon Strip detectors for a layer in the TKR,
 where a layer may cross tower boundaries.
 Layer is defined as in the UCSC BFEM User's Guide - where a layer is
 a single X or Y
 Apr 2000 Daniel Flath - Cleaned up to look like other files in package
 Jan 1999 Daniel Flath - ROOT HTML Documentation added
 Dec 1999 Daniel Flath - Creation
*/

class TkrLayer : public TObject
{
 private:
  /// Tracker plane number
  UChar_t m_plane;
  /// does this layer measure x or y?
  UChar_t m_xy; 
  /// ToT values for each of the 2 readout boards
  UInt_t m_ToT[2];
  /// ERRF values for each of the 2 readout boards
  UInt_t m_errf[2];	  
  TObjArray *m_strips; //->
  
  Int_t m_numStrips;

 public:
     // X will denote strips that measure X
     // Y denotes strips that measure Y
  typedef enum {
    X = 0,
    Y
  } TKRAxes;
  
  TkrLayer();
  TkrLayer(TObjArray *strips);
  virtual ~TkrLayer();
  void Clean(Option_t *option="");

  TObjArray* getStrips() { return m_strips; };
  void setToT(UInt_t right, UInt_t left);
  void setErrf(UInt_t right, UInt_t left);
  void setPlaneNum(UChar_t planeVal) { m_plane = planeVal; };
  void setXY(TKRAxes xyVal) { m_xy = (xyVal == X ? 0 : 1); };
  
  /// Retrieve the ToT for this layer for a specific controller
  Int_t getToT(UChar_t ctrlNum);
  Int_t getErrf(UChar_t ctrlNum);
  /// Plane is used as defined in the UCSC BFEM User's Guide - refers to an XY Pair.
  inline UShort_t getPlaneNum() const { return UShort_t(m_plane); };
  /// Layer is used as defined in the UCSC BFEM User's Guide - referring to a single X OR Y.
  UShort_t getLayerNum() const;
  /// Returns the orientation of the layer
  TKRAxes getXY() const { return (m_xy ? Y : X); };
  /// Returns the number of hits in the list of strips for this layer
  UInt_t getNumHits() { return m_strips->GetEntries(); };
  /// Retrieve a specific hit strip from the list
  StripId* getHit(int i) { return ( (i < m_strips->GetEntries()) ? (StripId*)(m_strips->At(i)) : 0); };
      
  /// Root >= 3.0 is now const correct for the Compare function
  Int_t Compare(const TObject *obj) const;
  Bool_t IsSortable() const;
  
  ClassDef(TkrLayer,3)	  // Information on a single tracker layer	   
};

#endif
