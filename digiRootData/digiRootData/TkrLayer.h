
#ifndef TKRLAYER_H
#define TKRLAYER_H

#include "TObject.h"
#include "TClonesArray.h"
#include "StripID.h"

/*! \class TkrHit
\brief
 Apr 2000 Daniel Flath - Cleaned up to look like other files in package
 Jan 1999 Daniel Flath - ROOT HTML Documentation added
 Dec 1999 Daniel Flath - Creation
*/

class TkrLayer : public TObject
{
 private:
  UChar_t m_layer;	  // Tracker layer number
  UChar_t m_xy; 	  // Is this an x or y tile?
  UInt_t m_TOT[2];	  // TOT values for each of the 2 readout boards
  UInt_t m_errf[2];	  // ERRF values for each of the 2 readout boards
  TClonesArray *m_strips; //->
  
  Int_t m_numStrips;

 public:
  typedef enum {
    X = 0,
    Y
  } TKRAxes;
  
  TkrLayer();
  TkrLayer(TClonesArray *strips);
  virtual ~TkrLayer();
  void Clean(Option_t *option="");

  TClonesArray* getStrips() { return m_strips; };
  void setTOT(UInt_t right, UInt_t left);
  void setErrf(UInt_t right, UInt_t left);
  void setLayer(UChar_t layerVal) { m_layer = layerVal; };
  void setXY(TKRAxes xyVal) { m_xy = (xyVal == X ? 0 : 1); };
  
  Int_t getTOT(UChar_t ctrlNum);
  Int_t getErrf(UChar_t ctrlNum);
  inline UChar_t getLayer() const { return m_layer; };
  TKRAxes getXY() const { return (m_xy ? Y : X); };
      
  StripID* AddStrip();

  /// Root >= 3.0 is now const correct for the Compare function
  Int_t Compare(const TObject *obj) const;
  Bool_t IsSortable() const;
  
  ClassDef(TkrLayer,3)	  // Information on a single tracker layer	   
};

#endif
