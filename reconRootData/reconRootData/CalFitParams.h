
#ifndef reconRootData_CalFitParams_H
#define reconRootData_CalFitParams_H

#include <TObject.h>
#include <TVector3.h>
#include "CalParams.h"

/** 
* @class CalFitParams
*
* @brief Root clone for the TDS CalFitParams
*
* Note that not all the constructors implemented in the tds class are re-implemented
* here (rather, we put in only the ones that are actually used in the code,
* corresponding to the old class which was not inheriting from CalParams).
*
* @author Bill Atwood, Tracy Usher, David Chamont, Luca Baldini
*
* $Header$
*/

class CalFitParams : public CalParams
{
public:

  /// Default constructor
  CalFitParams() { Clear() ; }

  /// Direct construction from all the elements (the old fashioned way)
  CalFitParams
    ( Double_t energy, Double_t eneError,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz,
      Int_t numFitLayers, Double_t chiSquare ) ;

  void init
    ( Double_t energy, Double_t eneError,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz,
      Int_t numFitLayers, Double_t chiSquare ) ;
  
  /// Direct construction from all the elements (the old fashioned way).
  /// This reflect the way the CalFitParams class looked like before it inherited from
  /// the base class CalParams.
  CalFitParams
    ( Int_t numFitLayers, Double_t chiSquare,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz ) ;
  
  void init
    ( Int_t numFitLayers, Double_t chiSquare,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz ) ;
  
  virtual ~CalFitParams() {}
  
  /// Retrieve the energy
  inline const Int_t    getFitLayers() const {return m_nFitLayers;}
  inline const Double_t getChiSquare() const {return m_chiSquare;}
  
  void Clear( Option_t * option ="" ) ;
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const CalFitParams &, const std::string & name = "" ) const ; // for tests
  void Print( Option_t * option ="" ) const ;
  
 private:
  
  /// Energy and associated error
  Int_t      m_nFitLayers;
  Double32_t m_chiSquare;
  
  ClassDef(CalFitParams,3)
   
} ;

#endif

