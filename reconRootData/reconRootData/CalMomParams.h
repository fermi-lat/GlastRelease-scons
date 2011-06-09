
#ifndef reconRootData_CalMomParams_H
#define reconRootData_CalMomParams_H

#include <TObject.h>
#include <TVector3.h>
#include "CalParams.h"

/** 
* @class CalMomParams
*
* @brief Root clone for the TDS CalMomParams
* @author Luca Baldini, Johan Bregeon
*
*/

class CalMomParams : public CalParams
{
 public:
  
  /// Default constructor
  CalMomParams() { Clear() ; }
  
  /// Direct construction from all the elements (the old fashioned way)
  CalMomParams
    ( Double_t energy, Double_t eneError,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz,
      Int_t numIterations, Int_t numCoreXtals, Int_t numXtals,
      Double_t transRms, Double_t longRms, Double_t longRmsAsym, Double_t longSkewness,
      Double_t coreEnergyFrac, Double_t fullLength, Double_t dEdxSpread,
      Double_t minGhostDoca ) ;
  
  /// init
  void init
    ( Double_t energy, Double_t eneError,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz,
      Int_t numIterations, Int_t numCoreXtals, Int_t numXtals,
      Double_t transRms, Double_t longRms, Double_t longRmsAsym, Double_t longSkewness,
      Double_t coreEnergyFrac, Double_t fullLength, Double_t dEdxSpread,
      Double_t minGhostDoca) ;
  
  virtual ~CalMomParams() {}
  
  /// Retrieve the class members.
  inline const Int_t getNumIterations()     const { return m_numIterations; }
  inline const Int_t getNumCoreXtals()      const { return m_numCoreXtals; }
  inline const Int_t getNumXtals()          const { return m_numXtals; }
  inline const Double_t getTransRms()       const { return m_transRms; }
  inline const Double_t getLongRms()        const { return m_longRms; }
  inline const Double_t getLongRmsAsym()    const { return m_longRmsAsym; }
  inline const Double_t getLongSkewness()   const { return m_longSkewness; }
  inline const Double_t getCoreEnergyFrac() const { return m_coreEnergyFrac; }
  inline const Double_t getFullLength()     const { return m_fullLength; }
  inline const Double_t getdEdxSpread()     const { return m_dEdxSpread; }
  inline const Double_t getMinGhostDoca()   const { return m_minGhostDoca; }
  
  Double_t getElongation() const;
  Double_t getdEdxAverage() const;
  
  void Clear( Option_t * option ="" ) ;
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const CalMomParams &, const std::string & name = "" ) const ; // for tests
  void Print( Option_t * option ="" ) const ;
  
 private:
  
  Int_t m_numIterations;
  Int_t m_numCoreXtals;
  Int_t m_numXtals;
  Double_t m_transRms;
  Double_t m_longRms;
  Double_t m_longRmsAsym;
  Double_t m_longSkewness;
  Double_t m_coreEnergyFrac;
  Double_t m_fullLength;
  Double_t m_dEdxSpread;
  Double_t m_minGhostDoca;
  
  ClassDef(CalMomParams,3)

} ;


#endif

