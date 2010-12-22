#ifndef reconRootData_CalXtalsParams_H
#define reconRootData_CalXtalsParams_H

#include <TObject.h>

/** 
* @class CalXtalsParams
*
* @brief Root clone for the TDS CalXtalsParams
* @author Luca Baldini, Johan Bregeon
*
*/

class CalXtalsParams
{
 public:
  
  /// Default constructor
  CalXtalsParams();
  
  /// Constructor from all members.
  CalXtalsParams(Int_t numXtals, Int_t numTruncXtals, Int_t numSaturatedXtals,
		 Double_t xtalRawEneSum, Double_t xtalCorrEneSum,
		 Double_t xtalEneRms, Double_t xtalEneSkewness);
  
  /// Destructor.
  virtual ~CalXtalsParams() {}
  
  /// Initialize the class members.
  void init(Int_t numXtals, Int_t numTruncXtals, Int_t numSaturatedXtals,
	    Double_t xtalRawEneSum, Double_t xtalCorrEneSum,
	    Double_t xtalEneRms, Double_t xtalEneSkewness);
  
  /// Retrieve class parameters...
  inline Int_t getNumXtals()             const { return m_numXtals; }
  inline Int_t getNumTruncXtals()        const { return m_numTruncXtals; }
  inline Int_t getNumSaturatedXtals()    const { return m_numSaturatedXtals; }
  inline Double_t getXtalRawEneSum()     const { return m_xtalRawEneSum; }
  inline Double_t getXtalCorrEneSum()    const { return m_xtalCorrEneSum; }
  inline Double_t getXtalEneRms()        const { return m_xtalEneRms; }
  inline Double_t getXtalEneSkewness()   const { return m_xtalEneSkewness; }
  
  void Clear( Option_t * option ="" ) ;
  void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
  Bool_t CompareInRange( const CalXtalsParams &, const std::string & name = "" ) const ; // for tests
  void Print( Option_t * option ="" ) const ;
  
 private:
  
  Int_t m_numXtals;
  Int_t m_numTruncXtals;
  Int_t m_numSaturatedXtals;
  Double_t m_xtalRawEneSum;
  Double_t m_xtalCorrEneSum;
  Double_t m_xtalEneRms;
  Double_t m_xtalEneSkewness;
  
  ClassDef(CalXtalsParams,1)
} ;




#endif
