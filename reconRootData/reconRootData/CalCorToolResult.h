
#ifndef reconRootData_CalCorToolResult_H
#define reconRootData_CalCorToolResult_H

#include "CalParams.h"
#include <TObjArray.h>
#include <TString.h>

class CalCorToolResult : public TObject
{
public:

    CalCorToolResult() ;
    // it takes ownership of free params
    CalCorToolResult
      ( const TString & correctionName,
        UInt_t statusBits,
        Double_t chiSquare,
        const CalParams & params,
        const TObjArray & freeAttributes ) ;
    // it takes ownership of free params
    void init
      ( const TString & correctionName,
        UInt_t statusBits,
        Double_t chiSquare,
        const CalParams & params,
        const TObjArray & freeAttributes ) ;
    virtual ~CalCorToolResult() ;

    const TString & getCorrectionName() const
     { return m_correctionName ; }
    UInt_t  getStatusBits() const
     { return m_statusBits ; }
    const Double_t getChiSquare() const
     { return m_chiSquare ; }
    const CalParams & getParams() const
     { return m_params ; }
    const TObjArray & getFreeAttributes() const
     { return m_freeAttributes ; }
     
    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalCorToolResult &, const std::string & name = ""  ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    CalCorToolResult & operator=( const CalCorToolResult & ) ;

    TString m_correctionName ;
    UInt_t m_statusBits ;
    Double32_t m_chiSquare ;
    CalParams m_params ;
    TObjArray m_freeAttributes ;

    ClassDef(CalCorToolResult,2)

} ;
    
#endif        






