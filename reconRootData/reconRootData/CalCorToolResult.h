
#ifndef reconRootData_CalCorToolResult_H
#define reconRootData_CalCorToolResult_H

#include "CalParams.h"
#include <TObjArray.h>
#include <TString.h>

class CalCorToolResult : public TObjArray
{
public:

    CalCorToolResult() ;
    CalCorToolResult
      ( const TString & correctionName,
        UInt_t statusBits,
        Double_t chiSquare,
        const CalParams & params,
        const TObjArray & freeParams ) ;
    void init
      ( const TString & correctionName,
        UInt_t statusBits,
        Double_t chiSquare,
        const CalParams & params,
        const TObjArray & freeParams ) ;
    virtual ~CalCorToolResult() ;

    const TString & getCorrectionName() const
     { return m_correctionName ; }
    UInt_t  getStatusBits() const
     { return m_statusBits ; }
    const Double_t getChiSquare() const
     { return m_chiSquare ; }
    const CalParams & getParams() const
     { return m_params ; }
     
    void Clear( Option_t * option ="" ) ;
    void Fake( UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalCorToolResult & ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    CalCorToolResult & operator=( const CalCorToolResult & ) ;

    TString m_correctionName ;
    UInt_t m_statusBits ;
    Double_t m_chiSquare ;
    CalParams m_params ;

    ClassDef(CalCorToolResult,1)

} ;
    
#endif	






