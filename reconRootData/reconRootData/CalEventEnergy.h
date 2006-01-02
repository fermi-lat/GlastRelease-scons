#ifndef ROOT_CalEventEnergy_H
#define ROOT_CalEventEnergy_H

#include "CalCorToolResult.h"
#include "CalParams.h"
#include <TObject.h>
#include <TVector3.h>

class CalEventEnergy : public TObject
{ 
public:
    
    CalEventEnergy() ;
    // it takes ownership of tools results
    CalEventEnergy
     ( const CalParams & params,
       UInt_t statusBits,
       const TObjArray & toolsResults ) ;
    // it takes ownership of tools results
    void init
     ( const CalParams & params,
       UInt_t statusBits,
       const TObjArray & toolsResults ) ;
    virtual ~CalEventEnergy() ;

    const CalParams & getParams() const
     { return m_params ; }
    UInt_t  getStatusBits() const
     { return m_statusBits ; }
    const TObjArray & getToolsResults() const
     { return m_toolsResults ; }
     
    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalEventEnergy &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    CalParams m_params ;
    UInt_t m_statusBits ;
    TObjArray m_toolsResults ;

    ClassDef(CalEventEnergy,1)
};


#endif 


