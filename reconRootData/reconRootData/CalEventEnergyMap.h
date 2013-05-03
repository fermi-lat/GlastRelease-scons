#ifndef ROOT_CalEventEnergyMap_H
#define ROOT_CalEventEnergyMap_H

#include <TMap.h>

class CalEventEnergyMap : public TMap
{ 
public:
    
    CalEventEnergyMap() ;
    virtual ~CalEventEnergyMap() ;
     
    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalEventEnergyMap &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    ClassDef(CalEventEnergyMap,1)
};


#endif 


