
#include <reconRootData/CalEventEnergyMap.h>
#include <commonRootData/RootDataUtil.h>
#include "Riostream.h"

ClassImp(CalEventEnergyMap)

void CalEventEnergyMap::Clear( Option_t * )
{
    this->Clear();
}

// dummy data, just for tests
void CalEventEnergyMap::Fake( Int_t ievent, UInt_t rank, Float_t randNum )
{
}

void CalEventEnergyMap::Print( Option_t * ) const
{
    this->Print() ;
}

Bool_t CalEventEnergyMap::CompareInRange( const CalEventEnergyMap & c, const std::string & name ) const 
{
    bool result = true ;
    return result ;
}

CalEventEnergyMap::CalEventEnergyMap()
{ Clear() ; }

CalEventEnergyMap::~CalEventEnergyMap()
{ 
    this->DeleteAll() ; 
}

