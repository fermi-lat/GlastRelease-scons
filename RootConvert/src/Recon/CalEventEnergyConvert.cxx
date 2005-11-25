
#include "GaudiKernel/ObjectVector.h"
#include <RootConvert/Recon/CalEventEnergyConvert.h>
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include "TParameter.h"

namespace RootPersistence {

void convert( const Event::CalCorToolResult & tdsResult, CalCorToolResult & rootResult ) {

    CalParams rootParams ;
    convert(tdsResult.getParams(),rootParams) ;
      
    TObjArray rootFreeAttributes ;
    Event::CalCorEneValueMap::const_iterator freeAttribute ;
    for ( freeAttribute=tdsResult.begin() ;
          freeAttribute!=tdsResult.end() ;
          ++freeAttribute ) {
        const std::string & name = freeAttribute->first ;
        const double & value = freeAttribute->second ;
        rootFreeAttributes.Add(new TParameter<Double_t>(name.c_str(),value)) ;
    }
      
    rootResult.init
      ( tdsResult.getCorrectionName(),tdsResult.getStatusBits(),
        tdsResult.getChiSquare(),rootParams,rootFreeAttributes ) ;
 ;
}
 
void convert( const CalCorToolResult & rootResult, Event::CalCorToolResult & tdsResult ) {
    
    tdsResult.setCorrectionName(rootResult.getCorrectionName().Data()) ;

    UInt_t rootStatusBits = rootResult.getStatusBits() ;
    tdsResult.setStatusBits(rootStatusBits) ;

    tdsResult.setChiSquare(rootResult.getChiSquare()) ;

    Event::CalParams tdsParams ;
    convert(rootResult.getParams(),tdsParams) ;    
    tdsResult.setParams(tdsParams) ;
    
    tdsResult.clear() ;
    TIter rootFreeAttributesIter(&rootResult.getFreeAttributes()) ;
    TParameter<Double_t> * rootFreeAttribute ;
    while ((rootFreeAttribute=(TParameter<Double_t> *)rootFreeAttributesIter.Next())) {
        tdsResult.insert(Event::CalCorEneValuePair
          ( rootFreeAttribute->GetName(),
            rootFreeAttribute->GetVal() )) ;
    }

}

void convert( const Event::CalEventEnergy & tdsEventEnergy, CalEventEnergy & rootEventEnergy ) {

    TObjArray rootResults ;
    ObjectVector<Event::CalCorToolResult>::const_iterator tdsResult ;
    for ( tdsResult=tdsEventEnergy.begin() ;
          tdsResult!=tdsEventEnergy.end() ;
          ++tdsResult ) {
        CalCorToolResult * rootResult = new CalCorToolResult ;
        convert(**tdsResult,*rootResult) ;
        rootResults.Add(rootResult) ;
    }
      
    CalParams rootParams ;
    convert(tdsEventEnergy.getParams(),rootParams) ;
      
    rootEventEnergy.init
      (rootParams,tdsEventEnergy.getStatusBits(),rootResults) ;
}
 
void convert( const CalEventEnergy & rootEventEnergy, Event::CalEventEnergy & tdsEventEnergy ) {
    
    tdsEventEnergy.clear() ;
    
    Event::CalParams tdsParams ;
    convert(rootEventEnergy.getParams(),tdsParams) ;    
    tdsEventEnergy.setParams(tdsParams) ;
    
    UInt_t rootStatusBits = rootEventEnergy.getStatusBits() ;
    tdsEventEnergy.setStatusBits(rootStatusBits) ;

    TIter rootResultsIter(&rootEventEnergy.getToolsResults()) ;
    CalCorToolResult * rootResult ;
    while ((rootResult=(CalCorToolResult *)rootResultsIter.Next())) {
        Event::CalCorToolResult * tdsResult = new Event::CalCorToolResult ;
        convert(*rootResult,*tdsResult) ;
        tdsEventEnergy.push_back(tdsResult) ;
    }

}

}
