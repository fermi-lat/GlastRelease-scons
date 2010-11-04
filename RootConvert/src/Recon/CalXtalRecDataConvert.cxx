
#include <RootConvert/Recon/CalXtalRecDataConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalXtalRecData & tdsXtalRecData, CalXtalRecData & rootXtalRecData ) {
   
    rootXtalRecData.Clear() ;
    UInt_t    statusBits = tdsXtalRecData.getStatusWord();
    CalXtalId rootId     = RootPersistence::convert(tdsXtalRecData.getPackedId()) ; 
    idents::CalXtalId::CalTrigMode tdsMode = tdsXtalRecData.getMode() ;   
    if ( tdsMode == idents::CalXtalId::BESTRANGE ) {
        
        rootXtalRecData.initialize(statusBits,CalXtalId::BESTRANGE,rootId) ;   
        const Event::CalXtalRecData::CalRangeRecData * tdsXtalRange  
         = tdsXtalRecData.getRangeRecData(0) ;   
        CalRangeRecData rootRangeRecData(   
            tdsXtalRange->getRange(idents::CalXtalId::POS),    
            tdsXtalRange->getEnergy(idents::CalXtalId::POS),    
            tdsXtalRange->getRange(idents::CalXtalId::NEG),   
            tdsXtalRange->getEnergy(idents::CalXtalId::NEG)
        ) ;   
        rootRangeRecData.initialize(RootPersistence::convertPoint(tdsXtalRange->getPosition())) ;   
        rootXtalRecData.addRangeRecData(rootRangeRecData) ;
        
    } else {
           
        rootXtalRecData.initialize(statusBits,CalXtalId::ALLRANGE,rootId) ;         
        int range ;   
        for ( range = idents::CalXtalId::LEX8 ;
              range <= idents::CalXtalId::HEX1 ;
              ++range ) {
                  
            const Event::CalXtalRecData::CalRangeRecData * tdsXtalRange  
             = tdsXtalRecData.getRangeRecData(range) ;   
            if (tdsXtalRange) {
                CalRangeRecData rootRangeRecData(   
                    tdsXtalRange->getRange(idents::CalXtalId::POS),    
                    tdsXtalRange->getEnergy(idents::CalXtalId::POS),    
                    tdsXtalRange->getRange(idents::CalXtalId::NEG),   
                    tdsXtalRange->getEnergy(idents::CalXtalId::NEG) );   
                rootRangeRecData.initialize(RootPersistence::convertPoint(tdsXtalRange->getPosition())) ;   
                rootXtalRecData.addRangeRecData(rootRangeRecData) ;
            }    
        }   
        
    }   
}
 
void convert( const CalXtalRecData & rootXtalRecData, Event::CalXtalRecData & tdsXtalRecData ) {

        unsigned int      statusBits = rootXtalRecData.getStatusWord();
        idents::CalXtalId tdsId      = RootPersistence::convert(rootXtalRecData.getPackedId()) ;        
        if (rootXtalRecData.getMode()==CalXtalId::ALLRANGE) {
            
            tdsXtalRecData.initialize(statusBits,idents::CalXtalId::ALLRANGE,tdsId) ;
            unsigned int range ;
            for ( range = idents::CalXtalId::LEX8 ;
                  range < idents::CalXtalId::HEX1 ;
                  ++range) {    
                const CalRangeRecData * rootXtalRange
                  = rootXtalRecData.getRangeRecData(range) ;
                if (rootXtalRange) {
                    Event::CalXtalRecData::CalRangeRecData tdsXtalRange(
                            rootXtalRange->getRange(CalXtalId::POS),
                            rootXtalRange->getEnergy(CalXtalId::POS),
                            rootXtalRange->getRange(CalXtalId::NEG),
                            rootXtalRange->getEnergy(CalXtalId::NEG)
                        );
                    tdsXtalRange.setPosition(RootPersistence::convertPoint(rootXtalRange->getPosition()));
                    tdsXtalRecData.addRangeRecData(tdsXtalRange);
               }
            }
            
        } else if (rootXtalRecData.getMode() == CalXtalId::BESTRANGE) {
            
            tdsXtalRecData.initialize(statusBits,idents::CalXtalId::BESTRANGE,tdsId) ;
            const CalRangeRecData * rootXtalRange
              = rootXtalRecData.getRangeRecData(0) ;   
            
            Event::CalXtalRecData::CalRangeRecData tdsXtalRange(
                    rootXtalRange->getRange(CalXtalId::POS),
                    rootXtalRange->getEnergy(CalXtalId::POS),
                    rootXtalRange->getRange(CalXtalId::NEG),
                    rootXtalRange->getEnergy(CalXtalId::NEG)
                ) ;
            tdsXtalRange.setPosition(RootPersistence::convertPoint(rootXtalRange->getPosition()));
            tdsXtalRecData.addRangeRecData(tdsXtalRange) ;
            
        }
    }

}
