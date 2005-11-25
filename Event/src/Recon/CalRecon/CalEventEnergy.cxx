
#include "Event/Recon/CalRecon/CalEventEnergy.h"

void Event::CalEventEnergy::clear() {
    
    m_statusBits = 0 ;
    m_params.clear() ;
    CalCorToolResultCol::iterator corIter ;
    for ( corIter=begin() ; corIter!=end() ; ++corIter ) {
        delete (*corIter) ;
    }
    CalCorToolResultCol::clear() ;

}

const Event::CalCorToolResult * Event::CalEventEnergy::findLast( const std::string & correctionName ) const {

    CalCorToolResultCol::const_reverse_iterator corIter ;
    for (
    corIter = rbegin() ; 
    corIter != rend() ; 
    corIter++ ) {
        const Event::CalCorToolResult * corResult = *corIter ;
        if ( corResult->checkStatusBit(CalCorToolResult::VALIDPARAMS) &&
             ( corResult->getCorrectionName() == correctionName ) ) {
            return corResult ;
        }
    }
    return 0 ;
    
}


