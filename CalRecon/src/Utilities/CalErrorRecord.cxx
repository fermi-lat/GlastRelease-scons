
// for class definition
#include "CalErrorRecord.h"
#include "GaudiKernel/MsgStream.h"
#include <iomanip>
#include <sstream>

void CalErrorRecord::print( MsgStream & log ) const {
    
    if (log.isActive()) {
        
        std::ostringstream lastTime ;
        lastTime
          <<std::setprecision(17)<<std::setw(25) 
          <<std::setiosflags(std::ios::scientific) 
          <<m_lastTime ;

        log<<endreq<<"**********************"<<endreq ;
        log<<"Run "<<m_run<<" Event "<<m_event 
            <<" Last Time "<<lastTime.str()<<endreq ;
        log<<m_catcherName<<": "<<m_comment<<endreq ;
        log<<"**********************"<<endreq ;

    }
    
}
