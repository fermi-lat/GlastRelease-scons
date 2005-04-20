
#ifndef __CalTransvOffsetCorr_H
#define __CalTransvOffsetCorr_H 1

#include "CalEnergyCorr.h"

/**   
* @class CalTransvOffsetCorr
*
* @brief Correction of TransvOffset when tracker info available
*
*/


class CalTransvOffsetCorr : public CalEnergyCorr
 {
  public :
    
    //! destructor
    CalTransvOffsetCorr( const std::string& type, const std::string& name, const IInterface* parent);
    ~CalTransvOffsetCorr() {} ; 
    
    StatusCode doEnergyCorr( Event::CalCluster * ) ;
 } ;

#endif



