
#ifndef __CalTransvOffsetTool_H
#define __CalTransvOffsetTool_H 1

#include "EnergyCorr.h"

/**   
* @class CalTransvOffsetTool
*
* @brief Correction of TransvOffset when tracker info available
*
*/


class CalTransvOffsetTool : public EnergyCorr
 {
  public :
    
    //! destructor
    CalTransvOffsetTool( const std::string& type, const std::string& name, const IInterface* parent);
    ~CalTransvOffsetTool() {} ; 
    
    StatusCode doEnergyCorr( Event::CalCluster * ) ;
 } ;

#endif



