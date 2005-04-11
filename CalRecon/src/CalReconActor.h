
#ifndef __CalReconActor_H
#define __CalReconActor_H 1

#include "CalReconKernel.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/ISvcLocator.h"

/**   
* @class CalReconActor
* @brief Common behavior of CalRecon algorithms and tools
*
* Commonalities between all active entities from the package,
* namely algorithms, tools. Commonalities are :
*  -  access to CalReconKernel.
* 
* \author david Chamont (LLR)
* \todo if wize, enlarge its use to XtalRec
*/


class CalReconActor
 {
  public :
    
    CalReconActor() ;
    ~CalReconActor() ;
    
    // to be called within derived class at initialize step
    StatusCode initialize( ISvcLocator * ) ;

    // available after call to initialize()
    CalReconKernel * getKernel() ;
        
  protected :
  
    CalReconKernel * m_kernel ;
    
 } ;
 
#endif



