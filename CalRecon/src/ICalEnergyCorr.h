
#ifndef __ICalEnergyCorr_H
#define __ICalEnergyCorr_H 1

#include "GaudiKernel/IAlgTool.h"

/**   
* @class ICalEnergyCorr
*
* base class for energy leakage corrections 
*
*
* $Header$
*/

static const InterfaceID IID_ICalEnergyCorr("ICalEnergyCorr", 1 , 0) ;

class ICalEnergyCorr : virtual public IAlgTool {

public:

    // retrieve interface ID
    static const InterfaceID & interfaceID() { return IID_ICalEnergyCorr ; }
    
    //! constructor
	ICalEnergyCorr() {}
    //! destructor
    virtual ~ICalEnergyCorr() {}
    
    //! Calculate corrections on all cal clusters
    virtual StatusCode doEnergyCorr() =0 ;
    
} ;

#endif



