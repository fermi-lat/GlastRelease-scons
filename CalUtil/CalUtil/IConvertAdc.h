#ifndef _IConvertAdc_H
#define _IConvertAdc_H 1


// Include files
#include "GaudiKernel/IAlgTool.h"
#include "idents/CalXtalId.h"
#include <string>
#include <vector>

/** @class ConvertAdc
 * @brief Base class for CAL  functions to convert deposited energy to ADC
 *
 * Author:  R.Dubois
 *
*/

static const InterfaceID IID_ConvertAdc("ConvertAdc", 1 , 0);


class IConvertAdc : virtual public IAlgTool {
    
public:
    
    IConvertAdc() {}; 
    virtual ~IConvertAdc() {};

    static const InterfaceID& interfaceID() { return IID_ConvertAdc; }

    virtual short unsigned int calculateAdc(idents::CalXtalId id,
		idents::CalXtalId::XtalFace face,
		idents::CalXtalId::AdcRange range,
        double* depositedEnergy) = 0;
    
    virtual float calculateEnergy(idents::CalXtalId id,
        idents::CalXtalId::XtalFace face,
		idents::CalXtalId::AdcRange range,
        short unsigned int adc) = 0;

	virtual idents::CalXtalId::AdcRange findRange(idents::CalXtalId id,
		idents::CalXtalId::XtalFace face,
        double* depositedEnergy) = 0;
};


#endif // _IConvertAdc_H
