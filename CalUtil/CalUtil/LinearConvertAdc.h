#ifndef _LinearConvertAdc_H
#define _LinearConvertAdc_H 1


// Include files

#include "IConvertAdc.h"
#include "GaudiKernel/AlgTool.h"

/** @class LinearConvertAdc
 * @brief Linear ConvertAdc model implementation for crystal digitization
 *
 * Author:  R.Dubois
 *
*/

class LinearConvertAdc : public AlgTool, virtual public IConvertAdc {
    
public:
    
    LinearConvertAdc( const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~LinearConvertAdc() {};

    virtual short unsigned int calculateAdc(idents::CalXtalId id,
        idents::CalXtalId::XtalFace face,
		idents::CalXtalId::AdcRange range,
        double* depositedEnergy);
    
    virtual float calculateEnergy(idents::CalXtalId id,
        idents::CalXtalId::XtalFace face,
		idents::CalXtalId::AdcRange range,
        short unsigned int adc);

    ///Implementation of the method provided by the base class AlgTool.
    virtual StatusCode initialize();

	virtual idents::CalXtalId::AdcRange findRange(idents::CalXtalId id,
		idents::CalXtalId::XtalFace face,
        double* depositedEnergy);


private:

    /// overall gain factor: MeV/channel
    double m_gain[2][4];
    /// pedestal: channel
    double m_pedestal;
    /// max energy responses per gain range
    double m_maxResponse[4];
    /// max number of ADC channels
	double m_maxAdc;
};


#endif // _GlastDigi_LinearConvertAdc_H
