#include "CalUtil/LinearConvertAdc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/MsgStream.h"


static ToolFactory<LinearConvertAdc> s_factory;
const IToolFactory& LinearConvertAdcFactory = s_factory;

LinearConvertAdc::LinearConvertAdc( const std::string& type, 
								   const std::string& name, 
								   const IInterface* parent)
								   : AlgTool(type,name,parent)
{
	// declare base interface for all consecutive concrete classes
	declareInterface<IConvertAdc>(this);
    // Declare the properties that may be set in the job options file
	
}


StatusCode LinearConvertAdc::initialize() {
    // Purpose and Method: initialize linear ConvertAdc model function .
    // Inputs: light attenuation parameter from detModel service.
	
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    StatusCode sc;
    
    IGlastDetSvc* detSvc;
	
    if( serviceLocator() ) {
		sc = serviceLocator()->service( "GlastDetSvc", detSvc, true );
    }
    if(sc.isFailure())
    {
		log << MSG::ERROR << "Could not find eventSvc" << endreq;
		return sc;
    }
	
	
    if(!detSvc->getNumericConstByName("cal.pedestal",&m_pedestal)) {
		log << MSG::ERROR << " constant " << " cal.pedestal not defined" << endreq;
		return StatusCode::FAILURE;
	} 
	
    
	if(!detSvc->getNumericConstByName("cal.maxResponse0",&m_maxResponse[0])) {
		log << MSG::ERROR << " constant " << " cal.maxResponse0 not defined" << endreq;
		return StatusCode::FAILURE;
	} 
	if(!detSvc->getNumericConstByName("cal.maxResponse1",&m_maxResponse[1])) {
		log << MSG::ERROR << " constant " << " cal.maxResponse1 not defined" << endreq;
		return StatusCode::FAILURE;
	} 
	if(!detSvc->getNumericConstByName("cal.maxResponse2",&m_maxResponse[2])) {
		log << MSG::ERROR << " constant " << " cal.maxResponse2 not defined" << endreq;
		return StatusCode::FAILURE;
	} 
	if(!detSvc->getNumericConstByName("cal.maxResponse3",&m_maxResponse[3])) {
		log << MSG::ERROR << " constant " << " cal.maxResponse3 not defined" << endreq;
		return StatusCode::FAILURE;
	} 
	
	if(!detSvc->getNumericConstByName("cal.maxAdcValue",&m_maxAdc)) {
		log << MSG::ERROR << " constant " << " cal.maxAdcValue not defined" << endreq;
		return StatusCode::FAILURE;
	} 
	
	for (int face=0; face < 2; face++) {
		for (int range = 0; range < 4; range++) {
			m_gain[face][range] = (m_maxAdc-m_pedestal)/m_maxResponse[range];
		}
	}
	
	
	return StatusCode::SUCCESS;
	
}


short unsigned int LinearConvertAdc::calculateAdc(idents::CalXtalId id,
												  idents::CalXtalId::XtalFace face,
												  idents::CalXtalId::AdcRange range,
												  double* depositedEnergy) {
	
	// range/2 is diode-type [0] is Large diode; [1] is small
	unsigned short int adc = (int)(depositedEnergy[range/2]*m_gain[face][range])
		+m_pedestal;
	if (adc > m_maxAdc) adc = m_maxAdc;

	return adc;
	
	
};

float LinearConvertAdc::calculateEnergy(idents::CalXtalId id,
										idents::CalXtalId::XtalFace face,
										idents::CalXtalId::AdcRange range,
										short unsigned int adc) {
	
	return (float)(adc-m_pedestal)/m_gain[face][range];
	
	
};

idents::CalXtalId::AdcRange LinearConvertAdc::findRange(idents::CalXtalId id,
														idents::CalXtalId::XtalFace face,
														double* depositedEnergy) {
	// loop over the 2 diodes and fetch the response. Find the first one below the 
	// max energy for its range.

	int br;
	for (int r = 0;r<4;r++){

		// Large = 0; Small = 1.
		int diode_type = r/2;
		double resp = depositedEnergy[diode_type];
		br=r;
		if( resp < m_maxResponse[r]) break;														
	}
	
	return (idents::CalXtalId::AdcRange)br;
}
