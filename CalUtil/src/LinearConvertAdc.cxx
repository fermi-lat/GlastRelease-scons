#include "CalUtil/LinearConvertAdc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDetDataSvc.h"
#include "CLHEP/Random/RandGauss.h"

#include "CalibData/CalibModel.h"
#include "CalibData/Cal/CalCalibPed.h"
#include "CalibData/Cal/CalCalibGain.h"

static ToolFactory<LinearConvertAdc> s_factory;
const IToolFactory& LinearConvertAdcFactory = s_factory;

//using namespace Event;

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
  
  //default m_maxResponse value 
  if(!detSvc->getNumericConstByName("cal.maxResponse0",&m_maxResponse[0])) {
    log<<MSG::ERROR<< " constant " << " cal.maxResponse0 not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  if(!detSvc->getNumericConstByName("cal.maxResponse1",&m_maxResponse[1])) {
    log<<MSG::ERROR<< " constant " << " cal.maxResponse1 not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  if(!detSvc->getNumericConstByName("cal.maxResponse2",&m_maxResponse[2])) {
    log<<MSG::ERROR<< " constant " << " cal.maxResponse2 not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  if(!detSvc->getNumericConstByName("cal.maxResponse3",&m_maxResponse[3])) {
    log<<MSG::ERROR<< " constant " << " cal.maxResponse3 not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  
  if(!detSvc->getNumericConstByName("cal.maxAdcValue",&m_maxAdc)) {
    log<<MSG::ERROR<< " constant " << " cal.maxAdcValue not defined" << endreq;
    return StatusCode::FAILURE;
  } 
     
  if(!detSvc->getNumericConstByName("cal.noiseLarge",m_noise)) {
    log<<MSG::ERROR<< " constant " << " cal.noiseLarge not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  
  if(!detSvc->getNumericConstByName("cal.noiseSmall",m_noise+1)) {
    log<<MSG::ERROR<< " constant " << " cal.noiseSmall not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  
  if(!detSvc->getNumericConstByName("cal.ePerMevLarge",m_ePerMeV)) {
    log<<MSG::ERROR<< " constant " << " cal.ePerMeVLarge not defined" << endreq;
    return StatusCode::FAILURE;
  } 
  
  if(!detSvc->getNumericConstByName("cal.ePerMeVSmall",m_ePerMeV+1)) {
    log<<MSG::ERROR<< " constant " << " cal.ePerMeVSmall not defined" << endreq;
    return StatusCode::FAILURE;
  } 

  if(!detSvc->getNumericConstByName("cal.zeroSuppressEnergy",&m_thresh)) {
    log<<MSG::ERROR<< " constant " << " cal.ePerMeVSmall not defined" << endreq;
    return StatusCode::FAILURE;
  } 
 
  for (int face=0; face < 2; face++) {
    for (int range = 0; range < 4; range++) {
      m_gain[face][range] = (m_maxAdc-m_pedestal)/m_maxResponse[range];
    }
  }
  
  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    log << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }
  
  return StatusCode::SUCCESS;
}

short unsigned int LinearConvertAdc::calculateAdc(idents::CalXtalId id,
                          idents::CalXtalId::XtalFace face,
                          idents::CalXtalId::AdcRange range,
                          double* depositedEnergy, 
                          CalibData::CalCalibPed *peds, 
                          CalibData::CalCalibGain *gains ) {
  
  //default adc value
  // range/2 is diode-type [0] is Large diode; [1] is small
  unsigned short int adc = (int)(depositedEnergy[range/2]*m_gain[face][range]
    +m_pedestal);

  //calibrated adc value
   if( gains && peds ){
    CalibData::RangeBase* pRangeGain = gains->getRange( id, 0, face );
    CalibData::Gain* pGain = dynamic_cast<CalibData::Gain * >(pRangeGain);

    CalibData::RangeBase* pRangePed = peds->getRange( id, range, face );
    CalibData::Ped* pPed = dynamic_cast<CalibData::Ped * >(pRangePed);
    
    adc = (int)(depositedEnergy[range/2]/pGain->getGain()+pPed->getAvr());
  }  

  if (adc > m_maxAdc) adc = (int) m_maxAdc;
  return adc;
};

idents::CalXtalId::AdcRange LinearConvertAdc::calculateAdcAndNoise(
                idents::CalXtalId id, idents::CalXtalId::XtalFace face,
                double* depositedEnergy, 
                CalibData::CalCalibPed *peds, CalibData::CalCalibGain *gains ){
  //calibrated adc value
  
  if( gains && peds ){
    CalibData::RangeBase* pRangePed=0;
    CalibData::Ped* pPed=0;
    CalibData::RangeBase* pRangeGain=0;
    CalibData::Gain* pGain=0;
    
    float noise[2]={ RandGauss::shoot(), RandGauss::shoot() };
    float rms[2];
    float cosAngle;

    // convert to adc in LARGE diode and add pedestal average
    for( short int range=1; range!=-1; --range ){
      pRangeGain = gains->getRange( id, range, face );
      pGain = dynamic_cast<CalibData::Gain * >(pRangeGain);
      depositedEnergy[range] /= pGain->getGain();

      pRangePed = peds->getRange( id, range, face );
      pPed = dynamic_cast<CalibData::Ped * >(pRangePed);

      rms[range]= pPed->getSig();

      depositedEnergy[range]+= pPed->getAvr();
    }
    
    // add pedestal in LARGE diode
    cosAngle= pPed->getCosAngle();
    if( cosAngle == 2 ){
      //complete correlation: backward compatibility
      depositedEnergy[0]+= rms[0]*noise[0];
      depositedEnergy[0]+= rms[1]*noise[0];
    } else {
      depositedEnergy[0]+= rms[0]*noise[0]*cosAngle;
      depositedEnergy[1]+= rms[1]*noise[1]*cosAngle;
      cosAngle= sqrt(1-cosAngle*cosAngle);
      depositedEnergy[0]-= rms[1]*noise[1]*cosAngle;
      depositedEnergy[0]+= rms[0]*noise[0]*cosAngle;
    }

    // apply zero_suppress
    if( (m_thresh>0) && 
        ( depositedEnergy[0]< (m_thresh/pGain->getGain()+pPed->getAvr()) )){
      depositedEnergy[0]=0;
      return idents::CalXtalId::LEX8;
    }

    // look for potential best range  for LARGE diode
    for( short int range=0; range<2; ++range ){
      depositedEnergy[range]= (int) depositedEnergy[range];
      if( depositedEnergy[range]<m_maxAdc )
        return (idents::CalXtalId::AdcRange) range;
      //if( depositedEnergy[range]>m_maxAdc ) depositedEnergy[range]= m_maxAdc;
    }    
    
    // convert to adc in SMALL diode and add pedestal average
    for( short int range=3; range!=1; --range ){
      CalibData::RangeBase* pRangeGain = gains->getRange( id, range, face );
      CalibData::Gain* pGain = dynamic_cast<CalibData::Gain * >(pRangeGain);
      depositedEnergy[range] /= pGain->getGain();
    
      pRangePed = peds->getRange( id, range, face );
      pPed = dynamic_cast<CalibData::Ped * >(pRangePed);
      rms[range]*= pPed->getSig();
      depositedEnergy[range]+= pPed->getAvr();
    }
    
    // add pedestal  in SMALL diode
    cosAngle= pPed->getCosAngle();
    if( cosAngle==2 ){
      // backward compatibility with xml files without pedestal correlations
      depositedEnergy[2]+= rms[0]*noise[1]*cosAngle;
      depositedEnergy[3]+= rms[1]*noise[1]*cosAngle;
    } else {
      depositedEnergy[2]+= rms[0]*noise[0]*cosAngle;
      depositedEnergy[3]+= rms[1]*noise[1]*cosAngle;

      cosAngle= sqrt(1-cosAngle*cosAngle);
      depositedEnergy[2]-= rms[1]*noise[1]*cosAngle;
      depositedEnergy[3]+= rms[0]*noise[0]*cosAngle;
    }

    //look for potential best range  
    depositedEnergy[2]= (int) depositedEnergy[2];
    if( depositedEnergy[2]<m_maxAdc ) return idents::CalXtalId::HEX8;
    //if( depositedEnergy[2]>m_maxAdc ) depositedEnergy[2]= m_maxAdc;
    depositedEnergy[3]= (int) depositedEnergy[3];
    if( depositedEnergy[3]>m_maxAdc ) depositedEnergy[3]= m_maxAdc;
    return idents::CalXtalId::HEX1;
  }  
  
  // default adc value
  // range/2 is diode-type [0] is Large diode; [1] is small
  for( short int diode=0; diode<2; ++diode ){
    depositedEnergy[2*diode+1]+=  RandGauss::shoot()*m_noise[diode]/m_ePerMeV[diode];
    //zero suppression 
    if( diode==0 && depositedEnergy[1]<m_thresh  ){
      depositedEnergy[0]=0.;
      return idents::CalXtalId::LEX8;
    }
      
    for( short int range=2*diode; range<2+2*diode; ++range ){
      depositedEnergy[range]=(int)depositedEnergy[2*diode+1]*m_gain[face][range];
      depositedEnergy[range]+=m_pedestal;
      if( depositedEnergy[range]<m_maxAdc ) 
        return (idents::CalXtalId::AdcRange) range;
    }
  }
  depositedEnergy[3]= m_maxAdc;
  return idents::CalXtalId::HEX1;
};

float LinearConvertAdc::calculateEnergy(idents::CalXtalId id,
                    idents::CalXtalId::XtalFace face,
                    idents::CalXtalId::AdcRange range,
                    short unsigned int adc,
                    CalibData::CalCalibPed *peds, 
                    CalibData::CalCalibGain *gains ){
  //calibrated adc value
   if(gains && peds){
    CalibData::RangeBase* pRangeGain = gains->getRange( id, range, face );
    CalibData::Gain* pGain = dynamic_cast<CalibData::Gain * >(pRangeGain);

    CalibData::RangeBase* pRangePed = peds->getRange( id, range, face );
    CalibData::Ped* pPed = dynamic_cast<CalibData::Ped * >(pRangePed);
    
    return (float)(adc-pPed->getAvr())*pGain->getGain();
  }  

  //default value
  return (float)(adc-m_pedestal)/m_gain[face][range];
};

idents::CalXtalId::AdcRange LinearConvertAdc::findRange(idents::CalXtalId id,
                            idents::CalXtalId::XtalFace face,
                            double* depositedEnergy,
                            CalibData::CalCalibPed *peds, 
                            CalibData::CalCalibGain *gains ) {
  // loop over the 2 diodes and fetch the response. Find the first one below the
  // max energy for its range.
  
  int br;
  
  if( (gains==0) || (peds==0) ){
    //default value
    for (int range = 0;range<4;range++){
      // Large = 0; Small = 1.
      int diode_type = range/2;
      double resp = depositedEnergy[diode_type];
      br=range;
    
      if( resp < m_maxResponse[range]) break;                            
    }
  } else {
    //calibration value
    for ( int range = 0; range < 4; range++) {
      // Large = 0; Small = 1.
      int diode_type = range/2;
      double resp = depositedEnergy[diode_type];
      br=range;

      CalibData::RangeBase* pRangeGain = gains->getRange( id, range, face );
      CalibData::Gain* pGain = dynamic_cast<CalibData::Gain * >(pRangeGain);
  
      CalibData::RangeBase* pRangePed = peds->getRange( id, range, face );
      CalibData::Ped* pPed = dynamic_cast<CalibData::Ped * >(pRangePed);
    
      if( resp < (m_maxAdc-pPed->getAvr())/pGain->getGain() ) break;
    }
  }
  return (idents::CalXtalId::AdcRange)br;
}
