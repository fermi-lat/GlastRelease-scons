
//#include "Event/dataManager.h"
//#include "Event/messageManager.h"
#include "CalRecon/CalBase.h"
#include "CalRecon/CalPedCalib.h"
#include "CalRecon/CalCalibLogs.h"
#include "CalRecon/CalRecLogsAlg.h"
#include "CalRecon/CalADCLogs.h"
#include "CalRecon/CalRecLogs.h"
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

static const AlgFactory<CalRecLogsAlg>  Factory;
const IAlgFactory& CalRecLogsAlgFactory = Factory;


	// constructor
CalRecLogsAlg::CalRecLogsAlg(const std::string& name, ISvcLocator* pSvcLocator):
	  Algorithm(name, pSvcLocator) {
			declareProperty("PedFile",m_PedFileName);
			declareProperty("GainFile",m_GainFileName);
			declareProperty("IntlinFile",m_IntlinFileName);
			declareProperty("RailFile",m_RailFileName);
			declareProperty("SlopeFile",m_SlopeFileName);
			
	  }


//################################################
StatusCode CalRecLogsAlg::initialize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
/*	
	m_CalPedLogs   = dataManager::instance()->getData("CalPedCalib",m_CalPedLogs);
	m_CalCalibLogs = dataManager::instance()->getData("CalCalibLogs",m_CalCalibLogs);
	m_CalGeoLogs   = dataManager::instance()->geo()->cal();

	m_CalRawLogs   = dataManager::instance()->getData("CalRawLogs",m_CalRawLogs);
	m_CalRecLogs   = dataManager::instance()->getData("CalRecLogs",m_CalRecLogs);
*/
     sc = service("CalGeometrySvc", m_CalGeo);
  

	m_CalPedLogs  = new CalPedCalib(); 
	m_CalPedLogs->setFileName(m_PedFileName);
	m_CalPedLogs->make();

	m_CalCalibLogs = new CalCalibLogs();
	m_CalCalibLogs->setFileNames(m_IntlinFileName, m_GainFileName, m_RailFileName, m_SlopeFileName);
	m_CalCalibLogs->make();


//	m_CalRecLogs->clear();

	return sc;
}
//################################################
StatusCode CalRecLogsAlg::execute()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
	sc = retrieve();
	int nLogs = m_CalRawLogs->num();
	for (int jlog = 0; jlog < nLogs ; jlog++) {
		
		CalADCLog* ADCLog = m_CalRawLogs->Log(jlog);

		int ilayer          = ADCLog->layer();
		detGeo::axis view   = ADCLog->view();
	    int icol            = ADCLog->column();

//		detGeo*    geoLog = m_CalGeoLogs->getLog(ilayer,view,icol);
		detGeo geoLog = m_CalGeo->getLog(ilayer,view,icol);

		
		CalADCLog* pedLog = m_CalPedLogs->getLogID(CalLogID::ID(ilayer,view,icol));
		CalCalibLog* calibLog = m_CalCalibLogs->getLogID(CalLogID::ID(ilayer,view,icol));

		CalRecLog* recLog = m_CalRecLogs->getLogID(CalLogID::ID(ilayer,view,icol));
		computeEnergy(recLog, ADCLog, pedLog, calibLog);
		computePosition(recLog,&geoLog, calibLog);
//		std::cout << " ilayer = " << ilayer << " view=" << view << " icol=" << icol << std::endl;
	}
//	m_CalRecLogs->writeOut();

	return sc;
}
//################################################
StatusCode CalRecLogsAlg::finalize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;

//	m_CalRecLogs->writeOut();

	return sc;
}

//################################################
StatusCode CalRecLogsAlg::retrieve()
//################################################
{
	
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	m_CalRecLogs = 0;
	m_CalRecLogs = new CalRecLogs();

DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( "/Event/CalRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/CalRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }

    

	m_CalRawLogs = SmartDataPtr<CalADCLogs>(eventSvc(),"/Event/CalRecon/CalADCLogs"); 


	 sc = eventSvc()->registerObject("/Event/CalRecon/CalRecLogs",m_CalRecLogs);
	return sc;
}

//----------------- private ----------------------
//################################################
void CalRecLogsAlg::computeEnergy(CalRecLog* recLog, const CalADCLog* ADCLog,
								  const CalADCLog* pedLog, const CalCalibLog* calibLog)
//################################################
{
	
	double MAXPED = 4095.;
	double ADCSATURATION = 3850;
	bool eneSet = false;
	for (int irange = 0; irange < CALNRANGES ; irange++) {
		CalBase::RANGE r = (irange == 0? CalBase::LEX : CalBase::LE);
		double eneNeg = 0.;
    	double enePos = 0.;
		double adcNeg = 0.;
		double adcPos = 0.;
		double adcSatNeg = 0.;
		double adcSatPos = 0.;
		for (int iside = 0; iside < CALNSIDES; iside++) {
			CalBase::SIDE s = (iside == 0? CalBase::NEG : CalBase::POS);
			if (irange > 1) r = (irange == 2? CalBase::HEX : CalBase::HE);
			double adc = ADCLog->ADC(s,r);
			double ped = pedLog->ADC(s,r);
			double adc_ped = adc - ped;
			//if (s== CalBase::NEG) recLog->setNegADC(r,(adc-ped)/(MAXPED-ped));
			//else recLog->setPosADC(r,(adc-ped)/(MAXPED-ped));
			if (s== CalBase::NEG) recLog->setNegADC(r,adc_ped);
			else recLog->setPosADC(r,adc_ped);

			double adcSat = 0.6*(calibLog->getRail(iside,irange));
			double ene = calibLog->adc_to_MeV(adc_ped,iside,irange);
			iside == 0? eneNeg = ene: enePos = ene;
			iside == 0? adcNeg = adc_ped: adcPos = adc_ped;
			iside == 0? adcSatNeg = adcSat: adcSatPos = adcSat;
		}
		recLog->setNegEnergy(r,eneNeg);
		recLog->setPosEnergy(r,enePos);
		if(adcNeg<50 || adcPos<50)eneNeg=enePos=0;

		if (!eneSet && (adcNeg < adcSatNeg && adcPos < adcSatPos)) {
			eneSet = true;
			recLog->setNegEnergy(eneNeg);
			recLog->setPosEnergy(enePos);
            recLog->setBestRange(r);
		}
	}
}
//################################################
void CalRecLogsAlg::computePosition(CalRecLog* recLog, const detGeo* geoLog,
									const CalCalibLog* calibLog)
//################################################
{
	Point pCenter = geoLog->position();
	Point pSize   = geoLog->size();
	
	detGeo::axis view = recLog->view();
	double xdir = 0.;
	double ydir = 0.;
	ydir = (view == detGeo::X? 1. : 0);
	xdir = (view == detGeo::Y? 1. : 0);

	Vector dirLog(xdir,ydir,0.); 
	
	
	double eneNeg = recLog->negEnergy();
	double enePos = recLog->posEnergy();
	double asym = recLog->asymmetry();
	double slope = calibLog->getSlope(recLog->bestRange());

	Point pLog = pCenter+dirLog*asym*slope;

	recLog->setPosition(pLog);
}
