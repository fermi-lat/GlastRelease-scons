
//#include "Event/dataManager.h"
//#include "Event/messageManager.h"
#include "CalRecon/CalBase.h"
#include "CalRecon/CalPedCalib.h"
#include "CalRecon/CalCalibLogs.h"
#include "CalRecon/CalDigiLogsAlg.h"
#include "CalRecon/CalADCLogs.h"
#include "CalRecon/CalRecLogs.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

static const AlgFactory<CalDigiLogsAlg>  Factory;
const IAlgFactory& CalDigiLogsAlgFactory = Factory;
#include "GlastEvent/Digi/CalDigi.h"


// constructor
CalDigiLogsAlg::CalDigiLogsAlg(const std::string& name, ISvcLocator* pSvcLocator):
Algorithm(name, pSvcLocator) {
    declareProperty("PedFile",m_PedFileName);
    declareProperty("GainFile",m_GainFileName);
    declareProperty("IntlinFile",m_IntlinFileName);
    declareProperty("RailFile",m_RailFileName);
    declareProperty("SlopeFile",m_SlopeFileName);
    declareProperty("MuPeaksFile",m_MuPeaksFileName);
    declareProperty("ChargePeaksFile",m_ChargePeaksFileName);
    
}


//################################################
StatusCode CalDigiLogsAlg::initialize()
//################################################
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    /*	
    m_CalPedLogs   = dataManager::instance()->getData("CalPedCalib",m_CalPedLogs);
    m_CalCalibLogs = dataManager::instance()->getData("CalCalibLogs",m_CalCalibLogs);
    m_CalGeoLogs   = dataManager::instance()->geo()->cal();
    
      m_CalRawLogs   = dataManager::instance()->getData("CalRawLogs",m_CalRawLogs);
      m_CalRecLogs   = dataManager::instance()->getData("CalRecLogs",m_CalRecLogs);
    */
    sc = service("CalGeometrySvc", m_CalGeo);

    if(sc.isFailure())
    {
        log <<MSG::ERROR << "Error ICalGeometrySvc not found" << endreq;
        return sc;
    }
    
    
    
    m_CalPedLogs  = new CalPedCalib(); 
    m_CalPedLogs->setFileName(m_PedFileName);
    m_CalPedLogs->make();
    
    m_CalCalibLogs = new CalCalibLogs();
    m_CalCalibLogs->setFileNames(m_IntlinFileName, m_GainFileName, 
                        m_RailFileName, m_SlopeFileName, 
                        m_ChargePeaksFileName, m_MuPeaksFileName);
    m_CalCalibLogs->make();
    
    
    //	m_CalRecLogs->clear();
    
    return sc;
}
//################################################
StatusCode CalDigiLogsAlg::execute()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
	sc = retrieve();
//	int nLogs = m_CalRawLogs->num();
        for (CalDigiVector::const_iterator it = m_CalRawLogs->begin(); 
           it != m_CalRawLogs->end(); it++) {
	  //	for (int jlog = 0; jlog < nLogs ; jlog++) {
		
	  //		CalADCLog* ADCLog = m_CalRawLogs->Log(jlog);
               idents::CalLogId ADCLog = (*it)->getPackedId();
	   int lyr = 7-ADCLog.getLayer();
		int towid = ADCLog.getTower();
	   idents::ModuleId mod(towid);
	   int  ilayer=lyr/2; 
       int viewKludge = lyr%2;
       CalDetGeo::axis view   = CalDetGeo::makeAxis(viewKludge);
	   int icol            = ADCLog.getColumn();
		
	   log << MSG::DEBUG << " tower=" << towid << " layer=" << lyr
		   << " col=" << icol << " mod=" << mod << endreq;

	   //		detGeo*    geoLog = m_CalGeoLogs->getLog(ilayer,view,icol);
	   CalDetGeo geoLog = m_CalGeo->getLog(ilayer,view,icol,mod);

		
	   CalADCLog* pedLog = m_CalPedLogs->getLogID(CalLogID::ID(ilayer,view,icol,mod));
	   CalCalibLog* calibLog = m_CalCalibLogs->getLogID(CalLogID::ID(ilayer,view,icol,mod));

	   CalRecLog* recLog = m_CalRecLogs->getLogID(CalLogID::ID(ilayer,view,icol,mod));
	   computeEnergy(recLog, (*it), pedLog, calibLog);
	   computePosition(recLog,&geoLog, calibLog);
	   //		std::cout << " ilayer = " << ilayer << " view=" << view << " icol=" << icol << std::endl;
	}
	m_CalRecLogs->writeOut();

	return sc;
}
//################################################
StatusCode CalDigiLogsAlg::finalize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;

//	m_CalRecLogs->writeOut();

	return sc;
}

//################################################
StatusCode CalDigiLogsAlg::retrieve()
//################################################
{
	
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	m_CalRecLogs = 0;

	int nModX = m_CalGeo->numModulesX();   
	int nModY = m_CalGeo->numModulesY();   
	int nLogs = m_CalGeo->numLogs();
	int nLayers = m_CalGeo->numLayers();
	int nViews = m_CalGeo->numViews();
 
	m_CalRecLogs = new CalRecLogs(nModX,nModY,nLogs,nLayers,nViews);

DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( "/Event/CalRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/CalRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }

    

	m_CalRawLogs = SmartDataPtr<CalDigiVector>(eventSvc(),"/Event/Digi/CalDigis"); 


	 sc = eventSvc()->registerObject("/Event/CalRecon/CalRecLogs",m_CalRecLogs);
	return sc;
}

//----------------- private ----------------------
//################################################
void CalDigiLogsAlg::computeEnergy(CalRecLog* recLog, const CalDigi* ADCLog,
								  const CalADCLog* pedLog, const CalCalibLog* calibLog)
//################################################
{
	MsgStream log(msgSvc(), name());
	
		double ped=100.;
		double maxadc=4095;
		double maxEnergy[]={200.,1600.,12800.,102400.};
		double eneNeg = 0.;
    	double enePos = 0.;
		double adcNeg = 0.;
		double adcPos = 0.;

		for (int iside = 0; iside < CALNSIDES; iside++) {
			CalLogReadout::LogFace s = (iside == 0? CalLogReadout::NEG : CalLogReadout::POS);
			double adc = ADCLog->getAdc(0,s);	
			int irange = ADCLog->getRange(0,s); 
			CalBase::RANGE r = (irange == 0? CalBase::LEX : CalBase::LE);
			if (irange > 1) r = (irange == 2? CalBase::HEX : CalBase::HE);
//            double ped = pedLog->ADC(static_cast<CalBase::SIDE>(s),r);

			log << MSG::DEBUG << " adc=" << adc << " irange=" << irange
				<< " iside=" << iside << endreq;
			
			double adc_ped = adc - ped;
			if (s== CalLogReadout::NEG) recLog->setNegADC(r,adc_ped);
			else recLog->setPosADC(r,adc_ped);

//			double adcSat = 0.6*(calibLog->getRail(iside,irange));
//			double ene = calibLog->adc_to_MeV(adc_ped,iside,irange);

			double ene = maxEnergy[irange]*adc_ped/(maxadc-ped);

			if (iside == 0){
				eneNeg=ene;
				adcNeg = adc_ped;
				recLog->setNegEnergy(r,eneNeg);
		 		recLog->setNegEnergy(eneNeg);
			} else {
				enePos = ene;
				adcPos = adc_ped;
				recLog->setPosEnergy(r,enePos);
				recLog->setPosEnergy(enePos);
			}
	        recLog->setBestRange(r);

		}
		

		
	
}
//################################################
void CalDigiLogsAlg::computePosition(CalRecLog* recLog, const CalDetGeo* geoLog,
									const CalCalibLog* calibLog)
//################################################
{
	Point pCenter = geoLog->position();
	Point pSize   = geoLog->size();
	
	CalDetGeo::axis view = recLog->view();
	double xdir = 0.;
	double ydir = 0.;
	ydir = (view == CalDetGeo::X? pSize.y() : 0);
	xdir = (view == CalDetGeo::Y? pSize.x() : 0);

	Vector dirLog(xdir,ydir,0.); 
	
	
	double eneNeg = recLog->negEnergy();
	double enePos = recLog->posEnergy();
	double asym = recLog->asymmetry();
//	double slope = calibLog->getSlope(recLog->bestRange());
	double latt = 0.65;
	double slope = (1+latt)/(1-latt);

	Point pLog = pCenter-dirLog*asym*slope;

	recLog->setPosition(pLog);
}
