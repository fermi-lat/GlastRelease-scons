// $Header$

// Include files
#include "CalRecon/CalIRFAlg.h"
#include "CalRecon/CalRecLogs.h"
#include "CalRecon/CalDetGeo.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


#include "GlastEvent/data/TdGlastData.h"
#include "data/CsIData.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
// #define TUPLE 1

static const AlgFactory<CalIRFAlg>  Factory;
const IAlgFactory& CalIRFAlgFactory = Factory;


//------------------------------------------------------------------------------
/// Algorithm parameters which can be set at run time must be declared.
/// This should be done in the constructor.
CalIRFAlg::CalIRFAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator), m_detSvc(0) {
}


//------------------------------------------------------------------------------
/*! The "functional" part of the class: For the EmptyAlgorithm example they do
nothing apart from print out info messages.
NB in the initialize method: you must explicitly initialize the base class
before using any services (message service, event data service etc.) otherwise 
the behaviour will be unpredictable at best.
*/
StatusCode CalIRFAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;


    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    //Look for the geometry service
    StatusCode sc = service("CalGeometrySvc", m_CalGeo);
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the CalGeometrySvc!" << endreq;
    }
    
    // now try to find the GlastDevSvc service
     sc = service("GlastDetSvc", m_detSvc);
    
    
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
    }
    
    

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalIRFAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;

	int nModX = m_CalGeo->numModulesX();   
	int nModY = m_CalGeo->numModulesY();   
	int nLogs = m_CalGeo->numLogs();
	int nLayers = m_CalGeo->numLayers();
	int nViews = m_CalGeo->numViews();
    CalRecLogs* crl = new CalRecLogs(nModX,nModY,nLogs,nLayers,nViews);
	DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( "/Event/CalRecon", pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject("/Event/CalRecon",new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }
	 sc = eventSvc()->registerObject("/Event/CalRecon/CalRecLogs",crl);

    // get the CsiData object from the TDS by a converter
//    SmartDataPtr<TdGlastData> (eventSvc(), "/Event/Raw/TdCsIDatas");


    SmartDataPtr<TdGlastData> glastData(eventSvc(),"/Event/TdGlastData");
    
    if( 0==glastData) { 
       log << MSG::ERROR << "could not find \""<< "/Event/TdGlastData" <<"\"" << endreq;
       return StatusCode::FAILURE;
    }
    

    const CsIData* csi = glastData->getCsIData();


    
    // see what is there
//    csi->printOn(std::cout);

	double ene = 0.0;
	for (int l=0; l < m_CalGeo->numLayers();l++){
		for (int v=0; v< m_CalGeo->numViews(); v++){
			CalDetGeo::axis view = CalDetGeo::makeAxis(v);
			int ilayer = l*(m_CalGeo->numViews())+v;
			int nhits = csi->nHits(ilayer);
			if(nhits == 0)continue;
			for ( int ihit = 0; ihit < nhits; ihit++){
				idents::ModuleId mod = csi->moduleId(ilayer, ihit);
				idents::XtalId xtalid = csi->xtalId(ilayer,ihit);				
				int icol = xtalid.xtal()-1;

//  to fix  incorrect log numbering for X layers (with view==Y, numbers should increase with Y)  
				if(v == 1) icol = m_CalGeo->numLogs()-1 - icol;

				float enePos = csi->Rresp(ilayer,ihit)*1000.;
				float eneNeg = csi->Lresp(ilayer,ihit)*1000.;
				ene += (enePos+eneNeg)/2;
				CalRecLog* recLog = crl->getLogID(CalLogID::ID(l,view,icol,mod));
				recLog->setNegEnergy(eneNeg);
				recLog->setPosEnergy(enePos);
				CalDetGeo geoLog = m_CalGeo->getLog(l,view,icol,mod);

				Point pCenter = geoLog.position();
				Point pSize   = geoLog.size();

				double xdir = 0.;
				double ydir = 0.;
				ydir = (view == CalDetGeo::X? pSize.y() : 0);
				xdir = (view == CalDetGeo::Y? pSize.x() : 0);

				Vector dirLog(xdir,ydir,0.); 
				double asym = recLog->asymmetry();
//				std::cout << " dirLog = " << dirLog << std::endl;
				Point pLog = pCenter+dirLog*asym*3.0
					;

				recLog->setPosition(pLog);

			}
		}
	}

//    crl->writeOut();
//	std::cout << std::endl << " ene = " << ene << std::endl << std::endl;    

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalIRFAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    
    return StatusCode::SUCCESS;
}





