// $Header$

// Include files
#include "CalRecon/CalIRFAlg.h"
#include "CalRecon/CalRecLogs.h"
#include "TkrRecon/detGeo.h"
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/AlgFactory.h"
#include "Gaudi/Interfaces/IDataProviderSvc.h"
#include "Gaudi/DataSvc/SmartDataPtr.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "xml/IFile.h"

#include "reconstruction/GlastTuple.h"
#include "reconstruction/PrintReconData.h"
#include "reconstruction/SummaryData.h"
#include "reconstruction/GlastTuple.h"

#include "GlastEvent/data/TdGlastData.h"
#include "data/CsIData.h"
#include "gui/DisplayControl.h"
#include "GuiSvc/IGuiSvc.h"
#include "gui/GuiMgr.h"
#include "CalRecon/calorimeterGeo.h"
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
    
    // now try to find the GlastDevSvc service
    StatusCode sc = service("GlastDetSvc", m_detSvc);
    
    
    if (!sc.isSuccess ()){
        log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
    }
    
    // test: get a constant from the ini file
    m_ini = const_cast<xml::IFile*>(m_detSvc->iniFile()); //OOPS!
    int nx = m_ini->getInt("glast", "xNum");
    
    
    // get the Gui service
    IGuiSvc* guiSvc=0;
    sc = service("GuiSvc", guiSvc);

    if (!sc.isSuccess ()){
        log << MSG::WARNING << "Couldn't find the GuiSvc!" << endreq;
        sc =StatusCode::SUCCESS; 
    }else
    {
       guiSvc->guiMgr()->display().add(m_recon->displayRep(), "Cal reco");
    }

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalIRFAlg::execute() {
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::INFO << "execute" << endreq;

    
    CalRecLogs* crl = new CalRecLogs;
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
	for (int l=0; l < calorimeterGeo::numLayers();l++){
		for (int v=0; v< calorimeterGeo::numViews(); v++){
			detGeo::axis view = detGeo::makeAxis(v);
			int ilayer = l*(calorimeterGeo::numViews())+v;
			int nhits = csi->nHits(ilayer);
			if(nhits == 0)continue;
			for ( int ihit = 0; ihit < nhits; ihit++){
				idents::ModuleId mod = csi->moduleId(ilayer, ihit);
				idents::XtalId xtalid = csi->xtalId(ilayer,ihit);				
				int icol = xtalid.xtal()-1;

//  to fix  incorrect log numbering for X layers (with view==Y, numbers should increase with Y)  
				if(v == 1) icol = calorimeterGeo::numLogs()-1 - icol;

				float enePos = csi->Rresp(ilayer,ihit);
				float eneNeg = csi->Lresp(ilayer,ihit);
				ene += (enePos+eneNeg)/2;
				CalRecLog* recLog = crl->getLogID(CalLogID::ID(l,view,icol,mod));
				recLog->setNegEnergy(eneNeg);
				recLog->setPosEnergy(enePos);
				detGeo geoLog = calorimeterGeo::getLog(l,view,icol,mod);

				Point pCenter = geoLog.position();
				Point pSize   = geoLog.size();

				double xdir = 0.;
				double ydir = 0.;
				ydir = (view == detGeo::X? pSize.y() : 0);
				xdir = (view == detGeo::Y? pSize.x() : 0);

				Vector dirLog(xdir,ydir,0.); 
				double asym = recLog->asymmetry();
//				std::cout << " dirLog = " << dirLog << std::endl;
				Point pLog = pCenter+dirLog*asym*3.0
					;

				recLog->setPosition(pLog);

			}
		}
	}

    crl->writeOut();
	std::cout << std::endl << " ene = " << ene << std::endl << std::endl;    

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalIRFAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize" << endreq;
    delete m_recon;
    //    delete m_summary;
    
    return StatusCode::SUCCESS;
}

/*! Can be used to check to see if the newstyle NTuple was written properly.
*/
StatusCode CalIRFAlg::printNewNTuple() {
    StatusCode status;
    MsgStream log(msgSvc(), name());
    
    NTuplePtr nt = m_gsummary->tuple()->getNTuple();
    if(nt)
    {
        NTuple::Item<float> data;
        status = nt->item("CsI_eLayer1" ,data);
        
        if(status.isSuccess())
        {
            log << MSG::INFO << "Test Value of New Ntuple :\n"  
                << "CsI_eLayer1: " << data << "\n" << endreq;
            return StatusCode::SUCCESS;
        } else {
            log << MSG::ERROR << "Didn't load the Item!" << endreq;
            return StatusCode::FAILURE;
        }
    } else {
        log << MSG::ERROR << "Dead NTuple !" << endreq;
        return StatusCode::FAILURE;
    }
}





