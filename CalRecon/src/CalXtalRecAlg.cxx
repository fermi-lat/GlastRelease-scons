#include "CalRecon/CalXtalRecAlg.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "geometry/Point.h"
#include "Event/Digi/CalDigi.h"
#include "Event/TopLevel/EventModel.h"
#include <map>


static const AlgFactory<CalXtalRecAlg>  Factory;
const IAlgFactory& CalXtalRecAlgFactory = Factory;

using namespace Event;

// constructor
CalXtalRecAlg::CalXtalRecAlg(const std::string& name, ISvcLocator* pSvcLocator):
Algorithm(name, pSvcLocator) { 

}


//################################################
StatusCode CalXtalRecAlg::initialize()
//################################################
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

        
    // extracting int constants
    double value;
    typedef std::map<int*,std::string> PARAMAP;
    PARAMAP param;
    param[&m_xNum]=        std::string("xNum");
    param[&m_yNum]=        std::string("yNum");
    param[&m_eTowerCAL]=   std::string("eTowerCAL");
    param[&m_eLATTowers]=  std::string("eLATTowers");
    param[&m_CALnLayer]=   std::string("CALnLayer");
    param[&m_nCsIPerLayer]=std::string("nCsIPerLayer");
    param[&m_nCsISeg]= std::string("nCsISeg");
    param[&m_eXtal]=       std::string("eXtal");
    param[&m_eDiodeMSmall]=std::string("eDiodeMSmall");
    param[&m_eDiodePSmall]=std::string("eDiodePSmall");
    param[&m_eDiodeMLarge]=std::string("eDiodeMLarge");
    param[&m_eDiodePLarge]=std::string("eDiodePLarge");
    param[&m_eMeasureX]=std::string("eMeasureX");
    param[&m_eMeasureY]=std::string("eMeasureY");
    param[m_noise]=std::string("cal.noiseLarge");
    param[m_noise+1]=std::string("cal.noiseSmall");
    param[m_ePerMeV+1]=std::string("cal.ePerMeVSmall");
    param[m_ePerMeV]=std::string("cal.ePerMevLarge");
    param[&m_pedestal]=std::string("cal.pedestal");
    param[&m_maxAdc]=std::string("cal.maxAdcValue");
    param[&m_thresh]=std::string("cal.zeroSuppressEnergy");
    
    // now try to find the GlastDevSvc service
    
//    IGlastDetSvc* detSvc;
    sc = service("GlastDetSvc", detSvc);
    
    
    for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
        if(!detSvc->getNumericConstByName((*it).second, &value)) {
            log << MSG::ERROR << " constant " <<(*it).second <<" not defined" << endreq;
            return StatusCode::FAILURE;
        } else *((*it).first)=value;
    }
    
    int nTowers = m_xNum * m_yNum;
    
    // extracting double constants
    
    typedef std::map<double*,std::string> DPARAMAP;
    DPARAMAP dparam;
    dparam[m_maxEnergy]=std::string("cal.maxResponse0");
    dparam[m_maxEnergy+1]=std::string("cal.maxResponse1");
    dparam[m_maxEnergy+2]=std::string("cal.maxResponse2");
    dparam[m_maxEnergy+3]=std::string("cal.maxResponse3");
    dparam[&m_lightAtt]=std::string("cal.lightAtt");
    dparam[&m_CsILength]=std::string("CsILength");
    
    for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++){
        if(!detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
            log << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
            return StatusCode::FAILURE;
        } 
    }
    
    for (int r=0; r<4;r++) m_maxEnergy[r] *= 1000.; // from GeV to MeV

	
	
	return sc;
}
//################################################
StatusCode CalXtalRecAlg::execute()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
	sc = retrieve();


	for (Event::CalDigiCol::const_iterator it = m_CalDigiCol->begin(); 
           it != m_CalDigiCol->end(); it++) {
               idents::CalXtalId xtalId = (*it)->getPackedId();
	   int lyr = xtalId.getLayer();
	   int towid = xtalId.getTower();
	   int icol  = xtalId.getColumn();
		
	   Event::CalXtalRecData* recData = new Event::CalXtalRecData((*it)->getMode(),xtalId);
	   
	   computeEnergy(recData, *it);
	   computePosition(recData);
	   m_CalXtalRecCol->push_back(recData);
	}

	return sc;
}
//################################################
StatusCode CalXtalRecAlg::finalize()
//################################################
{
	StatusCode sc = StatusCode::SUCCESS;

//	m_CalRecLogs->writeOut();

	return sc;
}

//################################################
StatusCode CalXtalRecAlg::retrieve()
//################################################
{
	
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

	m_CalXtalRecCol = 0;

 
	m_CalXtalRecCol = new CalXtalRecCol;

DataObject* pnode=0;

    sc = eventSvc()->retrieveObject( EventModel::CalRecon::Event /*"/Event/CalRecon"*/, pnode );
    
    if( sc.isFailure() ) {
        sc = eventSvc()->registerObject( EventModel::CalRecon::Event /*"/Event/CalRecon"*/,new DataObject);
        if( sc.isFailure() ) {
            
            log << MSG::ERROR << "Could not create CalRecon directory" << endreq;
            return sc;
        }
    }

    

	m_CalDigiCol = SmartDataPtr<Event::CalDigiCol>(eventSvc(),EventModel::Digi::CalDigiCol /*"/Event/Digi/CalDigis"*/); 


	 sc = eventSvc()->registerObject(EventModel::CalRecon::CalXtalRecCol /*"/Event/CalRecon/CalXtalRecCol"*/,m_CalXtalRecCol);
	return sc;
}

//----------------- private ----------------------
//################################################
void CalXtalRecAlg::computeEnergy(CalXtalRecData* recData, const Event::CalDigi* digi)
//################################################
{
	MsgStream log(msgSvc(), name());


	const Event::CalDigi::CalXtalReadoutCol& readoutCol = digi->getReadoutCol();
		
	for ( Event::CalDigi::CalXtalReadoutCol::const_iterator it = readoutCol.begin();
		      it !=readoutCol.end(); it++){
				int rangeP = it->getRange(idents::CalXtalId::POS); 
				int rangeM = it->getRange(idents::CalXtalId::NEG); 

				double adcP = it->getAdc(idents::CalXtalId::POS);	
				double adcM = it->getAdc(idents::CalXtalId::NEG);	

				double eneP = m_maxEnergy[rangeP]*(adcP-m_pedestal)/(m_maxAdc-m_pedestal);
				double eneM = m_maxEnergy[rangeM]*(adcM-m_pedestal)/(m_maxAdc-m_pedestal);
				
				CalXtalRecData::CalRangeRecData* rangeRec = new CalXtalRecData::CalRangeRecData(rangeP,eneP,rangeM,eneM);
				recData->addRangeRecData(*rangeRec);
		}		
	
}


//################################################
void CalXtalRecAlg::computePosition(CalXtalRecData* recData)
//################################################
{
	MsgStream log(msgSvc(), name());

	idents::CalXtalId xtalId = recData->getPackedId();	
	
	   int layer = xtalId.getLayer();
	   int tower = xtalId.getTower();
	   int col  = xtalId.getColumn();

			idents::VolumeIdentifier segm0Id;
			segm0Id.append(m_eLATTowers);
			segm0Id.append(tower/m_xNum);
			segm0Id.append(tower%m_xNum);
	  		segm0Id.append(m_eTowerCAL);
	  		segm0Id.append(layer);
	  		segm0Id.append(layer%2);
			segm0Id.append(col);
			segm0Id.append(m_eXtal);
			segm0Id.append(0);

            HepTransform3D transf;
			detSvc->getTransform3DByID(segm0Id,&transf);
			Vector vect0 = transf.getTranslation();

			
			idents::VolumeIdentifier segm11Id;
			for(int ifield = 0; ifield<fSegment; ifield++)segm11Id.append(segm0Id[ifield]);
			segm11Id.append(m_nCsISeg-1);

			detSvc->getTransform3DByID(segm11Id,&transf);
			Vector vect11 = transf.getTranslation();




	Point p0(0.,0.,0.);		
	Point pCenter = p0+(vect0+vect11)*0.5;
	Vector dirLog = 0.5*(vect11-vect0)*m_nCsISeg/(m_nCsISeg-1);	

	
	
	double eneNeg = recData->getEnergy(0,idents::CalXtalId::NEG);
	double enePos = recData->getEnergy(0,idents::CalXtalId::POS);
	double asym=0;
	if(enePos>0 && eneNeg>0)asym = (enePos-eneNeg)/(enePos+eneNeg);
	double slope = (1+m_lightAtt)/(1-m_lightAtt);

	Point pLog = pCenter+dirLog*asym*slope;
	
	(recData->getRangeRecData(0))->setPosition(pLog);

  
}
