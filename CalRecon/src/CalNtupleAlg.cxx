#define CalNtupleAlg_CPP 

// Include files

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"


#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/INTuple.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "ntupleWriterSvc/INtupleWriterSvc.h"
#include "CalRecon/CsIClusters.h"
#include "CalRecon/CalRecLogs.h"



//------------------------------------------------------------------------------
/*! \class CalNtupleAlg
\brief  alg to control writing of ntuple from calorimeter reconstruction

  */

class CalNtupleAlg : public Algorithm {
    
public:
    //! Constructor of this form must be provided
    CalNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator); 
    
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();
    
private:
    std::string m_tupleName;
	CsIClusterList* m_cls;
	CalRecLogs* m_crl;
    INTupleWriterSvc *m_ntupleWriteSvc;

};

//------------------------------------------------------------------------------
static const AlgFactory<CalNtupleAlg>  Factory;
const IAlgFactory& CalNtupleAlgFactory = Factory;
//------------------------------------------------------------------------------
/// 
CalNtupleAlg::CalNtupleAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator){

    declareProperty("tupleName",  m_tupleName="");

		m_cls = 0;
}

//------------------------------------------------------------------------------
/*! 
*/
StatusCode CalNtupleAlg::initialize() {
    
    StatusCode sc = StatusCode::SUCCESS;

    MsgStream log(msgSvc(), name());
        log << MSG::INFO << "Initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    setProperties();

    // get a pointer to our ntupleWriterSvc
    sc = service("ntupleWriterSvc", m_ntupleWriteSvc);

    if( sc.isFailure() ) {
        log << MSG::ERROR << "CalNtupleAlg failed to get the ntupleWriterSvc" << endreq;
        return sc;
    }

    return sc;
}



//------------------------------------------------------------------------------
StatusCode CalNtupleAlg::execute() {
     
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    m_cls  = SmartDataPtr<CsIClusterList>(eventSvc(),"/Event/CalRecon/CsIClusterList");

	if (m_cls == 0){ sc = StatusCode::FAILURE;
	
	        log << MSG::ERROR << "CalNtupleAlg failed to access CsIClusterList" << endreq;
			return sc;
	}


	double fit_ener,fitalpha,fitlambda,profchi2,eleak,start;
	float energy_sum;

	int nClust = m_cls->num();
	//cout << nClust << endl;
	for ( int icl = 0; icl<nClust;icl++){
		CsICluster* cl = m_cls->Cluster(icl);
		energy_sum = cl->energySum();

		float zpos = (cl->position()).z();
		float xpos = (cl->position()).x();
		float ypos = (cl->position()).y();

		const std::vector<double>& eneLayer = cl->getEneLayer();
		const std::vector<Vector>& posLayer = cl->getPosLayer();

		float trans_rms = cl->getRmsTrans();
		float long_rms = cl->getRmsLong();

		Vector caldir = cl->direction();
		float caltheta=acos(caldir.z());
		float calphi = 1.570796326794897;
		if(caldir.x()!=0) calphi = atan(caldir.y()/caldir.x());
		
		fit_ener = cl->getFitEnergy();
		profchi2 = cl->getProfChisq();
		fitalpha = cl->getCsiAlpha();
		fitlambda = cl->getCsiLambda();
		start = cl->getCsiStart();
		eleak = cl->energyLeak();
    
	
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_Energy_Deposit", energy_sum);
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_Energy_Inc_Prof", fit_ener);
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_Energy_Inc_LeakCorr", eleak);
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_Energy_Incident", fit_ener);  // for the moment
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Energy_Incident_CalTkr", energy_sum);  // to be updated

		const std::string name_eLayer = "Cal_eLayer";
		const char* digit[8]={"0","1","2","3","4","5","6","7"};
		for (int layer=0; layer<8; layer++)
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), (name_eLayer+digit[layer]).c_str(), eneLayer[layer]);
	
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_Z", zpos);
		
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_X", xpos);

		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_Y", ypos);

		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_transv_rms", trans_rms);
		
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_long_rms", long_rms);

		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_theta", caltheta);

		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_phi", calphi);

	}

	m_crl = SmartDataPtr<CalRecLogs>(eventSvc(),"/Event/CalRecon/CalRecLogs"); 
	if (m_crl == 0){ sc = StatusCode::FAILURE;
	
	        log << MSG::ERROR << "CalNtupleAlg failed to access CalRecLogs" << endreq;
			return sc;
	}
	
	int no_xtals=0;
	int no_xtals_trunc=0;
	int nLogs = m_crl->num();
	for (int jlog = 0; jlog < nLogs ; jlog++) {
		CalRecLog* recLog = m_crl->Log(jlog);

		double eneLog = recLog->energy();
		if(eneLog>0)no_xtals++;
		if(eneLog>0.01*energy_sum)no_xtals_trunc++;
	}

		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_No_Xtals", no_xtals);
		sc = m_ntupleWriteSvc->addItem(m_tupleName.c_str(), "Cal_No_Xtals_Trunc", no_xtals_trunc);

    return sc;
}


//------------------------------------------------------------------------------
StatusCode CalNtupleAlg::finalize() {
    StatusCode  sc = StatusCode::SUCCESS;
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "finalize CalNtupleAlg " << endreq;
 
    return StatusCode::SUCCESS;
}

//------------------------------------------------------------------------------



