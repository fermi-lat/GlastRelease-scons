#include "IGcrReconTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/SmartRefVector.h"

#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"


#include "TkrUtil/ITkrGeometrySvc.h"


#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h" 
#include "CLHEP/Vector/LorentzVector.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/GcrReconClasses.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/TopLevel/Event.h"

#include "Event/MonteCarlo/McIntegratingHit.h"
#include "Event/MonteCarlo/McPositionHit.h"
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/MonteCarlo/McParticle.h"

//needed to get CALEnergyRaw:
#include "Event/Recon/CalRecon/CalCluster.h"

//#include "TkrUtil/ITkrGeometrySvc.h"

#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"
#include "idents/CalXtalId.h"

#include "TkrRecon/../src/Filter/ITkrFilterTool.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"

#include "geometry/Vector.h"
#include "geometry/Ray.h"

#include "TMath.h"
#include "TObjArray.h"
#include "TObject.h"

#include "OnboardFilterTds/FilterStatus.h"
#include "OnboardFilterTds/ObfFilterStatus.h"
#include "OnboardFilterTds/FilterAlgTds.h"

#include "enums/TriggerBits.h"


/**   
 * @class GcrReconTool
 *
 */

//-----------------------------------------------------------------------------------------------------------------
class GcrReconTool : public IGcrReconTool,  public AlgTool 
{
public:    
  GcrReconTool(const std::string & type, const std::string & name, const IInterface * parent );
  virtual ~GcrReconTool() {};
  
  /// @brief Intialization of the tool
  virtual StatusCode initialize();

  virtual bool GcrReconTool::TriggerEngine4ON();
  //virtual bool GcrReconTool::OBF_HFCVetoExist();
  virtual bool GcrReconTool::checkFilters();

  virtual StatusCode GcrReconTool::findGcrXtals(std::string initDir);
  
  
private:

  /// PRIVATE METHODS
    
  DataSvc* getDataSvc(){return m_dataSvc;} 
  void setDataSvc(DataSvc* dataSvc){m_dataSvc = dataSvc;}
  Event::GcrXtalVec getGcrXtalVec(){return m_gcrXtalVec;} 
  void setGcrXtalVec(Event::GcrXtalVec gcrXtalVec){m_gcrXtalVec = gcrXtalVec;}

  ///loads Cal dimensions parameters: 
  StatusCode readGlastDet();
  
  ///gets the First McParticle launched for current event
  Event::McParticle* GcrReconTool::findFirstMcParticle();
  
  ///builds a 3D table of dimensions NTOw,NLAY,NCOL
  ///containing the central point for each Xtal
  void buildXtalCentersTab();
  void verifyXtalCentersTab();

  
  float GcrReconTool::getCALEnergyRaw();
  
  ///builds a 3D table of dimensions NbTow,NbLay,NbCol 
  ///containing a sequential number when a corresponding CalXtalRecData is found.
  ///If no corresponding CalXtalRecData, entry contains -1  
  ///@return number of hits found
  void GcrReconTool::buildHitsMap();
  
  ///clears GcrXtalVec19-09-2006
  void clearGcrXtalVec (Event::GcrXtalVec  *gcrXtalVec);
  
  ///builds GcrXtalsVec, a collection of elements GcrXtals.  References all Xtals
  /// that should theoretically (MonteCarlo) have been touched by first MC particle 
  /// associated to current Event. 
  void buildGcrXtalsVec();
  
  ///stores m_gcrXtalCol into the TDS, using the Event/Recon/CalRecon/GcrReconClasses information
  StatusCode  storeGcrXtals();
  ///stores m_gcrTrack into the TDS, using the Event/Recon/CalRecon/GcrReconClasses information  
  StatusCode storeGcrTrack (); 
  ///stores m_statusWord into the TDS, using the Event/Recon/CalRecon/GcrReconClasses information -- To be written! 
  StatusCode storeGcrReconVals (); 
  
  ///verifies reliability of informations in m_gcrXtalVec 
  void verifyGcrXtalsVec();
  ///verifies reliability of informations in m_hitsMap 
  void verifyHitsMap();
  
  //builds an integer containing information on crossed Xtal Faces as powers of 2
  void getCrossedFaces(int ilay, Point xtalCentralPoint, Point entryPoint, Point exitPoint, int &crossedFaces, double &minDist);
  int getClosestFace(int ilay, Point xtalCentralPoint, Point point, double& minDist);


  /// PRIVATE DATA MEMBERS
  static const float CALRawE_TH=15.0;  //Energy threshold to apply to OBF Gamma and HFC filters cuts
  
  static const int NTOW = 16;
  static const int NLAY = 8;
  static const int NCOL = 12;
  
  /// MsgStream member variable to speed up execution
  MsgStream          m_log;

  /// Pointer to the Gaudi data provider service
  DataSvc*           m_dataSvc;
    
  /// Pointer to the Geant4 propagator
  IPropagator * m_G4PropTool; 
    
  /// the GlastDetSvc used for access to detector info
  IGlastDetSvc*      m_detSvc;
    
  /// TkrGeometrySvc used for access to tracker geometry info
  ITkrGeometrySvc*   m_geoSvc;
  
  //data structures:
  
  /// collection of GcrXtals. Each GcrXtal contains a pointer
  /// to corresponding CalXtalRecData, if any, and the value
  /// of the pathlength defined by first MCParticle direction in
  /// the Xtal. Pointer points
  /// to null if no CalXtalRecData is found.
  Event::GcrXtalVec m_gcrXtalVec;     
  
  /// collection of GcrXtals, used to store info into TDS
  Event::GcrXtalCol*   m_gcrXtalCol; 

  /// GcrTrack, used to store info into TDS
  Event::GcrTrack*   m_gcrTrack; 
  
  /// Hits (Xtals with energy deposit) Map, pointer points on null if no XtalRecData for Xtal  
  Event::CalXtalRecData*  m_hitsMap[NTOW][NLAY][NCOL]; 
  
  /// Xtals Map, contains -1 if no energy deposit is expected from MCparticle trajectory extrapolation
  int  m_gcrXtalsMap[NTOW][NLAY][NCOL]; 
  
  /// number of Xtals touched by McFirstParticle trajectory extrapolation 
  int m_numGcrXtals; 
  
  /// xtal Centers Map
  Point m_xtalCentersTab[NTOW][NLAY][NCOL];
  
  /// variables defining the path of MCFirstParticle:  should these variables be stored in TDS in an Event::GcrTrack for exemple?
  Vector m_mcDir,m_initDir;
  Point m_initPos;
  Point m_calEntryPoint;
  Point m_calExitPoint;
  
  //Cal dimensions parameters: 
  double             m_calZTop;
  double             m_calZBot;
  double             m_calXLo;
  double             m_calXHi;
  double             m_calYLo;
  double             m_calYHi;
  //Xtal dimensions:
  double             m_xtalHeight,m_xtalWidth,m_xtalLength;
  //Other Xtal quantities:
  int m_xNum;    ///< x tower number
  int m_eLATTowers; ///< the value of fLATObjects field, defining LAT towers 
  int m_eTowerCAL;  ///< the value of fTowerObject field, defining calorimeter module 
  int m_eXtal;  ///< the value of fCellCmp field defining CsI crystal
  int m_nCsISeg;  ///< number of geometric segments per Xtal

  // variable that contains information of OBFFilters (Gamma, HFC, DGN, Mip) vetos 
  
  unsigned int m_gcrOBFStatusWord; 
  Event::GcrReconVals*   m_gcrReconVals; 

    //variable that indicates if we want to keep mcTrack direction or TrackReconTrack direction
  std::string m_propertyDir;


} ;

//-----------------------------------------------------------------------------------------------------------------
static ToolFactory<GcrReconTool> s_factory;
const IToolFactory& GcrReconToolFactory = s_factory;

//-----------------------------------------------------------------------------------------------------------------
GcrReconTool::GcrReconTool(const std::string & type, 
			   const std::string & nameIn,
			   const IInterface * parent ) : AlgTool( type, nameIn, parent ),
							 m_log(msgSvc(), name())
{ 
  declareInterface<IGcrReconTool>(this) ; 

  
  return;
}
    
//-----------------------------------------------------------------------------------------------------------------
StatusCode GcrReconTool::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;

  m_log.setLevel(outputLevel());

  //m_log << MSG::INFO << "GcrReconTool BEGIN initialize()" << endreq ;
   
  //Locate and store a pointer to the data service
  IService* iService = 0;
  if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    throw GaudiException("Service [EventDataSvc] not found", name(), sc);
  m_dataSvc = dynamic_cast<DataSvc*>(iService);

  // find TkrGeometrySvc service
  if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
    m_log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
    return StatusCode::FAILURE;
  }
    
  // find G4 propagation tool
  if(!toolSvc()->retrieveTool("G4PropagationTool", m_G4PropTool)) {
    m_log << MSG::ERROR << "Couldn't find the G4PropationTool!" << endreq;
    return StatusCode::FAILURE;
  }

  // find GlastDevSvc service
  if (service("GlastDetSvc", m_detSvc, true).isFailure()){
    m_log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
    return StatusCode::FAILURE;
  }

  // read geometry
  readGlastDet();

  //m_log << MSG::INFO << "GcrReconTool END initialize()" << endreq ;  
    
  

  return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------

/// clears variable gcrXtalVec given in parameter.  This should always correspond to m_gcrXtalVec member variable
void GcrReconTool::clearGcrXtalVec(Event::GcrXtalVec *gcrXtalVec)
{
  //    m_log << MSG::DEBUG << "clearGcrXtalVec in GcrReconTool" << endreq;
    
  gcrXtalVec->clear();
    
  //    m_log << MSG::DEBUG << "clearGcrXtalVec -end-" << endreq;
  return;
}


//-----------------------------------------------------------------------------------------------------------------
StatusCode GcrReconTool::readGlastDet()
{
  StatusCode sc = StatusCode::SUCCESS;
  //m_log << MSG::INFO << "GcrReconTool BEGIN readGlastDet()" << endreq ;  
 
  m_calZTop = m_geoSvc->calZTop();
  m_calZBot = m_geoSvc->calZBot();
  
  double towerPitch = m_geoSvc->towerPitch();
  int xNum = m_geoSvc->numXTowers();
  int yNum = m_geoSvc->numYTowers();
  double calXWidth = m_geoSvc->calXWidth();
  double calYWidth = m_geoSvc->calYWidth();
  double deltaX = 0.5*(xNum*towerPitch - calXWidth);
  double deltaY = 0.5*(yNum*towerPitch - calYWidth);


  m_calXLo = m_geoSvc->getLATLimit(0,LOW)  + deltaX;
  m_calXHi = m_geoSvc->getLATLimit(0,HIGH) - deltaX;
  m_calYLo = m_geoSvc->getLATLimit(1,LOW)  + deltaY;
  m_calYHi = m_geoSvc->getLATLimit(1,HIGH) - deltaY;
  
  //Xtal dimensions:
  m_detSvc->getNumericConstByName("CsIHeight", &m_xtalHeight);  
  m_detSvc->getNumericConstByName("CsIWidth", &m_xtalWidth);  
  m_detSvc->getNumericConstByName("CsILength", &m_xtalLength);  
  m_log << MSG::INFO << "CsIHeight,CsIWidth,CsILength=" << m_xtalHeight << "," << m_xtalWidth<< "," << m_xtalLength<< endreq ;  


  m_detSvc->getNumericConstByName("xNum"	 ,&m_xNum      );
  m_detSvc->getNumericConstByName("eLATTowers",&m_eLATTowers);
  m_detSvc->getNumericConstByName("eTowerCAL" ,&m_eTowerCAL );
  m_detSvc->getNumericConstByName("eXtal"	 ,&m_eXtal  );
 m_detSvc->getNumericConstByName("nCsISeg"	 ,&m_nCsISeg   );

     // m_calZTop,m_calZBot= -48.12,-217.47
     // m_calXLo,m_calXHi=  -728.22, 728.22
     // m_calYLo,m_calYHi=  -728.22, 728.22 

  buildXtalCentersTab();
  //DEBUG:verifyXtalCentersTab();


  //m_log << MSG::INFO << "GcrReconTool END readGlastDet()" << endreq ;  
  return sc;
}

bool GcrReconTool::TriggerEngine4ON(){
    //New definition of GltWord (GRv13r5p5):
    //  b_ROI =      1,  ///>  throttle
    //  b_Track=     2,  ///>  3 consecutive x-y layers hit
    //  b_LO_CAL=    4,  ///>  single log above low threshold
    //  b_HI_CAL=    8,  ///>  single log above high threshold
    //  b_ACDH =    16,  ///>  cover or side veto, high threshold ("CNO")
    //  b_ACDL =    64,  ///>  set if cover or side veto, low threshold

    // Trigger engine 4 (CNOTrigger here): 
    // Track, LoCal, ACDH("CNO"), THROTTLE ("ROI")
    
    // Recover EventHeader Pointer
    SmartDataPtr<Event::EventHeader> pEvent(m_dataSvc, EventModel::EventHeader);
    unsigned int word2 =  ( pEvent==0? 0 : pEvent->triggerWordTwo());

    unsigned int Trig_gemengine = ((word2 >> enums::ENGINE_offset) & enums::ENGINE_mask);
    
    bool engine4ON2 = (Trig_gemengine==4);

    //m_log << MSG::INFO << "engine4ON2= " << engine4ON2 << endreq;
    
    return engine4ON2;
}

float GcrReconTool::getCALEnergyRaw()
{
    StatusCode sc = StatusCode::SUCCESS;
    float CAL_EnergyRaw=0.0;
    // Recover pointer to CalClusters
    SmartDataPtr<Event::CalClusterCol> pCals(m_dataSvc,EventModel::CalRecon::CalClusterCol);
    //Make sure we have valid cluster data and some energy
    if (!pCals) return -1.0;
    if (pCals->empty()) return -2.0;
    Event::CalCluster* calCluster = pCals->front();
    CAL_EnergyRaw  = calCluster->getCalParams().getEnergy();
    if(CAL_EnergyRaw<1.0) return -3.0;

    m_log << MSG::DEBUG << "CAL_EnergyRaw=" << CAL_EnergyRaw << endreq;
    
    return CAL_EnergyRaw;
    
}


bool GcrReconTool::checkFilters(){
  m_log<<MSG::DEBUG<<"GcrReconTool::checkFilters Begin"<<endreq ;
  
  // determine if HFC, DGN, MIP vetos are set.  Returns true if any of them is NOT set, false if the event is vetoed by all these filters.
  // Transfert of statusWord to TDS needs still to be added at the end of this method 

  unsigned int filtersbGamma, filtersbHFC, filtersbDGN, filtersbMip;
  
  // get calEnergyRaw from TDS, to be considered for OBFGamma and OBFHFC cuts
  float calEnergyRaw =getCALEnergyRaw();

  bool debug=true; 
  
  int statusBytes =0; 
  if(debug)
    std::cout << "statusBytes=" << std::hex << statusBytes << std::dec << std::endl;

  SmartDataPtr<OnboardFilterTds::ObfFilterStatus> obfStatus(m_dataSvc, "/Event/Filter/ObfFilterStatus");
  if (obfStatus)
   {
       // Pointer to our retrieved objects
       const OnboardFilterTds::IObfStatus* obfResultGamma = 0;
       const OnboardFilterTds::IObfStatus* obfResultHFC = 0;
       const OnboardFilterTds::IObfStatus* obfResultDGN = 0;
       const OnboardFilterTds::IObfStatus* obfResultMip = 0;
       


       obfResultGamma = obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter);      
       if(obfResultGamma){
	 
	 filtersbGamma = obfResultGamma->getFiltersb() >> 4;
         statusBytes |= (obfResultGamma->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::GammaFilter);

	 if(debug) {
	   std::cout << "filtersbGamma=" << std::hex << filtersbGamma << std::dec << std::endl;
         }	 

       }
       else
           m_log << MSG::INFO <<  "no obfResultGAM" << endreq;


       obfResultHFC = obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::HFCFilter);      
       if(obfResultHFC){
	 filtersbHFC = obfResultHFC->getFiltersb() >> 4;;
         statusBytes |= (obfResultHFC->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::HFCFilter);

	 if(debug){
	   std::cout << "filtersbHFC=" << std::hex << filtersbHFC << std::dec << std::endl;
         }

       }
       else
           m_log << MSG::INFO <<  "no obfResultHFC" << endreq;


       obfResultMip = obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::MipFilter);//MipFilter: Filter key in ObfFilterStatus.h      
       if(obfResultMip){
	 filtersbMip = obfResultMip->getFiltersb() >> 4;;
         statusBytes |= (obfResultMip->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::MipFilter);
     
         if(debug) {
	   std::cout << "filtersbMIP=" << std::hex << filtersbMip << std::dec << std::endl;
	 }
       }
       else
           m_log << MSG::INFO <<  "no obfResultMIP" << endreq;


       obfResultDGN = obfStatus->getFilterStatus(OnboardFilterTds::ObfFilterStatus::DFCFilter); //DFCFilter: Filter key in ObfFilterStatus.h
       //DFCFilter is the old name of DgnFilter, as suggested in method initialize in OnboardFilter.cxx	  
       if(obfResultDGN){
	 filtersbDGN = obfResultDGN->getFiltersb() >> 4;;
         statusBytes |= (obfResultDGN->getFiltersb() >> 4) << (4 * OnboardFilterTds::ObfFilterStatus::DFCFilter);

	 if(debug){
	   std::cout << "filtersbDGN=" << std::hex << filtersbDGN << std::dec << std::endl;
         }

       }
       else
           m_log << MSG::INFO <<  "no obfResultDGN" << endreq;



   }
   else
     m_log << MSG::INFO << "no obfStatus"<< endreq;

   //{'Gam':0,'Hfc':1,'Mip':2,'Dfc':3}

   if(debug){
     std::cout << "statusBytes=" << std::hex << statusBytes << std::dec << std::endl;
     int contribSBGamma = (statusBytes >> 4*(OnboardFilterTds::ObfFilterStatus::GammaFilter))&0xF;
     int contribSBHFC = (statusBytes >> 4*(OnboardFilterTds::ObfFilterStatus::HFCFilter))&0xF;
     int contribSBMip = (statusBytes >> 4*(OnboardFilterTds::ObfFilterStatus::MipFilter))&0xF;
     int contribSBDGN = (statusBytes >> 4*(OnboardFilterTds::ObfFilterStatus::DFCFilter))&0xF;
     std::cout << "JC:contribGamma=" << std::hex << contribSBGamma << std::dec << std::endl;
     std::cout << "JC:contribHFC=" << std::hex << contribSBHFC << std::dec << std::endl;
     std::cout << "JC:contribMip=" << std::hex << contribSBMip << std::dec << std::endl;
     std::cout << "JC:contribDGN=" << std::hex << contribSBDGN << std::dec << std::endl;

   }

   bool cutGamma = (filtersbGamma==0) || (filtersbGamma==6);
   bool cutHFC = ((filtersbHFC==0) || (filtersbHFC==6)) && (calEnergyRaw > CALRawE_TH);
   bool cutMip = (filtersbMip==0) || (filtersbMip==6);
   bool cutDGN = ((filtersbDGN==0) || (filtersbDGN==6))  && (calEnergyRaw > CALRawE_TH);

   bool passFilter = cutGamma || cutHFC || cutMip || cutDGN;
   m_log << MSG::DEBUG << "passFilter:" << passFilter << endreq;

   //TDS variable containing the OBF filter flags:
   
   //ce calcul doit reprendre les 4 booleans et les ranger dans le meme ordre que m_statusBytes:
   m_gcrOBFStatusWord=cutGamma<<OnboardFilterTds::ObfFilterStatus::GammaFilter;
   m_gcrOBFStatusWord|=cutHFC<<OnboardFilterTds::ObfFilterStatus::HFCFilter;
   m_gcrOBFStatusWord|=cutMip<<OnboardFilterTds::ObfFilterStatus::MipFilter;
   m_gcrOBFStatusWord|=cutDGN<<OnboardFilterTds::ObfFilterStatus::DFCFilter;
     std::cout << "gcrOBFStatusWord=" << std::hex << m_gcrOBFStatusWord << std::dec << std::endl;
   
   // store gcrOBFStatus Word in TDS:
   storeGcrReconVals();

   m_log<<MSG::DEBUG<<"GcrReconTool::checkFilters End"<<endreq ;

   return(passFilter);

}


//-----------------------------------------------------------------------------------------------------------------

/**
   This method builds the map of the Xtals that are theoretically touched by the original MCParticle.
   The map is obtained by propagating the first MCparticle initial trajectory.
*/
StatusCode GcrReconTool::findGcrXtals(std::string initDir){

      m_log << MSG::DEBUG << "BEGIN findGcrXtals in GcrReconTool" << endreq;
      m_log << MSG::DEBUG << "initDir=" << initDir << endreq;

  StatusCode sc = StatusCode::SUCCESS;
  
  m_propertyDir = initDir;  
  

  
  // Task #1: build Hits (Xtals with energy deposit) Map 
  buildHitsMap();
  
  //DEBUG: verifyHitsMap();
  
  // Task #2: cleanup gcrXtalVec if any residual unreliable information  
  clearGcrXtalVec(&m_gcrXtalVec);

  // Task #3: initialization gcrXtalsMap
  for (int itow=0; itow<NTOW; itow++)
    for (int ilay=0; ilay<NLAY; ilay++)
      for (int icol=0; icol<NCOL; icol++){
	m_gcrXtalsMap[itow][ilay][icol]=-1;
	//log << MSG::DEBUG << "itow=" << itow << ",ilay=" << ilay << ",icol=" << icol << endreq;
		
      }
      
      
  //Task #4: build gcrXtalsVec
 
  buildGcrXtalsVec();
   
  m_log << MSG::INFO << "m_numGcrXtals= " << m_numGcrXtals << endreq;
  m_log << MSG::INFO << "number of entries in m_gcrXtalVec=" << m_gcrXtalVec.size() << endreq;

  //DEBUG: verifyGcrXtalsVec();

  //Task #5: store gcrXtalCol in TDS
  if (m_numGcrXtals>0)
    sc=storeGcrXtals();
    
  storeGcrTrack(); 
	
  //Task #6: clean up gcrXtalVec
  clearGcrXtalVec(&m_gcrXtalVec);

  //m_log << MSG::INFO << "END findGcrXtals in GcrReconTool" << endreq;
  
  return sc;

    
}


//-----------------------------------------------------------------------------------------------------------------
void GcrReconTool::verifyGcrXtalsVec(){
//m_log << MSG::INFO << "BEGIN verifyGcrXtalsVec in GcrReconTool" << endreq;

	int gcrXtalsVecSize = m_gcrXtalVec.size();

	for(int i=0; i<gcrXtalsVecSize; i++){
	   Event::GcrXtal currentGcrXtal = m_gcrXtalVec.at(i);
	   //m_log << MSG::INFO << "currentGcrXtal.getPathLength()= " << currentGcrXtal.getPathLength() << endreq;
	   //m_log << MSG::INFO << "currentGcrXtal.getXtalId()" << currentGcrXtal.getXtalId() << endreq;

	}

//m_log << MSG::INFO << "END verifyGcrXtalsVec in GcrReconTool" << endreq;

}

//-----------------------------------------------------------------------------------------------------------------
void GcrReconTool::buildGcrXtalsVec(){
  
  m_log << MSG::DEBUG << "BEGIN buildGcrXtalsVec in GcrReconTool" << endreq;
  
    
  /**m_log << MSG::INFO << " (m_calZBot,m_calZTop)= " << m_calZBot << "," << m_calZTop << endreq ;
  m_log << MSG::INFO << " (m_calXLo,m_calXHi)=" << m_calXLo << "," << m_calXHi << endreq ;
  m_log << MSG::INFO << " (m_calYLo,m_calYHi)= " << m_calYLo << "," << m_calYHi << endreq ;
  */
  
  //Point exitP;  // exit point from Calorimeter
  Vector totalTrack; // total propagation vector in the Calorimeter

  m_numGcrXtals=0;
    
  //====== Task #4.1: define initial Position and Direction =====================================
  
  //get MC direction:
  
  Event::McParticle* firstMcParticle = GcrReconTool::findFirstMcParticle();

  if(firstMcParticle != 0)
  	m_mcDir = Vector(firstMcParticle->initialFourMomentum().vect().unit());
  else 
        m_mcDir = Vector(-1e9,-1e9,-1e9);
 
 
  if(m_propertyDir=="MC"){ // if taking MCTrack, get theoretical initial position and direction

     // ******the first McParticle initial position (when launched):
      const HepPoint3D& initPosition =firstMcParticle->initialPosition(); //initialPosition of firstMCParticle 
      // casting to Point of initPosition, necessary to define propagator steps

      m_initPos.setX(initPosition.x());
      m_initPos.setY(initPosition.y());
      m_initPos.setZ(initPosition.z());

      // ******the first McParticle initial direction (when launched):
      m_initDir = Vector(m_mcDir.x(), m_mcDir.y(), m_mcDir.z());

  
  }
  
  else if(m_propertyDir=="TKR"){ // if taking MCTrack, get reconstructed initial position and direction
  
	SmartDataPtr<Event::TkrTrackCol>   pTracks(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

	if (pTracks){   
            // Count number of tracks
            int nTracks = pTracks->size();
            if(nTracks < 1) return;

            // Get the first Track - it should be the "Best Track"
            Event::TkrTrackColConPtr pTrack = pTracks->begin();

	    const Event::TkrTrack* track_1 = *pTrack;
             
	    m_initPos = track_1->getInitialPosition();
	    m_initDir = track_1->getInitialDirection().unit();

	}
	else 
    {
      m_log << MSG::INFO << "no TkrTrackCol found" << endreq;
    }
  } else if(m_propertyDir=="MOM")
    {
      MsgStream log(msgSvc(), name());
      StatusCode sc = StatusCode::SUCCESS;
      ITkrFilterTool* filterTool;
      if ((sc = toolSvc()->retrieveTool("TkrFilterTool", filterTool)).isFailure())
        {
          log << MSG::ERROR << " could not find TkrFilterTool " << endreq;
          return sc;
        }
      sc = filterTool->doFilterStep();
      //ok now the results need to be retrieved from the TDS:
      // Recover pointer to TkrEventParams
      Event::TkrEventParams* tkrEventParams = 
        SmartDataPtr<Event::TkrEventParams>(m_dataSvc, EventModel::TkrRecon::TkrEventParams);
      // First pass means no TkrEventParams
      if (tkrEventParams == 0)
        {
          log << MSG::ERROR << " could not find TkrEventParams in the TDS " << endreq;
          return StatusCode::FAILURE;
        }
      m_initPos = tkrEventParams->getEventPosition();
      m_initDir = -tkrEventParams->getEventAxis();    
    } else{
    m_log<<MSG::ERROR<<"Invalid property "<<m_propertyDir<<endreq;
    return;
  }
  
  double x0 = m_initPos.x();
  double y0 = m_initPos.y();
  double z0 = m_initPos.z();

  double ux0 = m_initDir.x();
  double uy0 = m_initDir.y();
  double uz0 = m_initDir.z();
 


    
  m_log << MSG::DEBUG << "initPos="<< '(' << x0 << ',' << y0 << ',' << z0 << ')'<<endreq;
  
  m_log << MSG::DEBUG << "initDir=" << '(' << ux0 << ',' << uy0 << ',' << uz0 << ')'<<endreq ;
  
  // ========================================================================================

    
   
  //====== Task #4.2 & 4.3: find entry and exit point from Calorimeter=========================================
  
  if (uz0!=0)
    {
      m_calEntryPoint = m_initPos+((m_calZTop-z0)/uz0)*m_initDir; 
      m_calExitPoint  = m_initPos+((m_calZBot-z0)/uz0)*m_initDir;        
      
    }
  else
    {
      if (ux0!=0)
        {
         
	  if (ux0>0){
	    m_calEntryPoint = m_initPos+((m_calXLo-x0)/ux0)*m_initDir; 
	    m_calExitPoint  = m_initPos+((m_calXHi-x0)/ux0)*m_initDir;
	  }
	  else{
	    m_calEntryPoint = m_initPos+((m_calXHi-x0)/ux0)*m_initDir; 
	    m_calExitPoint  = m_initPos+((m_calXLo-x0)/ux0)*m_initDir;
	    }
        }
      else if (uy0!=0)
        {
	  if (uy0>0){
	    m_calEntryPoint = m_initPos+((m_calYLo-y0)/uy0)*m_initDir; 
	    m_calExitPoint  = m_initPos+((m_calYHi-y0)/uy0)*m_initDir;
	    }
	  else{
	    m_calEntryPoint = m_initPos+((m_calYHi-y0)/uy0)*m_initDir; 
	    m_calExitPoint  = m_initPos+((m_calYLo-y0)/uy0)*m_initDir;
	    
	    }
        }
  
    }
    
  //m_log << MSG::INFO << "CalEntryPoint=       " << '(' << m_calEntryPoint.x() << ',' << m_calEntryPoint.y() << ',' << m_calEntryPoint.z() << ')'  <<endreq;
    
  //m_log << MSG::INFO << "CalExitPoint=       " << '(' << m_calExitPoint.x() << ',' << m_calExitPoint.y() << ',' << m_calExitPoint.z() << ')'  <<endreq;
  //exitP.printOn(std::cout); 

  // ========================================================================================
     
  //====== Task #4.4:  propagate and update m_gcrXtalsMap & gcrXtalVec:
  totalTrack=m_calEntryPoint-m_calExitPoint;

  double totalLength = sqrt(totalTrack*totalTrack); 
  m_log << MSG::DEBUG << "totalLength=" << totalLength <<  " ,XtalHeigth, XtalLength=" << m_xtalHeight << ","<< m_xtalLength<<endreq;

  if((totalLength<m_xtalHeight) || (totalLength>4*m_xtalLength))
    {
      m_log << MSG::DEBUG << "totalLength condition not fullfilled."<<endreq;
      return;
    }
  
  m_G4PropTool->setStepStart(m_calEntryPoint,m_initDir);
    
  //m_log << MSG::INFO << "\n totalTrack=" <<endreq;
  //totalTrack.printOn(std::cout);
    
  m_G4PropTool->step(totalLength);     
    
  int numSteps = m_G4PropTool->getNumberSteps();
    
  m_log << MSG::DEBUG << "number of Propagation Steps= " << numSteps <<endreq;
    
  idents::VolumeIdentifier volId;
  idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
  int itow=-1;
  int ilay=-1;
  int icol=-1;
  double arcLen_step=0.;
  Point stepPos,entryPoint,exitPoint;
    
 
  for(int istep=0; istep<numSteps; ++istep)
    {
      m_log << MSG::DEBUG << "G4 propagation step number= " << istep <<endreq;
      volId=m_G4PropTool->getStepVolumeId(istep);
      //m_log << MSG::INFO << "volId= " << volId.name() <<endreq;
      
      volId.prepend(prefix);
      
      arcLen_step = m_G4PropTool->getStepArcLen(istep);
      
      //Comment from Tracy Usher: 
      //getStepPosition(*int* stepIdx = -1) *const* = 0 actually returns the position at the end of a given step.
      entryPoint = m_G4PropTool->getStepPosition(istep-1);  
      exitPoint = m_G4PropTool->getStepPosition(istep);
      
      
      int crossedFaces;
      double minDist;
      
      
      
      // ********* build GcrXtalsMap and gcrXtalVec
      if(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0)// in Xtal ? 
        {
          itow=4*volId[1]+volId[2];
          ilay=volId[4];
          icol=volId[6];
          
          Point xtalCentralPoint = m_xtalCentersTab[itow][ilay][icol];
          
          m_log << MSG::DEBUG << "(itow,ilay,icol)=(" << itow << "," << ilay << "," << icol << ")" << endreq;
          m_log << MSG::DEBUG << "istep " << istep <<  ": entryPoint=" << entryPoint <<endreq;
          m_log << MSG::DEBUG << "istep " << istep <<  ": exitPoint= " << exitPoint << endreq;
          m_log << MSG::DEBUG << "arcLen_step, exitPoint - entryPoint " << arcLen_step << ","<<  sqrt((exitPoint - entryPoint) * (exitPoint - entryPoint)) << endreq;
          
          
          if(m_gcrXtalsMap[itow][ilay][icol]<0)
            { //if Xtal has not already been identified as part of the theoretical track
              //m_log << MSG::INFO << "new touched Xtal"<< endreq;
              
              Event::CalXtalRecData* xTalData = m_hitsMap[itow][ilay][icol];   // search for corresponding hits entry, if any
              
              if(xTalData)
                {//add and update entry on grcXtalVec, only if a corresponding CalXtalRecData is found
                  idents::CalXtalId xtalId = xTalData->getPackedId();
                  m_gcrXtalsMap[itow][ilay][icol] = m_numGcrXtals++; // the entry of m_gcrXtalsMap corresponds to m_gcrXtalVec index
                  m_gcrXtalVec.push_back(Event::GcrXtal());
                  Event::GcrXtal& gcrXtal = m_gcrXtalVec.back();
                  gcrXtal.setXtalId(xtalId);  
                  gcrXtal.setPathLength(arcLen_step);	 
                  gcrXtal.setEntryPoint(entryPoint);
                  gcrXtal.setExitPoint(exitPoint);
                  getCrossedFaces(ilay,xtalCentralPoint,entryPoint,exitPoint,crossedFaces,minDist);
                  gcrXtal.setCrossedFaces(crossedFaces);
                  gcrXtal.setClosestFaceDist(minDist);
                  
                  //TEST:
                  int n,m;
                  gcrXtal.getReadableXedFaces(crossedFaces,n,m);
                  m_log <<MSG::DEBUG<<"crossedFaces,n,m="<< crossedFaces << "," << n << "," << m << endreq;
                  m_log <<MSG::DEBUG<<"gcrXtalId,pathLength"<< gcrXtal.getXtalId() << "," << gcrXtal.getPathLength() << endreq;
                } else 
                { 
                  m_log<<MSG::DEBUG<< "track crosses crystal, but no corresponding CalXtal found"<<endreq;
                  m_gcrXtalsMap[itow][ilay][icol]=-2;  
                }

	   
	   // else
	    //     m_log << MSG::INFO << "new touched Xtal: no corresponding CalXtalRecData was found"<< endreq;
		 
            //m_log << MSG::INFO << "m_gcrXtalVec.size()= "<< m_gcrXtalVec.size() << endreq; 
	    
	  } else
            {// if Xtal has already been identified as part of the theoretical track
              //m_log << MSG::INFO << "Xtal already touched"<< endreq;
              
              int gcrXtalVecIndex = m_gcrXtalsMap[itow][ilay][icol];  // we get the index of gcrXtalVec entry
              //m_log << MSG::INFO << "gcrXtalVecIndex =" << gcrXtalVecIndex << endreq;
              // m_log << MSG::INFO << "m_gcrXtalVec.size()= "<< m_gcrXtalVec.size() << endreq; 
              
              if(m_hitsMap[itow][ilay][icol]){  // if a corresponding CalXtalRecData is found
                
                Event::GcrXtal currentGcrXtal = m_gcrXtalVec.at(gcrXtalVecIndex); // find corresponding entry of gcrXtalVec
                
                currentGcrXtal.setPathLength(currentGcrXtal.getPathLength() + arcLen_step); // pathLength update
                currentGcrXtal.setExitPoint(exitPoint);// exitPoint update
                getCrossedFaces(ilay,xtalCentralPoint,entryPoint,exitPoint,crossedFaces,minDist);
                currentGcrXtal.setCrossedFaces(crossedFaces);
                currentGcrXtal.setClosestFaceDist(minDist);
                
                m_gcrXtalVec.at(gcrXtalVecIndex) = currentGcrXtal;
                
                
                //TEST:
                m_log << MSG::DEBUG << " Xtal has already been identified as part of the theoretical track"<< endreq; 
                //int n,m;
                //currentGcrXtal.getReadableXedFaces(crossedFaces,n,m);
                m_log <<MSG::DEBUG<<"gcrXtalId,pathLength= "<< currentGcrXtal.getXtalId() << "," << currentGcrXtal.getPathLength() << endreq;
                m_log <<MSG::DEBUG<<"m_gcrXtalVec.at(gcrXtalVecIndex), gcrXtalId,pathLength= "<< m_gcrXtalVec.at(gcrXtalVecIndex).getXtalId() << "," << m_gcrXtalVec.at(gcrXtalVecIndex).getPathLength() << endreq;
                
              }
              //else     
              //m_log << MSG::INFO << "Xtal already touched: no corresponding CalXtalRecData was found"<< endreq; 
              
            }
          

        }else // end of: if in Xtal
        {
          m_log <<MSG::DEBUG<<"Volume Id not in crystal"<<endreq;
        }
    }// end of: for(int istep= ..)
  
  // ========================================================================================

  m_log << MSG::DEBUG << "END buildGcrXtalsVec in GcrReconTool"<< endreq;
}


//-----------------------------------------------------------------------------------------------------------------

Event::McParticle* GcrReconTool::findFirstMcParticle(){

  //m_log << MSG::INFO << "BEGIN findFirstMcParticle in GcrReconTool" << endreq;

  Event::McParticleCol* mcParticleCol = SmartDataPtr<Event::McParticleCol>(m_dataSvc, EventModel::MC::McParticleCol);
  // Event::McParticleCol::const_iterator mcParticleColIter=mcParticleCol->begin(); 

   
  Event::McParticle* firstMcParticle = 0; 
  if (mcParticleCol){  
    firstMcParticle = *mcParticleCol->begin();
    //m_log << MSG::INFO << "firstMcParticle particle ID= " << firstMcParticle->particleProperty()<< endreq;
  } 
  else 
    m_log << MSG::INFO << "no mcParticleCol" << endreq;
       
  //m_log << MSG::INFO << "firstMcParticle.x()= " << firstMcParticle->initialPosition().x() << endreq;
 
 
 
  //m_log << MSG::INFO << "END findFirstMcParticle in GcrReconTool" << endreq;
  
  return firstMcParticle;
   


}




//-----------------------------------------------------------------------------------------------------------------
/**
 *  This method builts a 3D map[NTOW][NLAY][NCOL], whose entries are pointers to XTalRecData, if any
 */
void GcrReconTool::buildHitsMap()
{ 

  //find Event info:

  Event::EventHeader* eventHeader = SmartDataPtr<Event::EventHeader>(m_dataSvc, EventModel::EventHeader); 
 
  m_log << MSG::INFO << "event= " << eventHeader->event() << endreq;

  //m_log << MSG::INFO << "BEGIN buildHitsMap in GcrReconTool" << endreq;
 
  int numHits=0;
  int itow, ilay, icol;
 
  //initialization of Xtal(Hits) Ids 
  for (itow=0; itow<NTOW; itow++)
    for (ilay=0; ilay<NLAY; ilay++)
      for (icol=0; icol<NCOL; icol++){
	m_hitsMap[itow][ilay][icol]=0;		
      }



  Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(m_dataSvc, EventModel::CalRecon::CalXtalRecCol); 

   
  if (calXtalRecCol){
    // loop over CalXtalRecdata
    //m_log << MSG::DEBUG << "juste before loop, calXtalRecCol->numberOfObjects()= " << calXtalRecCol->numberOfObjects() << endreq;
    for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
      {

	Event::CalXtalRecData* xTalData = *xTalIter;
           
	itow=xTalData->getPackedId().getTower();
	ilay=xTalData->getPackedId().getLayer();
	icol=xTalData->getPackedId().getColumn();
	    
	m_hitsMap[itow][ilay][icol] = xTalData;
	    
	// m_log << MSG::INFO << "m_hitsMap[itow][ilay][icol]->getPosition()= " << endreq;
	//(m_hitsMap[itow][ilay][icol])->getPosition().printOn(std::cout);
	// m_log << MSG::INFO << "\n" << endreq;
	    
	numHits++;
	  
      }
  }
  else 
    m_log << MSG::INFO << "no calXtalRecCol " << endreq;
	
    

  //m_log << MSG::INFO << "END buildHitsMap in GcrReconTool, numHits= " << numHits << endreq;
  
}

// ----------------------------------------------------------------------------
void GcrReconTool::buildXtalCentersTab(){

    //constants defining the position of the fields in VolumeIdentifier 
    enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
            fMeasure, fCALXtal,fCellCmp, fSegment};  


 for (int itow=0; itow<NTOW; itow++)
    for (int ilay=0; ilay<NLAY; ilay++)
      for (int icol=0; icol<NCOL; icol++){
      
      	idents::VolumeIdentifier segm0Id;
	segm0Id.append(m_eLATTowers);
	segm0Id.append(itow/m_xNum);
	segm0Id.append(itow%m_xNum);
	segm0Id.append(m_eTowerCAL);
	segm0Id.append(ilay);
	segm0Id.append(ilay%2); 
        segm0Id.append(icol);
	segm0Id.append(m_eXtal);
	segm0Id.append(0);
	
	HepTransform3D transf;

        //get 3D transformation for segment 0 of this crystal
	m_detSvc->getTransform3DByID(segm0Id,&transf);

        //get position of the center of the segment 0
        Vector vect0 = transf.getTranslation();

        // create Volume Identifier for the last segment of this crystal
        idents::VolumeIdentifier segm11Id;

        // copy all fields from segm0Id, except segment number
	for(int ifield = 0; ifield<fSegment; ifield++)segm11Id.append(segm0Id[ifield]);
	segm11Id.append(m_nCsISeg-1); // set segment number for the last segment

	//get 3D transformation for the last segment of this crystal
        m_detSvc->getTransform3DByID(segm11Id,&transf);

        //get position of the center of the last segment
	Vector vect11 = transf.getTranslation();

	Point p0(0.,0.,0.);		

        // position of the crystal center
	Point pCenter = p0+(vect0+vect11)*0.5; 

	m_xtalCentersTab[itow][ilay][icol]=pCenter;		
      }
      


}

void GcrReconTool::verifyXtalCentersTab(){

    for (int itow=0; itow<NTOW; itow++)
	for (int ilay=0; ilay<NLAY; ilay++)
	  for (int icol=0; icol<NCOL; icol++){
             Point centralPoint = m_xtalCentersTab[itow][ilay][icol];
	     m_log << MSG::INFO << "Xtal["<<itow<<","<<ilay << "," << icol <<"],CentralPoint="<< centralPoint<< endreq;
	  }


}

// ----------------------------------------------------------------------------
void GcrReconTool::verifyHitsMap(){

//m_log << MSG::INFO << "BEGIN verifyHitsMap in GcrReconTool"<< endreq;

for (int itow=0; itow<NTOW; itow++)
    for (int ilay=0; ilay<NLAY; ilay++)
      for (int icol=0; icol<NCOL; icol++){
        Event::CalXtalRecData* xTalData = m_hitsMap[itow][ilay][icol];
	if(xTalData)
          m_log << MSG::INFO << "m_hitsMap[" << itow << "][" << ilay << "][" << icol << "]->getPackedId()= " << xTalData->getPackedId()<< endreq;
			
      }

  //m_log << MSG::INFO << "END verifyHitsMap in GcrReconTool"<< endreq;
}

// ----------------------------------------------------------------------------

/**
 * @author CL 06/02/2006
 * This method allows to store GcrReconVals in TDS structure
 */
StatusCode GcrReconTool::storeGcrReconVals () {
  StatusCode sc = StatusCode::SUCCESS;
 
  
  m_log << MSG::DEBUG << "BEGIN storeGcrSelectVals in GcrSelectTool" << endreq;
  
  m_gcrReconVals = SmartDataPtr<Event::GcrReconVals>(m_dataSvc,EventModel::CalRecon::GcrReconVals);

  // If no pointer then create it
  if (m_gcrReconVals == 0)
  {
    m_gcrReconVals = new Event::GcrReconVals();
    sc = m_dataSvc->registerObject(EventModel::CalRecon::GcrReconVals, m_gcrReconVals);
    if (sc.isFailure()) throw GaudiException("Failed to create GCR ReconVals !", name(), sc);
  }


  m_log << MSG::DEBUG << "In storeGcrReconVals, m_gcrOBFStatusWord =" << m_gcrOBFStatusWord << endreq;

  m_gcrReconVals->setGcrOBFStatusWord(m_gcrOBFStatusWord);
  

  m_log << MSG::DEBUG << "END storeGcrReconVals in GcrSelectTool" << endreq;
 
  return sc;

  
}


// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

/**
 * @author CL 06/02/2006
 * This method allows to store GcrXtals in TDS structure
 */
StatusCode GcrReconTool::storeGcrXtals () {
  StatusCode sc = StatusCode::SUCCESS;

  //m_log << MSG::INFO << "BEGIN storeGcrXtals in GcrReconAlg" << endreq;

  m_gcrXtalCol = SmartDataPtr<Event::GcrXtalCol>(m_dataSvc,EventModel::CalRecon::GcrXtalCol);

  // If no pointer then create it
  if (m_gcrXtalCol == 0)
    {
      m_gcrXtalCol = new Event::GcrXtalCol();
      sc = m_dataSvc->registerObject(EventModel::CalRecon::GcrXtalCol, m_gcrXtalCol);
      if (sc.isFailure()) throw GaudiException("Failed to create Gcr Xtal Collection!", name(), sc);
    }

  int nbStoredGcrXtals=0;

  for(Event::GcrXtalVec::iterator gcrXtalIter=m_gcrXtalVec.begin(); gcrXtalIter != m_gcrXtalVec.end(); gcrXtalIter++)
    {
      Event::GcrXtal& gcrXtal = *gcrXtalIter;

      // Need to create a new GcrXtal which will be "owned" by the TDS
      Event::GcrXtal* newGcrXtal = new Event::GcrXtal();

      // Now copy to this new GcrXtal. This should properly copy all elements (including inherited vector)
      *newGcrXtal = gcrXtal;
       
	   m_log << MSG::DEBUG << "in storeGcrXtals, gcrXtalId,pathLength= " << gcrXtal.getXtalId() << "," << gcrXtal.getPathLength() << endreq;

      // Store in collection (and reliquish ownership of the object)
      m_gcrXtalCol->push_back(newGcrXtal); 
      nbStoredGcrXtals++;
    }

  //m_log << MSG::INFO << "END storeGcrXtals in GcrReconAlg" << endreq;
  m_log << MSG::INFO << "nbStoredGcrXtals=" << nbStoredGcrXtals << endreq;

  return sc;

  
}


// ----------------------------------------------------------------------------

/**
 * @author CL 06/02/2006
 * This method allows to store GcrTrack in TDS structure
 */
StatusCode GcrReconTool::storeGcrTrack () {

  StatusCode sc = StatusCode::SUCCESS;
 
  m_log << MSG::INFO << "BEGIN storeGcrTrack in GcrReconTool" << endreq;
  
 /**m_log << MSG::INFO << "m_initDir=" << m_initDir << endreq;
m_log << MSG::INFO << "m_calEntryPoint=" << m_calEntryPoint << endreq;
m_log << MSG::INFO << "m_calExitPoint=" << m_calExitPoint << endreq;
*/

  m_gcrTrack = SmartDataPtr<Event::GcrTrack>(m_dataSvc,EventModel::CalRecon::GcrTrack);

  // If no pointer then create it
  if (m_gcrTrack == 0)
    {
      m_gcrTrack = new Event::GcrTrack();
      sc = m_dataSvc->registerObject(EventModel::CalRecon::GcrTrack, m_gcrTrack);
      if (sc.isFailure()) throw GaudiException("Failed to create Gcr Track !", name(), sc);
    }

//  SETTING THE MC DIRECTION:

  //m_log << MSG::INFO << "m_mcDir=" << m_initDir << endreq;

  m_gcrTrack->setMcDir(m_mcDir);
  m_gcrTrack->setCalEntryPoint(m_calEntryPoint);
  m_gcrTrack->setCalExitPoint(m_calExitPoint);
  m_gcrTrack->setDirection(m_initDir);
  
  if(m_propertyDir=="MC"){// if keeping MC Track
	  m_log << MSG::INFO << "keeping MC Track" << endreq;	  
      m_gcrTrack->setDirError(Vector(0.0,0.0,0.0));
  }
  else{// if we are keeping the Tracker Track
    //SETTING RECONSTRUCTED TKR TRACK INFO:
	m_log << MSG::INFO << "keeping Tracker Recon Track" << endreq;


	SmartDataPtr<Event::TkrTrackCol>   pTracks(m_dataSvc,EventModel::TkrRecon::TkrTrackCol);

	if (pTracks){   
            // Count number of tracks
            int nTracks = pTracks->size();
            if(nTracks < 1) return sc;

            // Get the first Track - it should be the "Best Track"
            Event::TkrTrackColConPtr pTrack = pTracks->begin();

	    const Event::TkrTrack* track_1 = *pTrack;

		m_log << MSG::DEBUG << "m_TkrTrackDir=" << track_1->getInitialDirection() << endreq;
		m_log << MSG::DEBUG << "m_KalmanThetaMS=" << track_1->getKalThetaMS() << endreq;

	    m_gcrTrack->setDirError(Vector(0.0,0.0,cos(track_1->getKalThetaMS())));

	}
	else
            m_log << MSG::INFO << "no TkrTrackCol found" << endreq;
	
    }// end of "if we are keeping the Tracker Track"


  m_log << MSG::INFO << "END storeGcrTrack in GcrReconAlg" << endreq;
 
  return sc;

  
}

void GcrReconTool::getCrossedFaces(int ilay, Point xtalCentralPoint, Point entryPoint, Point exitPoint, int &crossedFaces, double &minDist){
  m_log << MSG::DEBUG << "GcrReconTool::getCrossedFaces BEGIN" << endreq;
  
  double d0,d1,dm; 
  
  m_log << MSG::DEBUG << "@@@@@ entryPoint: " << endreq;

  int face0 = getClosestFace(ilay,xtalCentralPoint,entryPoint,d0);

  m_log << MSG::DEBUG << "@@@@@ exitPoint: " << endreq;

  int face1 = getClosestFace(ilay,xtalCentralPoint,exitPoint,d1);

  m_log << MSG::DEBUG << "@@@@@ middlePoint: " << endreq;

  Point middlePoint = Point(0.5*(entryPoint.x()+exitPoint.x()),0.5*(entryPoint.y()+exitPoint.y()),0.5*(entryPoint.z()+exitPoint.z()));
  /*int facem =*/ getClosestFace(ilay,xtalCentralPoint,middlePoint,dm);
  
  m_log << MSG::DEBUG << "face0,face1=" << face0 <<"," << face1 << endreq;
  
  //  crossedFaces = std::pow(2.,double(face0)) + std::pow(2.,double(face1));
  crossedFaces = static_cast<int>(std::pow(2.,face0) + std::pow(2.,face1));

  m_log << MSG::DEBUG << "crossedFaces=" << crossedFaces << endreq;
 
  minDist=dm;
 
  m_log << MSG::DEBUG << "GcrReconTool::getCrossedFaces END" << endreq;

}

int GcrReconTool::getClosestFace(int ilay, Point xtalCentralPoint, Point point, double& minDist){

    m_log << MSG::DEBUG << "GcrReconTool::getClosestFace BEGIN" << endreq;

    enum{zTopFace,zBotFace,xLeftFace,xRightFace,yLeftFace,yRightFace};
    
    double xc,yc,zc,zTop,zBot,xLeft,xRight,yLeft,yRight, d_zTop, d_zBot, d_xLeft, d_xRight, d_yLeft, d_yRight;

    xc = xtalCentralPoint.x();
    yc = xtalCentralPoint.y();
    zc = xtalCentralPoint.z();

      
    zBot = zc-0.5*m_xtalHeight;
    zTop = zc+0.5*m_xtalHeight;
    
    if (ilay%2 == 0){
	yLeft = yc-0.5*m_xtalWidth;
	yRight = yc+0.5*m_xtalWidth;
	xLeft = xc-0.5*m_xtalLength;
	xRight = xc+0.5*m_xtalLength;
    }
    else {
    	xLeft = xc-0.5*m_xtalWidth;
	xRight = xc+0.5*m_xtalWidth;
	yLeft = yc-0.5*m_xtalLength;
	yRight = yc+0.5*m_xtalLength;
	
    }
    
	m_log << MSG::DEBUG << "Point,zTop,zBot,xLeft,xRight,yLeft,yRight="<< point<<","<<zTop<<","<<zBot<<","<<xLeft<<","<<xRight<<","<<yLeft<<","<<yRight<< endreq;
    
    
    d_zTop = std::fabs(point.z()-zTop);
    d_zBot = std::fabs(point.z()-zBot);
    d_xLeft = std::fabs(point.x()-xLeft);
    d_xRight = std::fabs(point.x()-xRight);
    d_yLeft = std::fabs(point.y()-yLeft);
    d_yRight = std::fabs(point.y()-yRight);

	m_log << MSG::DEBUG << "d_zTop,d_zBot,d_xLeft,d_x_Right,d_yLeft,d_yRight="<<d_zTop<< ","<<d_zBot<<","<<d_xLeft<<","<<d_xRight<<","<<d_yLeft<<","<<d_yRight << endreq;

    
    minDist = std::min(std::min(std::min(std::min(std::min(d_zBot,d_zTop),d_xLeft),d_xRight),d_yLeft),d_yRight);
    
    //m_log << MSG::INFO << "minDist=" << minDist<< endreq;


    if(minDist == d_zTop)
      return zTopFace;
      else if(minDist == d_zBot)
        return zBotFace;
	else if(minDist == d_xLeft)
	  return xLeftFace;
	  else if(minDist == d_xRight)
	    return xRightFace;
	    else if(minDist == d_yLeft)
	      return yLeftFace;
	      else if(minDist == d_yRight)
	        return yRightFace;
    else return -1000;
    
	m_log << MSG::DEBUG << "GcrReconTool::getClosestFace END" << endreq;

}

