// Include files
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/MsgStream.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "geometry/Point.h"
#include "idents/VolumeIdentifier.h"
#include "facilities/Timestamp.h"

#include "CalibData/Cal/CalCalibMuSlope.h"

#include "TestPosTool.h"

static ToolFactory<TestPosTool> s_factory;
const IToolFactory& TestPosToolFactory = s_factory;

TestPosTool::TestPosTool( const std::string& type, 
                          const std::string& name, 
                          const IInterface* parent)
  : AlgTool(type,name,parent) {
  declareInterface<ICalPosTool>(this);

  declareProperty( "startTime",  
                   m_startTimeAsc = "2003-1-10_00:20");
  declareProperty( "calibFlavor", m_calibFlavor = "none");

}

StatusCode TestPosTool::initialize() {
  MsgStream msglog(msgSvc(), name());
  msglog << MSG::INFO << "initialize" << endreq;

  StatusCode sc = StatusCode::SUCCESS;

  // extracting int constants
  double value;  // intermediate variable for reading constants from
  // xml file as doubles and converting them to interger 
  typedef std::map<int*,std::string> PARAMAP;
  PARAMAP param; // map containing pointers to integer constants to be read
  // with their symbolic names from xml file used as a key 

  //     filling the map with information on constants to be read 
    
  param[&m_xNum]=        std::string("xNum");
  param[&m_yNum]=        std::string("yNum");
  param[&m_eTowerCAL]=   std::string("eTowerCAL");
  param[&m_eLATTowers]=  std::string("eLATTowers");
  param[&m_nCsISeg]= std::string("nCsISeg");
  param[&m_eXtal]=       std::string("eXtal");
    
  // now try to find the GlastDevSvc service
    
  //    IGlastDetSvc* detSvc;
  sc = service("GlastDetSvc", detSvc);
    
  // loop over all constants information contained in the map
  for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){

    //  attempt to get teh constant value via the method of GlastDetSvc
    if(!detSvc->getNumericConstByName((*it).second, &value)) {

      // if not successful - give the error message and return
      msglog << MSG::ERROR << " constant " <<(*it).second
             <<" not defined" << endreq;
      return StatusCode::FAILURE;

      //  if successful - fill the constant using the pointer from the map
    } else *((*it).first)= int(value);
  }
        
  // get pointer to CalibDataSvc
  sc = service("CalibDataSvc", m_pCalibDataSvc, true);

  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR 
           << "Could not get IDataProviderSvc interface of CalibDataSvc" 
           << endreq;
    return sc;
  }
 
  // Query the IDetDataSvc interface of the calib data service
  sc = m_pCalibDataSvc->queryInterface(IID_IDetDataSvc, 
                                       (void**) &m_detDataSvc);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR 
           << "Could not query IDetDataSvc interface of CalibDataSvc" 
           << endreq;
    return sc;
  } else {
    msglog << MSG::DEBUG 
           << "Retrieved IDetDataSvc interface of CalibDataSvc" 
           << endreq;
  }

  // extracting double constants
    
  typedef std::map<double*,std::string> DPARAMAP;
  DPARAMAP dparam; // map containing pointers to double constants to be read
  // with their symbolic names from xml file used as a key 

  dparam[&m_lightAtt]=std::string("cal.lightAtt");
  dparam[&m_CsILength]=std::string("CsILength");
  dparam[&m_CsIWidth]=std::string("CsIWidth");
  dparam[&m_CsIHeight]=std::string("CsIHeight");
    
  for(DPARAMAP::iterator dit=dparam.begin(); dit!=dparam.end();dit++){
    if(!detSvc->getNumericConstByName((*dit).second,(*dit).first)) {
      msglog << MSG::ERROR << " constant " <<(*dit).second << " not defined" << endreq;
      return StatusCode::FAILURE;
    } 
  }
       
  // Get properties from the JobOptionsSvc
  sc = setProperties();
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR << "Could not set jobOptions properties" << endreq;
    return sc;
  }

  unsigned int underpos = m_startTimeAsc.find("_");
  if (underpos < m_startTimeAsc.size()) {
    m_startTimeAsc.replace(underpos, 1, " ");
  }
  m_startTime = facilities::Timestamp(m_startTimeAsc).getClibTime();

  msglog << MSG::DEBUG << "Properties were read from jobOptions" << endreq;
  msglog << MSG::INFO << "Time of first event: (ascii) "
         << m_startTimeAsc       << endreq; 
  msglog << MSG::INFO << "Time of first event: (seconds since 1970) "
         << m_startTime       << endreq; 


  // get cal energy tool
  sc = toolSvc()->retrieveTool(std::string("TestEnergyTool"),m_pCalEnergyTool);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create " << "TestEnergyTool" << endreq;
    return sc;
  }

  // set properties of calEnergyTool to match my properties.
  IProperty *property;
  // Query the IProperty interface of the cal energy tool
  sc = m_pCalEnergyTool->queryInterface(IID_IProperty, 
                                        (void**) &property);
  if ( !sc.isSuccess() ) {
    msglog << MSG::ERROR 
           << "Could not query IProperty interface of calenergytool" 
           << endreq;
    return sc;
  } else {
    msglog << MSG::DEBUG 
           << "Retrieved IProperty interface of calenergytool" 
           << endreq;
  }
    
  sc = property->setProperty(std::string("calibFlavor"), m_calibFlavor);
  if (sc.isFailure() ) {
    msglog << MSG::ERROR << "  Unable to create calibFlavor property for calEnergyTool" << endreq;
    return sc;
  }
  
  return StatusCode::SUCCESS;
}


// calculate position given the digital response on both faces
StatusCode TestPosTool::calculate(const idents::CalXtalId &xtalId,
                                  int adcP, 
                                  idents::CalXtalId::AdcRange rangeP,
                                  int adcN, 
                                  idents::CalXtalId::AdcRange rangeN,
                                  float &position                // output
                                  ) {

  MsgStream msglog(msgSvc(), name());

  DataObject *pObject;
  pMuSlopes = 0;
  std::string fullMuSlopePath = "/Calib/CAL_MuSlope/"+m_calibFlavor;
  if (m_pCalibDataSvc->retrieveObject(fullMuSlopePath, pObject) == StatusCode::SUCCESS) {
    pMuSlopes = dynamic_cast<CalibData::CalCalibMuSlope *> (pObject);
    if (!pMuSlopes) {
      msglog << MSG::ERROR << "Dynamic cast to CalCalibMuSlope failed" << endreq;
      return StatusCode::FAILURE;
    }
  } else {
    msglog << MSG::INFO << "Enable to retrieve mu slopes from calib database" << endreq;            
  }

  // unpack crystal identification into tower, layer and column number
  int layer = xtalId.getLayer();
  int tower = xtalId.getTower();
  int col   = xtalId.getColumn();

  // create Volume Identifier for segment 0 of this crystal
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

  //get 3D transformation for segment 0 of this crystal
  detSvc->getTransform3DByID(segm0Id,&transf);

  //get position of the center of the segment 0
  Vector vect0 = transf.getTranslation();

			
  // create Volume Identifier for the last segment of this crystal
  idents::VolumeIdentifier segm11Id;

  // copy all fields from segm0Id, except segment number
  for(int ifield = 0; ifield<fSegment; ifield++)segm11Id.append(segm0Id[ifield]);
  segm11Id.append(m_nCsISeg-1); // set segment number for the last segment

  //get 3D transformation for the last segment of this crystal
  detSvc->getTransform3DByID(segm11Id,&transf);

  //get position of the center of the last segment
  Vector vect11 = transf.getTranslation();

  Point p0(0.,0.,0.);		

  // position of the crystal center
  Point pCenter = p0+(vect0+vect11)*0.5; 

  //normalized vector of the crystal direction 
  Vector dirXtal = 0.5*(vect11-vect0)*m_nCsISeg/(m_nCsISeg-1);	

	
  //extract energy for the best range for both crystal faces 	
  float calEnergyP,calEnergyN;
  if (m_pCalEnergyTool->calculate(xtalId,(idents::CalXtalId::AdcRange)rangeP,
                                  (idents::CalXtalId::POS),
                                  adcP,.5,calEnergyP)!=StatusCode::SUCCESS) {
    msglog << MSG::ERROR << "Error on calenergytool::calculate" << endreq;
  }
  if (m_pCalEnergyTool->calculate(xtalId,(idents::CalXtalId::AdcRange)rangeN,
                                  (idents::CalXtalId::NEG),
                                  adcN,.5,calEnergyN)!=StatusCode::SUCCESS) {
    msglog << MSG::ERROR << "Error on calenergytool::calculate" << endreq;
  }
  msglog << MSG::INFO << " eneP " << calEnergyP << " eneM " << calEnergyN << std::endl;
   
  double eneNeg = calEnergyN;
  double enePos = calEnergyP;

  // calculate slope - coefficient to convert light asymmetry into the 
  // deviation from crystal center
	
  // default value of slope - if calibration data base isn't accessible
  double slope = (1+m_lightAtt)/(1-m_lightAtt)*0.5;
    
  // extract mu slopes for selected ranges from calibration objects
 
  if(pMuSlopes){
    CalibData::RangeBase* pRangeP = pMuSlopes->getRange(xtalId, rangeP,idents::CalXtalId::POS);
    CalibData::RangeBase* pRangeM = pMuSlopes->getRange(xtalId, rangeN,idents::CalXtalId::NEG);
    
    CalibData::MuSlope* pMuSlopeP = dynamic_cast<CalibData::MuSlope * >(pRangeP);
    CalibData::MuSlope* pMuSlopeM = dynamic_cast<CalibData::MuSlope * >(pRangeM);

    double slopeP = pMuSlopeP->getSlope();
    double slopeM = pMuSlopeM->getSlope();
        
    slope = 0.5*(slopeP+slopeM)/m_CsILength;
  }
  double asym=0;

  //calculate light asymmetry
  if(enePos>0 && eneNeg>0)

    asym =  log(enePos/eneNeg);
  //            asym = (enePos-eneNeg)/(enePos+eneNeg);


  // longitudinal position in relative units (this value equal to +1 or -1
  // at the end of a crystal
  double relpos = asym*slope;

  // if the calculated position is outside the crystal - set it
  // to the nearest crystal end
  if(relpos > 1.0) relpos = 1.0;
  if(relpos <-1.0) relpos = -1.0;
        
  // calculate average position of energy deposit in this crystal
  Point pXtal = pCenter+dirXtal*relpos;
	

  // store calculated position in the reconstructed data
  // for the best readout range
  position = relpos;
  //(recData->getRangeRecData(0))->setPosition(pXtal);

  
  return StatusCode::SUCCESS;
}
