// File and Version Information:
//      $Header$
//
// Description:
//      AcdReconAlgV2 is a Gaudi algorithm which performs the ACD reconstruction.
//      Using the AcdDigi collection available on the TDS to compute a number
//      of reconstruction quantities.
//          
// Author(s):
//      Heather Kelly           

#include "AcdReconAlgV2.h"

#include "AcdPocaTool.h"

#include "AcdUtil/AcdTileDim.h"
#include "AcdUtil/AcdRibbonDim.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "GlastSvc/MonteCarlo/IMcGetEventInfoTool.h"

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"
#include "Event/Recon/CalRecon/CalClusterMap.h"
#include "Event/Recon/CalRecon/CalParams.h"
#include "Event/Recon/AcdRecon/AcdTkrIntersection.h"
#include "Event/Recon/AcdRecon/AcdTkrHitPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrGapPoca.h"
#include "Event/Recon/AcdRecon/AcdTkrPoca.h"
#include "LdfEvent/Gem.h"

#include "Event/MonteCarlo/McParticle.h"
#include "Event/TopLevel/MCEvent.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "AcdUtil/AcdDetectorList.h"
#include "AcdUtil/RayDoca.h"

#include "geometry/Ray.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "../AcdRecon/AcdReconFuncs.h"
#include "../AcdRecon/AcdReconFuncsV2.h"

#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <math.h>

// Define the fiducial volume of the LAT
// FIXME -- this should come for some xml reading service
//
// top is defined by planes at + 754.6 -> up to stacking of tiles
// sides are defined by planes at +-840.14
// the bottom of the ACD is at the z=-50 plane

// Later we add 10 cm to make sure that we catch everything
AcdRecon::AcdVolume AcdReconAlgV2::s_acdVolume;

double AcdReconAlgV2::s_vetoThresholdMeV;

unsigned int AcdReconAlgV2::s_numSideRows;

// Rogue value returned for DOCAs and Active Dist. calcs. when not Tile present
static double maxDoca = 2000.0;

//------------------------------------------------------------------------

//static const AlgFactory<AcdReconAlgV2>  Factory;
//const IAlgFactory& AcdReconAlgV2Factory = Factory;
DECLARE_ALGORITHM_FACTORY(AcdReconAlgV2);


// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor
AcdReconAlgV2::AcdReconAlgV2(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator),
  m_patRecMap( new AcdRecon::AcdPatRecMap ){
  declareProperty("intersectionToolName", m_intersectionToolName="AcdTkrIntersectToolV2");
  declareProperty("hitToolName",m_hitToolName="AcdPha2MipTool");
  declareProperty("pocaToolName",m_pocaToolName="AcdPocaToolV2");  
  declareProperty("propToolName",m_propToolName="G4PropagationTool");
  declareProperty("doBackSplash",m_doBackSplash=false);
  declareProperty("doVertexExtrap",m_doVertexExtrap=false);
  declareProperty("doDownwardExtrap",m_doDownwardExtrap=false);
  declareProperty("doCRTrackExtrap",m_doCRTrackExtrap=true);
  declareProperty("nPocaSave",m_nPocaSave=3);
  declareProperty("nHitPocaSave",m_nHitPocaSave=2);
  declareProperty("nTrigPocaSave",m_nTrigPocaSave=2);
  declareProperty("Tolerance",m_patRecTol=50.);
}



StatusCode AcdReconAlgV2::initialize ( ) {
    StatusCode sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
    
    m_glastDetSvc = 0;
    sc = service("GlastDetSvc", m_glastDetSvc, true);
    if (sc.isSuccess() ) {
        sc = m_glastDetSvc->queryInterface(IID_IGlastDetSvc, (void**)&m_glastDetSvc);
    }
    
    if( sc.isFailure() ) {
        MsgStream log(msgSvc(), name());
        log << MSG::ERROR << "AcdReconAlgV2 failed to get the GlastDetSvc" << endreq;
        return sc;
    }

    sc = service("AcdGeometrySvc", m_acdGeoSvc, true);
    if (sc.isSuccess() ) {
        sc = m_acdGeoSvc->queryInterface(IID_IAcdGeometrySvc, (void**)&m_acdGeoSvc);
    }
    if (sc.isFailure()) {
        MsgStream log(msgSvc(), name());
        log << MSG::ERROR << "AcdReconAlgV2 failed to get the AcdGeometerySvc" 
            << endreq;
        return sc;
    }
  
    m_geomMap = &m_acdGeoSvc->geomMap();
    m_geomMap->setAcdGeomSvc(*m_acdGeoSvc);


    if (m_intersectionToolName == "") 
        m_intersectionTool = 0;
    else {
        sc = toolSvc()->retrieveTool(m_intersectionToolName,  m_intersectionTool);
        if (sc.isFailure() ) {
            MsgStream log(msgSvc(), name());
            log << MSG::ERROR << "  Unable to create " << m_intersectionToolName << endreq;
            return sc;
        }
    }

    if (m_hitToolName == "") 
      m_hitTool = 0;
    else {
      sc = toolSvc()->retrieveTool(m_hitToolName,  m_hitTool);
      if (sc.isFailure() ) {
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR << "  Unable to create " << m_hitToolName << endreq;
    return sc;
      }
    }
    
    if (m_pocaToolName == "") 
      m_pocaTool = 0;
    else {
      sc = toolSvc()->retrieveTool(m_pocaToolName,  m_pocaTool);
      if (sc.isFailure() ) {
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR << "  Unable to create " << m_pocaToolName << endreq;
    return sc;
      }
    }

    if (m_propToolName == "") 
      m_G4PropTool = 0;
    else {
      sc = toolSvc()->retrieveTool(m_propToolName,  m_G4PropTool);
      if (sc.isFailure() ) {
    MsgStream log(msgSvc(), name());
    log << MSG::ERROR << "  Unable to create " << m_propToolName << endreq;
    return sc;
      }
    }

    // the pat rec stuff
    sc = fillPatRecMap();
    if (sc.isFailure()) {
      MsgStream log(msgSvc(), name());
      log << MSG::ERROR << "Failed to fill the ACD Pattern Recognition Map" << endreq;
      return sc;
    }
    
    // The track vector
    sc = toolSvc()->retrieveTool("TkrTrackVecTool", m_pTrackVec);
    if (sc.isFailure()) {
      MsgStream log(msgSvc(), name());
      log << MSG::ERROR << "Failed to get the TkrTrackVectTool" << endreq;
      return sc;
    }   

    getParameters();    
    return sc;
}


StatusCode AcdReconAlgV2::execute() {
    // Purpose and Method:  Called once per event.  This routine calls the functions
    //        that do the ACD reconstruction.
    // TDS Inputs:  EventModel::Digi::AcdDigiCol
    // Outputs:  Gaudi StatusCode
    // Dependencies:  The DOCA and active calculations rely upon the existance of the 
    //   reconstructed tracks on the TDS.  AcdReconAlgV2 will not fail if the tracks are
    //   not available, however, it is recommended that AcdReconAlgV2 be executed after
    //   some other routine can create the tracks on the TDS, if possible.
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    SmartDataPtr<Event::EventHeader> header(eventSvc(), "/Event");
    unsigned long evtId = (header) ? header->event() : 0;

    static bool firstEvent(true);
    if ( firstEvent ) {
      sc = m_acdGeoSvc->findCornerGaps();
      if ( sc.isFailure() ) {
    log << MSG::ERROR << "Failed to load ACD corner gap rays" << endreq;
      }
      firstEvent = false;
    }
    
    SmartDataPtr<Event::AcdDigiCol> acdDigiCol(eventSvc(), EventModel::Digi::AcdDigiCol);
    if (!acdDigiCol) {
        log << MSG::INFO << "No AcdDigiCol found on the TDS" << endreq;
        return sc;
    }

    // reset all member variables to their defaults
    clear();
        
    // run the reconstruction
    sc = reconstruct(acdDigiCol);
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "AcdReconAlgV2::reconstruct failed" << endreq;
      return sc;
    }
    return sc;
}


StatusCode AcdReconAlgV2::finalize() {    
    clear();
    return StatusCode::SUCCESS;
}


void AcdReconAlgV2::getParameters () {
    // Purpose and Method:  Retrieves constants using the GlastDetSvc.
    
    MsgStream   log( msgSvc(), name() );
    StatusCode sc;
    
    sc = m_glastDetSvc->getNumericConstByName("acd.vetoThreshold", &s_vetoThresholdMeV);
    if (sc.isFailure()) {
        log << MSG::INFO << "Unable to retrieve threshold, setting the value to 0.4 MeV" << endreq;
        s_vetoThresholdMeV = 0.4;
    }
    
    double temp;
    sc = m_glastDetSvc->getNumericConstByName("numSideRows", &temp);
    if (sc.isFailure()) {
        log << MSG::INFO << "Unable to retrieve numSideRows, setting the value to 4" << endreq;
        temp = 4.0;
    }
    s_numSideRows = (unsigned int) temp;        
}


void AcdReconAlgV2::clear() {
    // Purpose and Method:  Initializes all member variables

    // For new ActDist calc.
    m_hitMap.clear();

}

StatusCode AcdReconAlgV2::fillPatRecMap( ) {  
  AcdRecon::AcdPatRecMap_Rev revMap;
  // loop on the map of detector ID
  const std::map<idents::AcdId, int>& acdIDMap = m_acdGeoSvc->getAcdIdVolCountCol();
  for ( std::map<idents::AcdId, int>::const_iterator itr = acdIDMap.begin();
    itr != acdIDMap.end(); itr++ ) {
    const idents::AcdId& acdId = itr->first;
    if ( acdId.ribbon() ) {
      const AcdRibbonDim* ribbon = m_geomMap->getRibbon(acdId,*m_acdGeoSvc);
      if ( ribbon == 0 ) return StatusCode::FAILURE;
      if ( ! AcdRecon::addRibbonToPatRecMaps(*ribbon,*m_patRecMap,revMap,m_patRecTol) ) return false;
    } else if ( acdId.tile() ) {
      const AcdTileDim* tile = m_geomMap->getTile(acdId,*m_acdGeoSvc);
      if ( tile == 0 ) return StatusCode::FAILURE;
      if ( ! AcdRecon::addTileToPatRecMaps(*tile,*m_patRecMap,revMap,m_patRecTol) ) return false;
    }
  }
  unsigned int nBinUsed(0);
  unsigned int nAssoc(0);
  for ( AcdRecon::AcdPatRecMap::const_iterator itrPrint = m_patRecMap->begin();
    itrPrint != m_patRecMap->end(); itrPrint++ ) {
    unsigned int nElem = itrPrint->second.size();
    if ( nElem == 0 ) continue;
    nBinUsed++;
    nAssoc += nElem;    
  }
  MsgStream log( msgSvc(), name() );
  log << MSG::WARNING << "Built ACD ROI lookup table.  " << nBinUsed << " active bins and " << nAssoc << " entries." << endreq;
  return StatusCode::SUCCESS;
}


StatusCode AcdReconAlgV2::reconstruct (const Event::AcdDigiCol& digiCol) {
    // Purpose and Method:  Actually performs the ACD reconstruction.
    //        Counts the number of hit tiles and determines the total energy 
    //        deposited in the ACD.
    // Inputs:  digiCol is a pointer to the TDS ACD detector data.
    // Outputs:  Gaudi StatusCode:  returns StatusCode::FAILURE if an error 
    //           occurs
    // TDS Output:  EventModel::AcdRecon
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    // is this a periodic trigger?
    bool isPeriodicEvent(false);
    unsigned gemDeltaEventTime(0);
    
    SmartDataPtr<LdfEvent::Gem> gemTds(eventSvc(), "/Event/Gem");    
    if (gemTds) {
      gemDeltaEventTime = gemTds->deltaEventTime();
      isPeriodicEvent = gemTds->periodicSet();
    }

    // make the hits (with MIP peak data) and fill the hitMap
    static Event::AcdHitCol acdHits;
    if (m_hitTool != 0) {

      sc = m_hitTool->makeAcdHits(digiCol,isPeriodicEvent,gemDeltaEventTime,acdHits,m_hitMap);
      if ( sc.isFailure() ) {
    log << MSG::WARNING << "AcdHitTool Failed - we'll bravely carry on" 
            << endreq;
    sc = StatusCode::SUCCESS;
      }
    }
    
    Event::AcdPocaSet acdPocaSet;
    static Event::AcdTkrAssocCol acdTkrAssocs;
    static Event::AcdCalAssocCol acdCalAssocs;
    
    sc = trackDistances(acdHits,acdTkrAssocs,true);
    if (sc.isFailure()) {
      log << MSG::ERROR << "AcdReconAlgV2::trackDistances(up) failed" << endreq;
      return sc;
    }  

    if ( m_doDownwardExtrap ) {
      sc = trackDistances(acdHits,acdTkrAssocs,false);
      if (sc.isFailure()) {
	log << MSG::ERROR << "AcdReconAlgV2::trackDistances(down) failed" << endreq;
	return sc;
      }    
    }

    sc = calClusterDistances(acdHits,acdCalAssocs);
    if (sc.isFailure()) {
      log << MSG::ERROR << "AcdReconAlgV2::calClusterDistances failed" << endreq;
      return sc;
    }  

    if ( m_doVertexExtrap ) {
      
      sc = vertexDistances(acdHits,acdTkrAssocs,true);
      if (sc.isFailure()) {
	log << MSG::ERROR << "AcdReconAlgV2::vertexDistances(up) failed" << endreq;
	return sc;
      }
      if ( m_doDownwardExtrap ) {
	sc = vertexDistances(acdHits,acdTkrAssocs,false);
	if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdReconAlgV2::vertexDistances(down) failed" << endreq;
	  return sc;
	}    
      }
    }

    static Event::AcdPocaMap acdPocaMap;
    for (  Event::AcdPocaSet::iterator itrPoca = acdPocaSet.begin(); 
           itrPoca != acdPocaSet.end();
       itrPoca++ ) {
      Event::AcdTkrHitPoca* sortPoca = const_cast<Event::AcdTkrHitPoca*>(*itrPoca);
      if ( sortPoca == 0 ) continue;
      acdPocaMap.add(*sortPoca);
    }

    SmartDataPtr<Event::AcdReconV2> checkAcdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    Event::AcdReconV2* acdRecon(0);
    if (checkAcdRecTds) {
      log << MSG::DEBUG;
      if (log.isActive()) log.stream() << "AcdReconV2 data already on TDS!";
      log << endreq;
      acdRecon = checkAcdRecTds.ptr();
      acdRecon->clear();
    } else {
      acdRecon = new Event::AcdReconV2;
      sc = eventSvc()->registerObject(EventModel::AcdReconV2::Event, acdRecon);
      if (sc.isFailure()) {
    log << MSG::ERROR << "Failed to register AcdRecon" << endreq;
    return StatusCode::FAILURE;
      }
    }

    acdRecon->getHitCol().init(acdHits);
    acdRecon->getTkrAssocCol().init(acdTkrAssocs);
    acdRecon->getCalAssocCol().init(acdCalAssocs);
    
    sc = fillAcdEventTopology(acdHits,acdRecon->getEventTopology());
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "Failed to register AcdRecon" << endreq;
    }

    // Do the MC if needed
    sc = doMC(acdHits);
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "Failed to register AcdRecon" << endreq;
    }

    // ownership handed to TDS, clear local copies  
    acdTkrAssocs.clear();
    acdCalAssocs.clear();
    acdHits.clear();    

    if ( log.level() <= MSG::DEBUG ) {         
      log << MSG::DEBUG << "AcdReconAlgV2::reconstruct() finished" << std::endl << std::endl << endreq;
    }
    return sc;
}


StatusCode AcdReconAlgV2::doMC (const Event::AcdHitCol& acdHits) {
    // Purpose and Method:  Actually performs the ACD reconstruction.
    //        Counts the number of hit tiles and determines the total energy 
    //        deposited in the ACD.
    // Inputs:  digiCol is a pointer to the TDS ACD detector data.
    // Outputs:  Gaudi StatusCode:  returns StatusCode::FAILURE if an error 
    //           occurs
    // TDS Output:  EventModel::AcdRecon::AcdMCPoints
    //              EventModel::AcdRecon::AcdMCPocas
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    static Event::AcdTkrAssocCol acdTkrAssocs;
    sc = mcDistances(acdHits,acdTkrAssocs);
    if (sc.isFailure()) return sc;

    SmartDataPtr<Event::AcdTkrAssocCol> checkAcdMcTds(eventSvc(), EventModel::MC::McAcdTkrAssocCol);
    if (checkAcdMcTds) {
      log << MSG::DEBUG;
      if (log.isActive()) log.stream() << "AcdRecon MC already on TDS!";
      log << endreq;
      checkAcdMcTds->clear();
      checkAcdMcTds->init(acdTkrAssocs);
    } else {
      // create the TDS location for the AcdRecon
      Event::AcdTkrAssocCol *acdTkrAssocCol = new Event::AcdTkrAssocCol(acdTkrAssocs);
      sc = eventSvc()->registerObject(EventModel::MC::McAcdTkrAssocCol, acdTkrAssocCol);
      if (sc.isFailure()) {
    log << "Failed to register " << EventModel::MC::McAcdTkrAssocCol << endreq;
    return StatusCode::FAILURE;
      }
    }
    acdTkrAssocs.clear();
  
    return sc;
}

StatusCode AcdReconAlgV2::trackDistances(const Event::AcdHitCol& acdHits, 
					 Event::AcdTkrAssocCol& tkrAssocs, 
					 bool upward) {

    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities. 
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCole
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
     
    // Retrieve the information on fit tracks
    std::vector<Event::TkrTrack*> trackVec = m_pTrackVec->getTrackVec();
    int nTrack = trackVec.size();

    // Places to store the track endpoint and direction
    AcdRecon::TrackData extend;
    
    // where does this track leave the LAT?
    AcdRecon::ExitData exitData;

    for ( int iTrack(0); iTrack < nTrack; iTrack++) {

      const Event::TkrTrack* trackTds  = trackVec[iTrack];       // The TDS track
      bool isCR = trackTds->getStatusBits() & Event::TkrTrack::COSMICRAY;
      
      if ( isCR && (!m_doCRTrackExtrap) ) {
	// skip it
	continue;
      }

      // grap the track direction information      
      extend.m_index = isCR ? iTrack + 100 : iTrack;
      extend.m_upward = upward;
      const Event::TkrTrackHit* hitToUse(0);
      if (upward) {
	hitToUse = (*trackTds)[0];	  
      } else {
	hitToUse = (*trackTds)[trackTds->getNumHits() - 1];	  
      }

      if (hitToUse == 0){
	log << MSG::WARNING << "AcdRecon found a null hit on track " << extend.m_index << "  We'll bravely carry on." << endreq;
	continue;
      }
      extend.m_energy = hitToUse->getEnergy();
      AcdRecon::ReconFunctions::convertToAcdRep(hitToUse->getTrackParams(Event::TkrTrackHit::SMOOTHED),
						hitToUse->getZPlane(),extend);
      
      // get the LAT exit points
      if ( ! AcdRecon::ReconFunctions::exitsLat(extend,s_acdVolume,exitData) ) {
	log << MSG::WARNING << "AcdRecon::exitsLat() failed on " << (upward ? "Upward" : "Downward") 
	    << " end for track " << extend.m_index << endreq;
	return StatusCode::SUCCESS;
      }
      
      // keep track of all the pocas to hit tiles
      AcdRecon::PocaDataMap pocas;
      AcdRecon::PocaDataMap downwardPocas;
      
      // calculate all the distances to the hit tiles at once
      sc = hitDistances(extend,acdHits,exitData,pocas);
      if (sc.isFailure()) {
	log << MSG::ERROR << "AcdReconAlgV2::hitDistances (" << (upward ? "Upward" : "Downward") 
	    << ") failed for track " << extend.m_index << endreq;
	return sc;
      }
      
      // filter the lists for further procsessing
      AcdRecon::PocaDataPtrMap pocasCut;
      
      if ( m_pocaTool != 0 ) {
	sc = m_pocaTool->filter(pocas,pocasCut);
	if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdPocaTool::filter (" << (upward ? "Upward" : "Downward") 
	      << ") failed for track " << extend.m_index << endreq;
	  return sc;	  
	}
      }
      
      if ( log.level() <= MSG::DEBUG ) {
	log << MSG::DEBUG << "AcdReconAlgV2::trackDistances(" << iTrack << ") poca calculations finished." << std::endl << endreq;
      }
      
      // Now extrapolate the track as far as needed, 
      // this makes the AcdTkrHitPoca, AcdTkrGapPoca, AcdTkrPoint objects	
      std::vector<Event::AcdTkrHitPoca*> hitPocae;
      std::vector<Event::AcdTkrGapPoca*> gapPocae;
      int ssdVeto(0);
      float cornerDoca(0.);
      Event::AcdTkrPoint* trackPoint(0);
      
      // extrapolate the track upwards
      const Event::TkrTrackParams& trackPars = hitToUse->getTrackParams(Event::TkrTrackHit::SMOOTHED); 
      sc = extrapolateTrack(trackPars, extend, exitData, 
			    pocasCut, ssdVeto, hitPocae, gapPocae, trackPoint);
      if (sc.isFailure()) {
	log << MSG::ERROR << "AcdPocaTool::extrapolateTrack(" << (upward ? "Upward" : "Downward") 
	    << ") failed for track " << extend.m_index << endreq;
	return sc;
      }
      
      sc = calcCornerDoca(extend,cornerDoca);
      if (sc.isFailure()){
	log << MSG::ERROR << "AcdPocaTool::calcCornerDoca(" << (upward ? "Upward" : "Downward") 
	    << ") failed for track " << extend.m_index << endreq;
	return sc;
      }
      
      // sort the POCAe
      sortPocae(hitPocae,gapPocae);

      // Get the cone energies
      float energy15(0), energy30(0.), energy45(0.);
      float triggerEnergy15(0), triggerEnergy30(0.), triggerEnergy45(0.);
      
      sc = m_pocaTool->getConeDistances(hitPocae,
					energy15,energy30,energy45,
					triggerEnergy15,triggerEnergy30,triggerEnergy45);
      if (sc.isFailure()){
	log << MSG::ERROR << "AcdPocaTool::getConeDistances(" << (upward ? "Upward" : "Downward") 
	    << ") failed for track " << extend.m_index  << endreq;
	return sc;
      }

      // select only the best POCA
      sc = m_pocaTool->selectPocae(hitPocae,m_nPocaSave,m_nHitPocaSave,m_nTrigPocaSave);
      if (sc.isFailure()){
	log << MSG::ERROR << "AcdPocaTool::selectPocae(" << (upward ? "Upward" : "Downward") 
	    << ") failed for track"  << extend.m_index << endreq;
	return sc;
      }
      

      // Make the assoc object
      HepVector3D propVect = extend.m_current - extend.m_point;
      
      Event::AcdAssoc* theAssoc = 
	new Event::AcdAssoc(extend.m_index,upward,extend.m_energy,
			    extend.m_point,extend.m_dir,propVect.mag(),
			    extend.m_cov_orig,extend.m_cov_prop,
			    ssdVeto,cornerDoca,
			    energy15,energy30,energy45,
			    triggerEnergy15,triggerEnergy30,triggerEnergy45);
      
      sc = fillTkrAssoc(*theAssoc,hitPocae,gapPocae,trackPoint);
      if (sc.isFailure()){
	log << MSG::ERROR << "AcdPocaTool::fillTkrAssoc(" << (upward ? "Upward" : "Downward") 
	    << ") failed for track"  << extend.m_index << endreq;
	return sc;
      }
      tkrAssocs.push_back(theAssoc);
    }
   
    return sc;
    
}



StatusCode AcdReconAlgV2::vertexDistances(const Event::AcdHitCol& acdHits, 
					  Event::AcdTkrAssocCol& tkrAssocs,
					  bool upward) {

    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities. 
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCole
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
     
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrVertexCol> vertexTds(eventSvc(), EventModel::TkrRecon::TkrVertexCol);
    
    if (!vertexTds) {
      log << MSG::DEBUG << "No reconstructed vertex collection found on the TDS" 
      << endreq;
        return StatusCode::SUCCESS;
    } else {
      log << MSG::DEBUG << "AcdReconAlgV2::vertexDistances using " << vertexTds->size() << " vertices." << endreq;
    }
    
    int nVtx = vertexTds->size();
    if ( nVtx == 0 ) {      
      log << MSG::DEBUG << "No reconstructed vertices found on the TDS" 
      << endreq;
      return StatusCode::SUCCESS;
    }
    Event::TkrVertex* theVertex = (*vertexTds)[0];
    if ( theVertex == 0 ) {
      log << MSG::DEBUG << "Missed a vertex" 
      << endreq;
        return StatusCode::SUCCESS;
    }

    // Places to store the track endpoint and direction
    AcdRecon::TrackData extend;
    
    // where does this track leave the LAT?
    AcdRecon::ExitData exitData;

    // grap the vertex information
    extend.m_energy = theVertex->getEnergy();
    extend.m_index = -1;
    extend.m_upward = upward;
    AcdRecon::ReconFunctions::convertToAcdRep(theVertex->getVertexParams(),
					      theVertex->getPosition().z(),extend);

    // get the LAT exit points
    if ( ! AcdRecon::ReconFunctions::exitsLat(extend,s_acdVolume,exitData) ) {
      log << MSG::WARNING << "AcdRecon::exitsLat() failed on " << (upward ? "Upward" : "Downward") << " end." << endreq;
      return StatusCode::SUCCESS;
    }
    

    // keep track of all the pocas to hit tiles
    AcdRecon::PocaDataMap pocas;
    
    // calculate all the distances to the hit tiles at once
    sc = hitDistances(extend,acdHits,exitData,pocas);
    if (sc.isFailure()) return sc;

    // filter the lists for further procsessing
    AcdRecon::PocaDataPtrMap pocasCut;
    
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->filter(pocas,pocasCut);
      if (sc.isFailure()) return sc;
    }

    if ( log.level() <= MSG::DEBUG ) {
      log << MSG::DEBUG << "AcdReconAlgV2::vertexDistances() poca calculations finished." << std::endl << endreq;
    }

    std::vector<Event::AcdTkrHitPoca*> hitPocae;
    std::vector<Event::AcdTkrGapPoca*> gapPocae;
    Event::AcdTkrPoint* tkrPoint(0);
    float cornerDoca(0.);
    int ssdVeto(0);


    // extrapolate the track upwards
    sc = extrapolateVertex(extend, exitData, pocasCut, 
			   ssdVeto, hitPocae, gapPocae, tkrPoint);
    if (sc.isFailure()) {
      log << MSG::ERROR << "AcdPocaTool::extrapolateVertex(" << (upward ? "Upward" : "Downward") << ") failed." << endreq;
      return sc;
    }

    sc = calcCornerDoca(extend,cornerDoca);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::calcCornerDoca(" << (upward ? "Upward" : "Downward") << ") failed." << endreq;
      return sc;
    }

    
    // sort the POCAe
    sortPocae(hitPocae,gapPocae);
    
    // Get the cone energies
    float energy15(0), energy30(0.), energy45(0.);
    float triggerEnergy15(0), triggerEnergy30(0.), triggerEnergy45(0.);
    
    sc = m_pocaTool->getConeDistances(hitPocae,
				      energy15,energy30,energy45,
				      triggerEnergy15,triggerEnergy30,triggerEnergy45);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::getConeDistances(" << (upward ? "Upward" : "Downward") << ") failed" << endreq;
      return sc;
    }
    
    // select only the best POCA
    sc = m_pocaTool->selectPocae(hitPocae,m_nPocaSave,m_nHitPocaSave,m_nTrigPocaSave);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::getConeDistances(" << (upward ? "Upward" : "Downward") << ") failed" << endreq;
      return sc;
    }
    
    HepVector3D propVect = extend.m_current - extend.m_point;    
    Event::AcdAssoc* theAssoc = 
      new Event::AcdAssoc(-1,upward,extend.m_energy,
			  extend.m_point,extend.m_dir,propVect.mag(),
			  extend.m_cov_orig,extend.m_cov_prop,
			  ssdVeto,cornerDoca,
			  energy15,energy30,energy45,
			  triggerEnergy15,triggerEnergy30,triggerEnergy45);
    sc = fillTkrAssoc(*theAssoc,hitPocae,gapPocae,tkrPoint);

    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::fillTkrAssoc(" << (upward ? "Upward" : "Downward") << ") failed." << endreq;
      return sc;
    }

    tkrAssocs.push_back(theAssoc);

    return sc;
    
};


StatusCode AcdReconAlgV2::mcDistances(const Event::AcdHitCol& acdHits, 
                    Event::AcdTkrAssocCol& tkrAssocs) {

    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities.  Updates the
    // local data members m_doca, m_rowDocaCol, m_act_dist, m_rowActDistCol
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCole
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
     
    // Retrieve the information on mc tracks
    SmartDataPtr<Event::McParticleCol> pMcParticle(eventSvc(), EventModel::MC::McParticleCol);
    if (pMcParticle == 0) {
      return StatusCode::SUCCESS;
    } else {
      log << MSG::DEBUG << "AcdReconAlgV2::mcDistances()" << endreq;
    }
    Event::McParticleCol::const_iterator mcFirst = pMcParticle->begin();
    Event::McParticle* mcPart = *mcFirst;
    if ( mcPart == 0 ) {
      return StatusCode::SUCCESS;
    }
    
    // Places to store the track endpoint and direction
    AcdRecon::TrackData extend;

    // where does this track enter the LAT?
    AcdRecon::ExitData enter;

    // grap the vertex information
    extend.m_energy = mcPart->initialFourMomentum().e();
    extend.m_index = -2;
    extend.m_upward = (extend.m_dir.z() > 0);
    extend.m_point = extend.m_current = mcPart->initialPosition();
    extend.m_dir   = mcPart->initialFourMomentum().vect().unit();

    // get the LAT exit points
    if ( ! AcdRecon::ReconFunctions::entersLat(extend,s_acdVolume,enter) ) {
      log << MSG::DEBUG << "AcdRecon::entersLat() failed on MC track - we'll bravely carry on" << endreq;
      return StatusCode::SUCCESS;
    }

    // keep track of all the pocas to hit tiles
    AcdRecon::PocaDataMap pocas;
    
    // calculate all the distances to the hit tiles at once
    sc = hitDistances(extend,acdHits,enter,pocas);
    if (sc.isFailure()) return sc;

    // filter the lists for further procsessing
    AcdRecon::PocaDataPtrMap pocasCut;
    
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->filter(pocas,pocasCut);
      if (sc.isFailure()) return sc;
    }

    std::vector<Event::AcdTkrHitPoca*> hitPocae;
    std::vector<Event::AcdTkrGapPoca*> gapPocae;
    Event::AcdTkrPoint* point(0);
    float cornerDoca(0.);
    int ssdVeto(0);  

    // extrapolate the track upwards
    sc = extrapolateVertex(extend, enter, pocasCut, 
               ssdVeto,hitPocae,gapPocae,point);
    if (sc.isFailure()) {
      log << MSG::ERROR << "AcdPocaTool::extrapolateVertex(mc) failed" << endreq;
      return sc;
    }
    sc = calcCornerDoca(extend,cornerDoca);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::calcCornerDoca(mc) failed" << endreq;
      return sc;
    }

    
    // sort the POCAe
    sortPocae(hitPocae,gapPocae);
    
    // Get the cone energies
    float energy15(0), energy30(0.), energy45(0.);
    float triggerEnergy15(0), triggerEnergy30(0.), triggerEnergy45(0.);
    
    sc = m_pocaTool->getConeDistances(hitPocae,
				      energy15,energy30,energy45,
				      triggerEnergy15,triggerEnergy30,triggerEnergy45);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::getConeDistances(" << (extend.m_upward ? "Upward" : "Downward") << ") failed" << endreq;
      return sc;
    }
    
    // select only the best POCA
    sc = m_pocaTool->selectPocae(hitPocae,m_nPocaSave,m_nHitPocaSave,m_nTrigPocaSave);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::getConeDistances(" << (extend.m_upward ? "Upward" : "Downward") << ") failed" << endreq;
      return sc;
    }
      

    HepVector3D propVect = extend.m_current - extend.m_point;    
    Event::AcdAssoc* assoc = 
      new Event::AcdAssoc(-2,true,extend.m_energy,
			  extend.m_point,extend.m_dir,propVect.mag(),
			  extend.m_cov_orig,extend.m_cov_prop,
			  ssdVeto,cornerDoca,
			  energy15,energy30,energy45,
			  triggerEnergy15,triggerEnergy30,triggerEnergy45);
    sc = fillTkrAssoc(*assoc,hitPocae,gapPocae,point);
    if (sc.isFailure()){
      log << MSG::ERROR << "AcdPocaTool::fillTkrAssoc(up) failed" << endreq;
      return sc;
    }
    tkrAssocs.push_back(assoc);

    log << MSG::DEBUG << "AcdReconAlgV2::mcDistances() finished" << endreq;

    return sc;  
}


StatusCode AcdReconAlgV2::hitDistances(const AcdRecon::TrackData& aTrack, 
                       const Event::AcdHitCol& acdHits, 
                       const AcdRecon::ExitData& data,
                       AcdRecon::PocaDataMap& pocaMap) {
  /// get the all the distances to hit tiles for track in one direction

  StatusCode sc = StatusCode::SUCCESS;

  MsgStream   log( msgSvc(), name() );
  log << MSG::DEBUG << "AcdReconAlgV2::hitDistances for track " << aTrack.m_index 
      << " going " << (aTrack.m_upward ? "up" : "down" ) <<  " with " << acdHits.size() << " hits." << endreq;
  
  // Get the list of relevent detector elements
  std::set<idents::AcdId> idSet;
  if ( data.m_face != 5 ) {
    Vector tDir(aTrack.m_dir.x(),aTrack.m_dir.y(),aTrack.m_dir.z());
    if ( ! AcdRecon::buildElementSet( data.m_x, tDir, data.m_arcTol, idSet, *m_patRecMap ) ) {
      log << MSG::ERROR << "Failed to build element set" << endreq;
      return StatusCode::FAILURE;
    }
  }

  for ( std::set<idents::AcdId>::const_iterator itr = idSet.begin(); itr != idSet.end(); itr++ ) {
    log << MSG::DEBUG << "ACD Element near point " << data.m_x.x() << ',' <<  data.m_x.y() << ',' <<  data.m_x.z() << " :" 
    << itr->id() << endreq;
      // get the data object to store all the computations
    AcdRecon::PocaData& pocaData = pocaMap[*itr];
    sc = elemDistances(aTrack,*itr,pocaData);
    if ( sc.isFailure() ) {
      // Already gave error message, just return
      return sc;
    }
    if ( pocaData.m_arcLengthPlane < 0. ) {
      log << MSG::DEBUG << "ACD Element near point " << data.m_x.x() << ',' <<  data.m_x.y() << ',' <<  data.m_x.z() << " :" 
      << itr->id() << " Has plane intersect at " << pocaData.m_arcLengthPlane << endreq;
      pocaMap.erase(*itr);
    }
  }

  for (Event::AcdHitCol::const_iterator acdHitIt = acdHits.begin(); acdHitIt != acdHits.end(); acdHitIt++) {
    idents::AcdId acdId = (*acdHitIt)->getAcdId();
    // check to see if we have it already
    AcdRecon::PocaDataMap::iterator itr = pocaMap.find(acdId);
    // it if is already there just tag the fact that it has a hit
    if ( itr != pocaMap.end() ) { 
      itr->second.m_hasHit = true;
      continue;
    }
    // not there, make it
    AcdRecon::PocaData pocaData;
    sc = elemDistances(aTrack,acdId,pocaData);
    if ( sc.isFailure() ) {
      // Already gave error message, just return
      return sc;
    } 
    if ( pocaData.m_arcLengthPlane >= 0. ) {
        pocaMap[acdId] = pocaData;
    }
    pocaData.m_hasHit = true;

  }
  return sc;
}


StatusCode AcdReconAlgV2::elemDistances(const AcdRecon::TrackData& aTrack, 
                      const idents::AcdId& acdId,
                      AcdRecon::PocaData& pocaData) {
    StatusCode sc = StatusCode::SUCCESS;

    /// get the distances for a single element and a single track in one detector
    if (acdId.na()) { 
      MsgStream log( msgSvc(), name() );
      log << MSG::ERROR << "Skipping NA hit" << endreq;
      return StatusCode::SUCCESS;
    } else if (acdId.tile()) {      
      // Tiles
      const AcdTileDim* tileDim = m_geomMap->getTile(acdId,*m_acdGeoSvc);
      sc = tileDim->statusCode();
      if ( sc.isFailure() ) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() 
        << endreq;
    return sc;
      }
      if ( m_pocaTool != 0 ) {
    sc = m_pocaTool->tileDistances(*tileDim,aTrack,pocaData);
      } else {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "No Poca Tool" << endreq;
    return StatusCode::FAILURE;
      }
      if ( sc.isFailure() ) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to get hit distances for a tile" 
        << acdId.id() << endreq;
    return sc;
      }
    } else if ( acdId.ribbon() ) {
      // Ribbons
      const AcdRibbonDim* ribbonDim = m_geomMap->getRibbon(acdId,*m_acdGeoSvc);
      sc = ribbonDim->statusCode();
      if ( sc.isFailure() ) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to get geom for a ribbon " << acdId.id() 
        << endreq;
    return sc;
      }
      if ( m_pocaTool != 0 ) {
    sc = m_pocaTool->ribbonDistances(*ribbonDim,aTrack,pocaData);
      } else {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "No Poca Tool" << endreq;  
    return StatusCode::FAILURE;
      }
      if ( sc.isFailure() ) {
    MsgStream log( msgSvc(), name() );
    log << MSG::ERROR << "Failed to get hit distances for a ribbon" 
        << acdId.id() << endreq;
    return sc;
      }
    } else {
      MsgStream log( msgSvc(), name() );
      log << MSG::ERROR << "Neither NA, nor tile, nor ribbon" << endreq;
      return StatusCode::FAILURE;
    }
    return StatusCode::SUCCESS;    
}


StatusCode AcdReconAlgV2::calClusterDistances(const Event::AcdHitCol& acdHits, 
                          Event::AcdCalAssocCol& calAssocs) {

    // Purpose and Method:  Calculates the DOCA and Active Distance quantities for
    // the CAL cluster axes. Similar to the calculation done for tracks and vertices.
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
     
    // Retrieve the information on the CAL cluster(s)
    SmartDataPtr<Event::CalClusterMap> calClusterMapTds(eventSvc(), 
                             EventModel::CalRecon::CalClusterMap);
    
    if (!calClusterMapTds) {
      log << MSG::DEBUG << "No reconstructed cluster(s) found on the TDS" << std::endl
      << endreq;
      return StatusCode::SUCCESS;
    } else {
      log << MSG::DEBUG << "AcdReconAlgV2::calClusterDistances using " << calClusterMapTds->getRawClusterVec().size() << " clusters." << endreq;
    }   


    // Places to store the track endpoint and direction
    AcdRecon::TrackData upwardExtend;

    // where does this track leave the LAT volume?
    AcdRecon::ExitData upwardExit;

    int iCluster(-1);

    // Converting from the old scheme of looping over all clusters in the CalClusterCol (which included ubers)
    // to the new scheme of getting clusters from the CalClusterMap. If we want to still include the ubers here,
    // this is not obviously the most efficient way to do this, but it does bring the code up to date in the 
    // event that CalClusterCol goes away in the future. 

    // Start with outer loop over the elements of the map
    for(Event::CalClusterMap::iterator calMapItr  = calClusterMapTds->begin();
                                       calMapItr != calClusterMapTds->end();
                                       calMapItr++)
    {
        Event::CalClusterVec calClusterVec = calMapItr->second;

        // Inner loop is now over clusters
        for(Event::CalClusterVec::iterator calVecItr  = calClusterVec.begin();
                                           calVecItr != calClusterVec.end();
                                           calVecItr++)
        {
            Event::CalCluster* calClusterTds = *calVecItr;
            iCluster++;

            // Check to be sure the moments analysis has run (and there is a valid axis)
            if (!calClusterTds->checkStatusBit(Event::CalCluster::MOMENTS)) continue;

            Event::CalParams clusterParams = calClusterTds->getMomParams();

            // grab the cluster direction information
            upwardExtend.m_energy = calClusterTds->getMomParams().getEnergy();
            upwardExtend.m_index = iCluster;
            upwardExtend.m_upward = true;
            upwardExtend.m_tracker = false;

            // Ph.Bruel: use the cal position corrected for the hodoscopic effect (using the cal direction)
            Vector cal_dir  = clusterParams.getAxis();
            Point cal_pos  = clusterParams.getCentroid();
            if(cal_dir.z()!=0) cal_pos = clusterParams.getCorCentroid(cal_dir);

            upwardExtend.m_dir.set(cal_dir.x(), 
                       cal_dir.y(), 
                       cal_dir.z());

            upwardExtend.m_point.set(cal_pos.x(), 
                       cal_pos.y(), 
                       cal_pos.z());

            // check if cluster enters the acd 
            bool entersLat = AcdRecon::ReconFunctions::entersLat(upwardExtend,s_acdVolume,upwardExit);

            // overwrittes values of upwardExtend
            AcdRecon::ReconFunctions::convertToAcdRep(clusterParams, upwardExtend);

            // get the LAT exit points 
            if ( ! AcdRecon::ReconFunctions::exitsLat(upwardExtend,s_acdVolume,upwardExit) ) {
                if ( !entersLat ) {
                    // Don't care if exitLat fails if cluster never entered the ACD volume
                    log << MSG::DEBUG << "AcdRecon::exitsLat() CalCluster doesn't enter LAT - we'll bravely carry on" << endreq;
                    return StatusCode::SUCCESS;
                } else {
                    // Something else happened...
                    log << MSG::WARNING << "AcdRecon::exitsLat() failed for CalCluster - we'll bravely carry on " 
			<< "v = (" << cal_dir.x() << "," << cal_dir.y() << "," << cal_dir.z() << ") " 		      
			<< "x = (" << cal_pos.x() << "," << cal_pos.y() << "," << cal_pos.z() << ")" 
			<< endreq;
                    return StatusCode::SUCCESS;
                }
            }
            
            // keep track of all the pocas to hit tiles
            AcdRecon::PocaDataMap upwardPocas;

            // calculate all the distances to the hit tiles at once
            sc = hitDistances(upwardExtend,acdHits,upwardExit,upwardPocas);
            if (sc.isFailure()) {
	        log << MSG::ERROR << "AcdReconAlgV2::hitDistances(up) failed for cluster " << upwardExtend.m_index << endreq;
                return sc;
            }
        
            // filter the lists for further procsessing
            AcdRecon::PocaDataPtrMap upPocasCut;
        
            if ( m_pocaTool != 0 ) {
                sc = m_pocaTool->filter(upwardPocas,upPocasCut);
                if (sc.isFailure()) {
                    log << MSG::ERROR << "AcdPocaTool::filter(up) failed for cluster " << upwardExtend.m_index << endreq;
                    return sc;    
                }
            }
       
            if ( log.level() <= MSG::DEBUG ) {
              log << MSG::DEBUG << "AcdReconAlgV2::calClusterDistances(" << iCluster << ") poca calculations finished." << std::endl << endreq;
            }
            
            // Now extrapolate the track as far as needed, 
            // this makes the AcdTkrHitPoca, AcdTkrGapPoca, AcdTkrPoint objects 
            std::vector<Event::AcdTkrHitPoca*> upHitPocae;
            std::vector<Event::AcdTkrGapPoca*> upGapPocae;
            int ssdVetoUp(0);
            float cornerDocaUp(0.);
            Event::AcdTkrPoint* upPoint(0);

            // extrapolate the track upwards
            sc = m_intersectionTool->makeTkrPoint(upwardExtend,upwardExit,upPoint);
            if ( sc.isFailure() ){
              log << MSG::ERROR << "AcdTkrIntersectionTool::makeTkrPoint failed for cluster " << upwardExtend.m_index << endreq;
              return sc;
            }

            if ( upPoint != 0 ) {
              if ( log.level() <= MSG::DEBUG ) {    
                upPoint->writeOut(log);
              }
            }

            static Event::TkrTrackParams trackParams;
            AcdRecon::ReconFunctions::convertToTkrRep(upwardExtend, trackParams);
            sc = extrapolateTrack(trackParams, upwardExtend, upwardExit, 
                          upPocasCut, ssdVetoUp, upHitPocae, upGapPocae, upPoint);
            if (sc.isFailure()) {
              log << MSG::ERROR << "AcdPocaTool::extrapolateTrack(up) failed for cluster " << upwardExtend.m_index << endreq;
              return sc;
            }

            sc = calcCornerDoca(upwardExtend,cornerDocaUp);
            if (sc.isFailure()){
              log << MSG::ERROR << "AcdPocaTool::calcCornerDoca(down) failed for cluster " << upwardExtend.m_index << endreq;
              return sc;
            }

            // sort the POCAe
            sortPocae(upHitPocae,upGapPocae);

	    // Get the cone energies
	    float energy15(0), energy30(0.), energy45(0.);
	    float triggerEnergy15(0), triggerEnergy30(0.), triggerEnergy45(0.);
	    
	    sc = m_pocaTool->getConeDistances(upHitPocae,
					      energy15,energy30,energy45,
					      triggerEnergy15,triggerEnergy30,triggerEnergy45);
	    if (sc.isFailure()){
	      log << MSG::ERROR << "AcdPocaTool::getConeDistances( ) failed for cluster " << upwardExtend.m_index << endreq;
	      return sc;
	    }
	    
	    // select only the best POCA
	    sc = m_pocaTool->selectPocae(upHitPocae,m_nPocaSave,m_nHitPocaSave,m_nTrigPocaSave);
	    if (sc.isFailure()){
	      log << MSG::ERROR << "AcdPocaTool::selectPocae( ) failed for cluster " << upwardExtend.m_index << endreq;
	      return sc;
	    }
 	    
            HepVector3D propVect = upwardExtend.m_current - upwardExtend.m_point;

            Event::AcdAssoc* upAssoc = 
              new Event::AcdAssoc(upwardExtend.m_index,true,upwardExtend.m_energy,
				  upwardExtend.m_point,upwardExtend.m_dir,propVect.mag(),
				  upwardExtend.m_cov_orig,upwardExtend.m_cov_prop,
				  ssdVetoUp,cornerDocaUp,
				  energy15,energy30,energy45,
				  triggerEnergy15,triggerEnergy30,triggerEnergy45);

            sc = fillCalAssoc(*upAssoc,upHitPocae,upGapPocae,upPoint);
            if (sc.isFailure()){
              log << MSG::ERROR << "AcdPocaTool::fillCalAssoc(up) failed for cluster " << upwardExtend.m_index << endreq;
              return sc;
            }
            calAssocs.push_back(upAssoc);
        }
    }

    return sc;   
}


StatusCode AcdReconAlgV2::extrapolateTrack(const Event::TkrTrackParams& trackParams,
                     const AcdRecon::TrackData& trackData,
                     const AcdRecon::ExitData& isectData,
                     AcdRecon::PocaDataPtrMap& pocaDataMap,
                     int& ssdVeto,
                     std::vector<Event::AcdTkrHitPoca*>& hitPocae,
                     std::vector<Event::AcdTkrGapPoca*>& gapPocae,
                     Event::AcdTkrPoint*& point){

  MsgStream   log( msgSvc(), name() );
  StatusCode sc = StatusCode::SUCCESS;  

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlgV2::extrapolateTrack(" 
    << trackData.m_index << ',' << (trackData.m_upward ? "up" : "down" ) << ')' << endreq;
  }

  ssdVeto = 0;

  // This is for any new Poca Data we need to make.  
  // They go out of scope as soon as we have converted them to TDS format
  std::list<AcdRecon::PocaData> ownedPocaData;
  
  // build all the intersections
  if ( m_intersectionTool != 0 ) {
    try {
      sc = m_intersectionTool->makeIntersections(*m_G4PropTool,trackParams,trackData,isectData,m_hitMap,
                         *m_geomMap,pocaDataMap,ownedPocaData,gapPocae);
    } catch (...) {
      log << MSG::ERROR << "Caught exception using propagator to make intersection with track " << trackData.m_index 
    << ".  No ACD intersections calculated for that track." <<   endreq;
      // Don't crash, just continue.
      return sc;
    }
  }
  if ( sc.isFailure() ) {
    log << MSG::ERROR << "AcdTkrIntersectionTool::makeIntersections failed" << endreq;
    return sc;
  }
  
  // build all the pocas
  Event::TkrTrackParams next_params;
  AcdRecon::PocaDataPtrMap::const_iterator itr;
  for ( itr = pocaDataMap.begin(); itr != pocaDataMap.end(); itr++ ) {
    
    AcdRecon::AcdHitMap::const_iterator findHit = m_hitMap.find(itr->first);
    const Event::AcdHit* theHit = findHit != m_hitMap.end() ? findHit->second : 0;

    AcdRecon::PocaData& pocaData = *(itr->second);
    Event::AcdTkrHitPoca* aPoca(0);
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->makePoca(trackData,pocaData,theHit,aPoca);
    }
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "AcdPocaTool::makePoca failed" << endreq;
      return sc;
    }
    
    if ( aPoca != 0 ) {      
      if ( log.level() <= MSG::DEBUG ) {    
    aPoca->writeOut(log);
      }
      hitPocae.push_back(aPoca);
    }
  }

  if ( log.level() <= MSG::DEBUG ) {        
    for ( std::vector<Event::AcdTkrGapPoca*>::const_iterator itrG = gapPocae.begin();
      itrG != gapPocae.end(); itrG++ ) {
      (*itrG)->writeOut(log);
    }
  }

  // build the TrkPoint  
  AcdRecon::ReconFunctions::propagateToArcLength(*m_G4PropTool,isectData.m_arcLength,trackData,next_params);
  Event::AcdTkrPoint* exitPoint(0);
  if ( m_intersectionTool != 0 ) {
    sc = m_intersectionTool->makeTkrPoint(trackData,isectData,exitPoint);
    if ( sc.isFailure() ){
      log << MSG::ERROR << "AcdTkrIntersectionTool::makeTkrPoint failed" << endreq;
      return sc;
    }
    if ( exitPoint != 0 ) {
      if ( log.level() <= MSG::DEBUG ) {    
    exitPoint->writeOut(log);
      }
      point = exitPoint;
    }
  }

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlgV2::extrapolateTrack() finished" << std::endl << endreq;
  }

  return sc;

}

StatusCode AcdReconAlgV2::extrapolateVertex(const AcdRecon::TrackData& trackData,
                      const AcdRecon::ExitData& isectData,
                      AcdRecon::PocaDataPtrMap& pocaDataMap,
                      int& ssdVeto,
                      std::vector<Event::AcdTkrHitPoca*>& hitPocae,
                      std::vector<Event::AcdTkrGapPoca*>& /* gapPocae */,
                      Event::AcdTkrPoint*& point) {

  MsgStream   log( msgSvc(), name() );
  StatusCode sc = StatusCode::SUCCESS;  
 
  ssdVeto = 0;

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlgV2::extrapolateVertex(" 
    << (trackData.m_upward ? "up" : "down" ) << ')' << endreq;
  }


  // first figure out how far to extrapolate track
  double maxArcLength(0.);

  // figure out which direction we are going in
  bool forward = isectData.m_arcLength > 0;
  
  // run the propagator out to the right arclength
  if ( !forward ) maxArcLength *= -1.;

  // build all the pocas
  for ( AcdRecon::PocaDataPtrMap::const_iterator itr = pocaDataMap.begin(); itr != pocaDataMap.end(); itr++ ) {
    AcdRecon::AcdHitMap::const_iterator findHit = m_hitMap.find(itr->first);
    const Event::AcdHit* theHit = findHit != m_hitMap.end() ? findHit->second : 0;
    
    AcdRecon::PocaData& pocaData = *(itr->second);
    //float pocaArcLength = forward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    //Event::TkrTrackParams next_params = m_G4PropTool->getTrackParams(pocaArcLength,startEnergy,true);
    //AcdRecon::projectErrorAtPoca(trackData,next_params,pocaData.m_poca,pocaData.m_pocaVector,pocaData.m_active3DErr);
    Event::AcdTkrHitPoca* aPoca(0);
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->makePoca(trackData,pocaData,theHit,aPoca);
    }
    if ( sc.isFailure() ) return sc;
    if ( aPoca != 0 ) {
      aPoca->writeOut(log);
      hitPocae.push_back(aPoca);
    }
  }

  // build the TrkPoint
  Event::AcdTkrPoint* exitPoint(0);
  if ( m_intersectionTool != 0 ) {
    sc = m_intersectionTool->makeTkrPoint(trackData,isectData,exitPoint);
    if ( sc.isFailure() ) return sc;
    if ( exitPoint != 0 ) {
      exitPoint->writeOut(log);
      point = exitPoint;    
    }
  }

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlgV2::extrapolateTrack() finished" << std::endl << endreq;
  }

  return sc;
  
}

///  Pick the best pocae for a track
StatusCode AcdReconAlgV2::sortPocae(std::vector<Event::AcdTkrHitPoca*>& hitPocae, 
                  std::vector<Event::AcdTkrGapPoca*>& gapPocae ) {

  std::sort( hitPocae.begin(), hitPocae.end(), hit_pointer_less() );
  std::sort( gapPocae.begin(), gapPocae.end(), gap_pointer_less() );
  return StatusCode::SUCCESS;
}


StatusCode AcdReconAlgV2::calcCornerDoca(const AcdRecon::TrackData& track, 
					 float &return_dist) {

    return_dist = maxDoca;

    // Use this flag to determine whether to apply sign or not.. if we never
    // find a corner DOCA where the intersection is within the limits of the 
    // side, the return distance should remain at maxDoca
    bool foundOne = false;
    double testArcLen(0.), testRayLen(0.), testDist(0.);
    Point testPoint;
    Vector testDir;
    int testRegion;

    // iterate over all corner gaps
    unsigned int iCorner;
    for (iCorner=0; iCorner<4; iCorner++) {
        const Ray& gapRay = m_acdGeoSvc->getCornerGapRay(iCorner);

        // Compute DOCA between the track and gap ray   
    AcdRecon::rayDoca_withCorner(track,gapRay,testArcLen,testRayLen,testDist,testPoint,testDir,testRegion);

    // only take forward intersections
    if ( testArcLen < 0. ) continue;
    if ( testDist > return_dist ) continue;
    foundOne = true;
    return_dist = (return_dist > testDist) ? testDist : return_dist;
    }

    // If no DOCAs were found, just return, and skip the sign calculation
    if(!foundOne) return StatusCode::SUCCESS;

    // Now we have DOCA to the corners
    // Next compute sign based on (Tkr1X0*Tkr1YDir - Tkr1Y0*Tkr1XDir)
    
    //std::cout << return_dist << std::endl;
    float sign = ( (track.m_point.x() * track.m_dir.y()) - (track.m_point.y() * track.m_dir.x()) ) > 0 ? 1. : -1;
    return_dist *= sign;    
    return StatusCode::SUCCESS;
  
}


/// Fill an AcdTkrAssoc with data
StatusCode AcdReconAlgV2::fillTkrAssoc(Event::AcdAssoc& assoc,
                     const std::vector<Event::AcdTkrHitPoca*>& hitPocae,
                     const std::vector<Event::AcdTkrGapPoca*>& gapPocae,
                     Event::AcdTkrPoint* point) {

  for ( std::vector<Event::AcdTkrHitPoca*>::const_iterator itrHitPoca = hitPocae.begin(); 
    itrHitPoca != hitPocae.end(); itrHitPoca++ ) {
    assoc.addHitPoca(*(*itrHitPoca));
  }
  
  for ( std::vector<Event::AcdTkrGapPoca*>::const_iterator itrGapPoca = gapPocae.begin(); 
    itrGapPoca != gapPocae.end(); itrGapPoca++ ) {
    assoc.addGapPoca(*(*itrGapPoca));
  }

  if ( point != 0 ) {
    assoc.setPoint(*point);
  }
  return StatusCode::SUCCESS;
}

/// Fill an AcdCalAssoc with data
StatusCode AcdReconAlgV2::fillCalAssoc(Event::AcdAssoc& assoc,
                     const std::vector<Event::AcdTkrHitPoca*>& hitPocae,
                     const std::vector<Event::AcdTkrGapPoca*>& gapPocae,
                     Event::AcdTkrPoint* point) {

  for ( std::vector<Event::AcdTkrHitPoca*>::const_iterator itrHitPoca = hitPocae.begin(); 
    itrHitPoca != hitPocae.end(); itrHitPoca++ ) {
    assoc.addHitPoca(*(*itrHitPoca));
  }
  
  for ( std::vector<Event::AcdTkrGapPoca*>::const_iterator itrGapPoca = gapPocae.begin(); 
    itrGapPoca != gapPocae.end(); itrGapPoca++ ) {
    assoc.addGapPoca(*(*itrGapPoca));
  }

  if ( point != 0 ) {
    assoc.setPoint(*point);
  }
  return StatusCode::SUCCESS;
}


/// Fill the AcdEventTopology object
StatusCode AcdReconAlgV2::fillAcdEventTopology(const Event::AcdHitCol& acdHits,
                         Event::AcdEventTopology& evtTopo) {
  
  unsigned tileCount(0),ribbonCount(0),vetoCount(0),tileVeto(0);
  float totalTileEnergy(0),totalRibbonEnergy(0);  
  float tileEnergy(0),ribbonEnergy(0);  
  float ghostTileEnergy(0),ghostRibbonEnergy(0);  
  float triggerTileEnergy(0),triggerRibbonEnergy(0);  
  unsigned nTilesTop(0);  
  unsigned nTilesSideRow[4] = {0,0,0,0};  
  unsigned nTilesSideFace[4] = {0,0,0,0};  
  unsigned nVetoTop(0);  
  unsigned nVetoSideRow[4] = {0,0,0,0};  
  unsigned nVetoSideFace[4] = {0,0,0,0};  
  float tileEnergyTop(0);  
  float tileEnergySideRow[4] = {0.,0.,0.,0.};    
  float tileEnergySideFace[4] = {0.,0.,0.,0.};  
  unsigned nSidesHit(0),nSidesVeto(0);

  std::set<int> sidesHit;
  std::set<int> sidesVeto;

  for ( Event::AcdHitCol::const_iterator itr = acdHits.begin(); itr != acdHits.end(); itr++ ) {
    Event::AcdHit* theHit = *itr;
    const idents::AcdId& id = theHit->getAcdId();
    bool hasGhost   = theHit->getGhost();
    bool hasTrigger = theHit->getTriggerVeto();
    if ( hasTrigger)   vetoCount++;
    if ( id.tile() ) {
      tileCount++;
      sidesHit.insert( id.face() );
      if ( hasTrigger ) {
    sidesVeto.insert( id.face() );
    tileVeto++;
    switch ( id.face() ) {
    case 0:
      nVetoTop++;
      break;
    case 1: case 2: case 3: case 4:
      nVetoSideFace[id.face()-1]++;
      nVetoSideRow[id.row()]++;
    }
      }
      float energy = theHit->tileEnergy();
      totalTileEnergy += energy;
      tileEnergy += energy*!hasGhost;
      ghostTileEnergy += energy*hasGhost;
      triggerTileEnergy += energy*hasTrigger;
      switch ( id.face() ) {
      case 0:
    nTilesTop++;
    tileEnergyTop += energy;
    break;
      case 1:
      case 2:
      case 3:
      case 4:
    nTilesSideFace[id.face() - 1]++;
    tileEnergySideFace[id.face() - 1] += energy;
    nTilesSideRow[id.row()]++;
    tileEnergySideRow[id.row()] += energy;
    break;
      }      
    } else if ( id.ribbon() ) {
      ribbonCount++;
      float energy = theHit->ribbonEnergy( Event::AcdHit::A ) + theHit->ribbonEnergy( Event::AcdHit::B );
      energy /= 2.;
      totalRibbonEnergy += energy;
      ribbonEnergy += energy*!hasGhost;
      ghostRibbonEnergy += energy*hasGhost;
      triggerRibbonEnergy += energy*hasTrigger;
    }
  }

  nSidesHit = sidesHit.size();
  nSidesVeto = sidesVeto.size();
  evtTopo.set( tileCount,  ribbonCount,  vetoCount, tileVeto,
           totalTileEnergy, totalRibbonEnergy, tileEnergy, ribbonEnergy,
               ghostTileEnergy, ghostRibbonEnergy, triggerTileEnergy, triggerRibbonEnergy,
           nTilesTop,  nTilesSideRow,  nTilesSideFace,
           nVetoTop,  nVetoSideRow,  nVetoSideFace,
           tileEnergyTop,  tileEnergySideRow,  tileEnergySideFace,
           nSidesHit,  nSidesVeto);

  return StatusCode::SUCCESS;

}


