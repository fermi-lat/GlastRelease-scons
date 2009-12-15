// File and Version Information:
//      $Header$
//
// Description:
//      AcdReconAlg is a Gaudi algorithm which performs the ACD reconstruction.
//      Using the AcdDigi collection available on the TDS to compute a number
//      of reconstruction quantities.
//          
// Author(s):
//      Heather Kelly			

#include "AcdReconAlg.h"

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
#include "../AcdRecon/AcdTkrParams.h"


#include <algorithm>
#include <cstdio>
#include <stdlib.h>


// Define the fiducial volume of the LAT
// FIXME -- this should come for some xml reading service
//
// top is defined by planes at + 754.6 -> up to stacking of tiles
// sides are defined by planes at +-840.14
// the bottom of the ACD is at the z=-50 plane

// Later we add 10 cm to make sure that we catch everything
AcdRecon::AcdVolume AcdReconAlg::s_acdVolume;

double AcdReconAlg::s_vetoThresholdMeV;

unsigned int AcdReconAlg::s_numSideRows;

// Rogue value returned for DOCAs and Active Dist. calcs. when not Tile present
static double maxDoca = 2000.0;

//------------------------------------------------------------------------

static const AlgFactory<AcdReconAlg>  Factory;
const IAlgFactory& AcdReconAlgFactory = Factory;

// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor
AcdReconAlg::AcdReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
  declareProperty("intersectionToolName", m_intersectionToolName="AcdTkrIntersectTool");
  declareProperty("hitToolName",m_hitToolName="AcdPha2MipTool");
  declareProperty("pocaToolName",m_pocaToolName="AcdPocaTool");  
  declareProperty("propToolName",m_propToolName="G4PropagationTool");
  declareProperty("doBackSplash",m_doBackSplash=false);
}



StatusCode AcdReconAlg::initialize ( ) {
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
        log << MSG::ERROR << "AcdReconAlg failed to get the GlastDetSvc" << endreq;
        return sc;
    }

    sc = service("AcdGeometrySvc", m_acdGeoSvc, true);
    if (sc.isSuccess() ) {
        sc = m_acdGeoSvc->queryInterface(IID_IAcdGeometrySvc, (void**)&m_acdGeoSvc);
    }
    if (sc.isFailure()) {
        MsgStream log(msgSvc(), name());
        log << MSG::ERROR << "AcdReconAlg failed to get the AcdGeometerySvc" 
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
	
    getParameters(); 	
    return sc;
}


StatusCode AcdReconAlg::execute() {
    // Purpose and Method:  Called once per event.  This routine calls the functions
    //        that do the ACD reconstruction.
    // TDS Inputs:  EventModel::Digi::AcdDigiCol
    // Outputs:  Gaudi StatusCode
    // Dependencies:  The DOCA and active calculations rely upon the existance of the 
    //   reconstructed tracks on the TDS.  AcdReconAlg will not fail if the tracks are
    //   not available, however, it is recommended that AcdReconAlg be executed after
    //   some other routine can create the tracks on the TDS, if possible.
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

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
      log << MSG::ERROR << "AcdReconAlg::reconstruct failed" << endreq;
      return sc;
    }
    return sc;
}


StatusCode AcdReconAlg::finalize() {    
    clear();
    return StatusCode::SUCCESS;
}


void AcdReconAlg::getParameters () {
    // Purpose and Method:  Retrieves constans using the GlastDetSvc.
	
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


void AcdReconAlg::clear() {
    // Purpose and Method:  Initializes all member variables
    m_totEnergy = 0.0;
    m_tileCount = 0;
    m_totRibbonEnergy = 0.0;
    m_ribbonCount = 0;
    m_gammaDoca = maxDoca;
    m_doca = maxDoca;
    m_rowDocaCol.clear();
    // one for each side, plus one for the top
    m_rowDocaCol.resize(s_numSideRows+1, maxDoca);  
    m_rowActDistCol.clear();
    // one for each side, plus one for the top
    m_rowActDistCol.resize(s_numSideRows+1, -maxDoca);

    // For new ActDist calc.
    m_rowActDist3DCol.clear();
    m_rowActDist3DCol.resize(s_numSideRows+1, -maxDoca);
    m_rowActDist3DCol_down.clear();
    m_rowActDist3DCol_down.resize(s_numSideRows+1, -maxDoca);

    m_energyCol.clear();
    m_idCol.clear();
    m_energyRibbonCol.clear();
    m_idRibbonCol.clear();
    m_act_dist = -maxDoca;
    m_act_dist3D = -maxDoca;  // active distance 3D upward tracks
    m_act_dist3D_down = -maxDoca;  // active distance downward tracks
    m_ribbon_act_dist = -maxDoca;

    m_cornerDoca = maxDoca;

    m_hitMap.clear();

    idents::AcdId resetId;
    resetId.na(1);
    m_minDocaId = resetId;
    m_ribbon_act_dist_id = resetId;
    m_maxActDistId = resetId;
    m_maxActDist3DId = resetId;
    m_maxActDist3DId_down = resetId;

    // Don't reset the map, Geometry should remain the same over the life of the job
    //m_geomMap.reset();
}


StatusCode AcdReconAlg::reconstruct (const Event::AcdDigiCol& digiCol) {
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
	
    Event::AcdDigiCol::const_iterator acdDigiIt;    
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        
        idents::AcdId id = (*acdDigiIt)->getId();

        if (id.tile()) {
	  m_tileCount++;
	  double tileEnergy = (*acdDigiIt)->getEnergy();
	  m_totEnergy += tileEnergy;
	  
	  // Temporarily populate reconstructed energy collection with digi energy
	  m_idCol.push_back(id);
	  m_energyCol.push_back(tileEnergy);
        } else { // otherwise this is a ribbon
	  m_ribbonCount++;
	  m_totRibbonEnergy += (*acdDigiIt)->getEnergy();
	  m_idRibbonCol.push_back(id);
	  m_energyRibbonCol.push_back((*acdDigiIt)->getEnergy());
        }
    }
	    
    if ( log.level() <= MSG::DEBUG ) {         
      log << MSG::DEBUG << "AcdReconAlg::reconstruct()." << std::endl
	  << "\tnTiles = " << m_tileCount 
	  << ", total energy = " << m_totEnergy
	  << ", numRibbons = " << m_ribbonCount
	  << ", total Ribbon energy = " << m_totRibbonEnergy << std::endl << endreq;
    }

    Event::AcdPocaSet acdPocaSet;
    static Event::AcdTkrIntersectionCol acdIntersections;
    static Event::AcdTkrGapPocaCol acdGapPocas;
    static Event::AcdTkrPocaCol acdPocas;
    static Event::AcdTkrHitPocaCol acdHitPocas;
    static Event::AcdTkrPointCol acdPoints;
    static Event::AcdSplashVarsCol acdSplashVars;

    sc = trackDistances(acdHits,acdPocaSet,acdIntersections,acdGapPocas,acdPoints);
    if (sc.isFailure()) {
      log << MSG::ERROR << "AcdReconAlg::trackDistances failed" << endreq;
      return sc;
    }    

    sc = vertexDistances(acdHits,acdPocaSet,acdPoints);
    if (sc.isFailure()) {
      log << MSG::ERROR << "AcdReconAlg::vertexDistances failed" << endreq;
      return sc;
    }

    static Event::AcdPocaMap acdPocaMap;
    for (  Event::AcdPocaSet::iterator itrPoca = acdPocaSet.begin(); 
           itrPoca != acdPocaSet.end();
	   itrPoca++ ) {
      Event::AcdTkrHitPoca* sortPoca = const_cast<Event::AcdTkrHitPoca*>(*itrPoca);
      if ( sortPoca == 0 ) continue;
      acdPocaMap.add(*sortPoca);
      acdHitPocas.push_back(sortPoca);
    }

    if ( m_doBackSplash ) {
      sc = doBacksplash(digiCol,acdSplashVars);
      if (sc.isFailure()) {
	log << MSG::ERROR << "AcdReconAlg::doBacksplash failed" << endreq;
	return sc;    
      }
    }

    if ( log.level() <= MSG::DEBUG ) {      
      log << MSG::DEBUG << "AcdActiveDistance: " << m_act_dist << std::endl << endreq;
    }

    SmartDataPtr<Event::AcdRecon> checkAcdRecTds(eventSvc(), 
                 EventModel::AcdRecon::Event);  
    if (checkAcdRecTds) {
        log << MSG::DEBUG;
        if (log.isActive()) log.stream() << "AcdRecon data already on TDS!";
        log << endreq;
        checkAcdRecTds->clear();
        checkAcdRecTds->init(m_totEnergy, m_totRibbonEnergy, 
			     m_tileCount, m_ribbonCount, 
			     m_gammaDoca, m_doca, m_minDocaId,  
			     m_act_dist, m_maxActDistId, 
			     m_rowDocaCol, m_rowActDistCol, m_idCol, 
                             m_energyCol,
			     m_act_dist3D, m_maxActDist3DId, m_rowActDist3DCol,
			     acdIntersections, acdPocas, acdHits, acdHitPocas, 
                             acdGapPocas, acdPoints, acdSplashVars, m_ribbon_act_dist, 
                             m_ribbon_act_dist_id, m_cornerDoca);

        checkAcdRecTds->initActDist3D_down(m_act_dist3D_down, 
                                           m_maxActDist3DId_down, 
                                           m_rowActDist3DCol_down);
    } else {
        // create the TDS location for the AcdRecon
        Event::AcdRecon *acdRecon = new Event::AcdRecon(m_totEnergy, 
                                    m_totRibbonEnergy, m_tileCount, 
                                    m_ribbonCount, m_gammaDoca, m_doca, 
                                    m_minDocaId, m_act_dist, m_maxActDistId, 
                                    m_rowDocaCol, m_rowActDistCol, m_idCol, 
                                    m_energyCol, m_ribbon_act_dist, 
                                    m_ribbon_act_dist_id, acdIntersections, 
                                    acdPocas, acdHits, acdHitPocas, 
                                    acdGapPocas, acdPoints, acdSplashVars, m_act_dist3D, 
                                    m_maxActDist3DId, m_rowActDist3DCol, 
                                    m_cornerDoca);

        acdRecon->initActDist3D_down(m_act_dist3D_down, m_maxActDist3DId_down, 
                                     m_rowActDist3DCol_down);

        sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, acdRecon);
        if (sc.isFailure()) {
            log << MSG::ERROR << "Failed to register AcdRecon" << endreq;
            return StatusCode::FAILURE;
        }
    }

    // ownership handed to TDS, clear local copies	
    acdIntersections.clear();
    acdPocas.clear();
    acdPocaMap.clear();
    acdGapPocas.clear();
    acdHitPocas.clear();
    acdPoints.clear();
    acdSplashVars.clear();

    // Do the MC if needed
    sc = doMC(acdHits);
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "Failed to register AcdRecon" << endreq;
    }
    acdHits.clear();	

    if ( log.level() <= MSG::DEBUG ) {         
      log << MSG::DEBUG << "AcdReconAlg::reconstruct() finished" << std::endl << std::endl << endreq;
    }
    return sc;
}


StatusCode AcdReconAlg::doMC (const Event::AcdHitCol& acdHits) {
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
		
    Event::AcdPocaSet acdPocaSet;
    static Event::AcdTkrHitPocaCol acdHitPocas;
    static Event::AcdTkrPointCol acdPoints;
    
    sc = mcDistances(acdHits,acdPocaSet,acdPoints);
    if (sc.isFailure()) return sc;

    static Event::AcdPocaMap acdPocaMap;
    for (  Event::AcdPocaSet::iterator itrPoca = acdPocaSet.begin(); 
           itrPoca != acdPocaSet.end();
	   itrPoca++ ) {
      Event::AcdTkrHitPoca* sortPoca = const_cast<Event::AcdTkrHitPoca*>(*itrPoca);
      if ( sortPoca == 0 ) continue;
      acdPocaMap.add(*sortPoca);
      acdHitPocas.push_back(sortPoca);
    }

    SmartDataPtr<Event::AcdTkrPointCol> checkAcdMcPointTds(eventSvc(), EventModel::MC::McAcdTkrPointCol);
    if (checkAcdMcPointTds) {
      log << MSG::DEBUG;
      if (log.isActive()) log.stream() << "AcdRecon data already on TDS!";
      log << endreq;
      checkAcdMcPointTds->clear();
      checkAcdMcPointTds->init(acdPoints);

    } else {
      // create the TDS location for the AcdRecon
      Event::AcdTkrPointCol *acdPointCol = new Event::AcdTkrPointCol(acdPoints);
      sc = eventSvc()->registerObject(EventModel::MC::McAcdTkrPointCol, acdPointCol);
      if (sc.isFailure()) {
	log << "Failed to register " << EventModel::MC::McAcdTkrPointCol << endreq;
	return StatusCode::FAILURE;
      }
    }

    SmartDataPtr<Event::AcdTkrHitPocaCol> checkAcdMcHitPocaTds(eventSvc(), EventModel::MC::McAcdTkrHitPocaCol);
    if (checkAcdMcHitPocaTds) {
      log << MSG::DEBUG;
      if (log.isActive()) log.stream() << "AcdRecon data already on TDS!";
      log << endreq;
      checkAcdMcHitPocaTds->clear();
      checkAcdMcHitPocaTds->init(acdHitPocas);

    } else {
      // create the TDS location for the AcdRecon
      Event::AcdTkrHitPocaCol *acdHitPocaCol = new Event::AcdTkrHitPocaCol(acdHitPocas);
      sc = eventSvc()->registerObject(EventModel::MC::McAcdTkrHitPocaCol, acdHitPocaCol);
      if (sc.isFailure()) {
	log << "Failed to register" << EventModel::MC::McAcdTkrHitPocaCol << endreq;
	return StatusCode::FAILURE;
      }
    }

    // ownership handed to TDS, clear local copies	
    acdPocaSet.clear();
    acdPocaMap.clear();
    acdHitPocas.clear();
    acdPoints.clear();
    return sc;
}

StatusCode AcdReconAlg::trackDistances(const Event::AcdHitCol& acdHits, 
				 Event::AcdPocaSet& pocaSet,
				 Event::AcdTkrIntersectionCol& acdIntersections,
				 Event::AcdTkrGapPocaCol& gapPocas,
				 Event::AcdTkrPointCol& exitPoints) {

    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities.  Updates the
    // local data members m_doca, m_rowDocaCol, m_act_dist, m_rowActDistCol
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCole
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
     
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrTrackCol> tracksTds(eventSvc(), 
                                          EventModel::TkrRecon::TkrTrackCol);
	
    if (!tracksTds) {
      log << MSG::DEBUG << "No reconstructed tracks found on the TDS" << std::endl
            << endreq;
        return StatusCode::SUCCESS;
    } else {
      log << MSG::DEBUG << "AcdReconAlg::trackDistances using " << tracksTds->size() << " tracks." << endreq;
    }	

    bool firstPassDone = false;

    // Places to store the track endpoint and direction
    AcdRecon::TrackData upwardExtend;
    AcdRecon::TrackData downwardExtend;
    
    // where does this track leave the LAT?
    AcdRecon::ExitData upwardExit;
    AcdRecon::ExitData downwardExit;

    int iTrack(-1);
    Event::TkrTrackColPtr trkPtr = tracksTds->begin();
    while(trkPtr != tracksTds->end())
    {

        const Event::TkrTrack* trackTds  = *trkPtr++;       // The TDS track
	iTrack++;

	// grap the track direction information
	const Event::TkrTrackHit* firstHit = (*trackTds)[0];
	upwardExtend.m_point = firstHit->getPoint(Event::TkrTrackHit::SMOOTHED);
	upwardExtend.m_dir   = -(firstHit->getDirection(Event::TkrTrackHit::SMOOTHED));
	upwardExtend.m_energy = firstHit->getEnergy();
	upwardExtend.m_index = iTrack;
	upwardExtend.m_upward = true;

	const unsigned int lastHitIdx = trackTds->getNumHits() - 1;
	const Event::TkrTrackHit* lastHit = (*trackTds)[lastHitIdx];
	downwardExtend.m_point = lastHit->getPoint(Event::TkrTrackHit::SMOOTHED);
	downwardExtend.m_dir   = lastHit->getDirection(Event::TkrTrackHit::SMOOTHED);
	downwardExtend.m_energy = lastHit->getEnergy();
	downwardExtend.m_index = iTrack;
	downwardExtend.m_upward = false;

	// get the LAT exit points
	if ( ! AcdRecon::exitsLat(upwardExtend,s_acdVolume,upwardExit) ) {
	  log << MSG::WARNING << "AcdRecon::exitsLat() failed on upward end - we'll bravely carry on" << endreq;
	  return StatusCode::SUCCESS;
	} else {
	  if ( log.level() <= MSG::DEBUG ) {
	    log << MSG::DEBUG << "  Acd Track exit(up).  ";
	    writeExitPoint(log.stream(),upwardExtend,upwardExit);
	    log << endreq;
	  }
	}
	
	if ( ! AcdRecon::exitsLat(downwardExtend,s_acdVolume,downwardExit) ) {
	  log << MSG::WARNING << "AcdRecon::exitsLat() failed on downward end - we'll bravely carry on" << endreq;
	  return StatusCode::SUCCESS;
	} else {
	  if ( log.level() <= MSG::DEBUG ) {
	    log << MSG::DEBUG << "  Acd Track exit(down).  ";
	    writeExitPoint(log.stream(),downwardExtend,downwardExit);
	    log << std::endl << endreq;
	  }
	}	  
	
	// keep track of all the pocas to hit tiles
	AcdRecon::PocaDataMap upwardPocas;
	AcdRecon::PocaDataMap downwardPocas;

	// calculate all the distances to the hit tiles at once
	sc = hitDistances(upwardExtend,acdHits,upwardPocas);
	if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdReconAlg::hitDistances(up) failed" << endreq;
	  return sc;
	}

	sc = hitDistances(downwardExtend,acdHits,downwardPocas);
	if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdReconAlg::hitDistances(down) failed" << endreq;
	  return sc;
	}

        // grab the best New Active Distance Calc (3D) values
        sc = tileActDist(upwardPocas, m_rowActDist3DCol, m_act_dist3D, m_maxActDist3DId);

        if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdReconAlg::tileActDist(up) failed" << endreq;
	  return sc;
	}

        // Now calculate using downward Pocas
        sc = tileActDist(downwardPocas, m_rowActDist3DCol_down,
                         m_act_dist3D_down, m_maxActDist3DId_down);
        if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdReconAlg::tileActDist(up) failed" << endreq;
	  return sc;
	}

	// grab the best "Active Distance" from ribbons
        sc = hitRibbonDist(upwardPocas, m_ribbon_act_dist, m_ribbon_act_dist_id);
        if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdReconAlg::hitRibbonDist(up) failed" << endreq;
	  return sc;
	}
	
	// filter the lists for further procsessing
	AcdRecon::PocaDataPtrMap upPocasCut;
	AcdRecon::PocaDataPtrMap downPocasCut;
	
	if ( m_pocaTool != 0 ) {
	  sc = m_pocaTool->filter(upwardPocas,upPocasCut);
	  if (sc.isFailure()) {
	    log << MSG::ERROR << "AcdPocaTool::filter(up) failed" << endreq;
	    return sc;	  
	  }
	  sc = m_pocaTool->filter(downwardPocas,downPocasCut);
	  if (sc.isFailure()) { 
	    log << MSG::ERROR << "AcdPocaTool::filter(down) failed" << endreq;
	    return sc;
	  }
	}
	
	if ( log.level() <= MSG::DEBUG ) {
	  log << MSG::DEBUG << "AcdReconAlg::trackDistances(" << iTrack << ") poca calculations finished." << std::endl << endreq;
	}
	
	// Now extrapolate the track as far as needed, 
	// this makes the AcdTkrPoca and AcdTkrIntersection objects
	Event::AcdTkrIntersectionCol upwardIntersections;
	Event::AcdTkrIntersectionCol downwardIntersections;

	Event::AcdTkrGapPocaCol upGapPocas;
	Event::AcdTkrGapPocaCol downGapPocas;

	// extrapolate the track upwards
	sc = extrapolateTrack(*trackTds, upwardExtend, upPocasCut, upwardExit, 
			      pocaSet, upwardIntersections, upGapPocas, exitPoints);
	if (sc.isFailure()) {
	  log << MSG::ERROR << "AcdPocaTool::extrapolateTrack(up) failed" << endreq;
	  return sc;
	}
	for ( Event::AcdTkrIntersectionCol::iterator itrU = upwardIntersections.begin(); 
	      itrU != upwardIntersections.end(); ++itrU ) {
	  acdIntersections.push_back(*itrU);
	}
	for ( std::vector<Event::AcdTkrGapPoca*>::iterator itrGU = upGapPocas.begin();
	      itrGU != upGapPocas.end(); ++itrGU ) {
	  gapPocas.push_back(*itrGU);
	}

	// extrapolate the track downwards
	sc = extrapolateTrack(*trackTds, downwardExtend, downPocasCut, 
                              downwardExit, pocaSet, downwardIntersections, 
                              downGapPocas, exitPoints);
	if (sc.isFailure()){
	  log << MSG::ERROR << "AcdPocaTool::extrapolateTrack(down) failed" << endreq;
	  return sc;
	}
	for ( Event::AcdTkrIntersectionCol::iterator itrD = downwardIntersections.begin(); 
	      itrD != downwardIntersections.end(); ++itrD ) {
	  acdIntersections.push_back(*itrD);
	}
	for ( std::vector<Event::AcdTkrGapPoca*>::iterator itrGD = downGapPocas.begin();
	      itrGD != downGapPocas.end(); ++itrGD ) {
	  gapPocas.push_back(*itrGD);
	}

	// clean up
	upwardIntersections.clear();
	upGapPocas.clear();	
	downwardIntersections.clear();
	downGapPocas.clear();
	
        if ((!firstPassDone) && (m_calcCornerDoca)) {
            // take the first track, since it is the best track ignore the rest
            calcCornerDoca(trackTds->getInitialPosition(), 
                          -(trackTds->getInitialDirection()), m_cornerDoca);
            // First track in the list, is the reconstructed gamma
            firstPassDone = true;
	    //std::cout << "CornerDoca: " << m_cornerDoca << std::endl;
        }
    }
    
    if ( log.level() <= MSG::DEBUG ) {
      log << MSG::DEBUG << "AcdCornerDoca = " << m_cornerDoca << std::endl << endreq;
    }
    return sc;
    
}



StatusCode AcdReconAlg::vertexDistances(const Event::AcdHitCol& acdHits, 
					Event::AcdPocaSet& pocaSet,
					Event::AcdTkrPointCol& exitPoints) {

    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities.  Updates the
    // local data members m_doca, m_rowDocaCol, m_act_dist, m_rowActDistCol
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
      log << MSG::DEBUG << "AcdReconAlg::vertexDistances using " << vertexTds->size() << " vertices." << endreq;
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
    AcdRecon::TrackData upwardExtend;
    AcdRecon::TrackData downwardExtend;
    
    // where does this track leave the LAT?
    AcdRecon::ExitData upwardExit;
    AcdRecon::ExitData downwardExit;

    // grap the vertex information
    upwardExtend.m_point = theVertex->getPosition();
    upwardExtend.m_dir   = -(theVertex->getDirection());
    upwardExtend.m_energy = theVertex->getEnergy();
    upwardExtend.m_index = -1;
    upwardExtend.m_upward = true;
    Point upPoint(upwardExtend.m_point.x(),upwardExtend.m_point.y(),upwardExtend.m_point.z());
    Vector upDir(upwardExtend.m_dir.x(),upwardExtend.m_dir.y(),upwardExtend.m_dir.z());   

    downwardExtend.m_point = theVertex->getPosition();
    downwardExtend.m_dir   = theVertex->getDirection();
    downwardExtend.m_energy = theVertex->getEnergy();
    downwardExtend.m_index = -1;
    downwardExtend.m_upward = false;
    Point downPoint(downwardExtend.m_point.x(),downwardExtend.m_point.y(),downwardExtend.m_point.z());
    Vector downDir(downwardExtend.m_dir.x(),downwardExtend.m_dir.y(),downwardExtend.m_dir.z());

    // get the LAT exit points
    if ( ! AcdRecon::exitsLat(upwardExtend,s_acdVolume,upwardExit) ) {
      log << MSG::WARNING << "AcdRecon::exitsLat() failed on upward end - we'll bravely carry on" << endreq;
      return StatusCode::SUCCESS;
    } else {
      if ( log.level() <= MSG::DEBUG ) {
	log << MSG::DEBUG << "  Acd Vertex exit(up).  ";
	writeExitPoint(log.stream(),upwardExtend,upwardExit);      
	log << endreq;
      }
    }
    
    if ( ! AcdRecon::exitsLat(downwardExtend,s_acdVolume,downwardExit) ) {
      log << MSG::WARNING << "AcdRecon::exitsLat() failed on downward end - we'll bravely carry on" << endreq;
      return StatusCode::SUCCESS;
    } else {
      if ( log.level() <= MSG::DEBUG ) {
	log << MSG::DEBUG << "  Acd Vertex exit(down).  ";
	writeExitPoint(log.stream(),downwardExtend,downwardExit);	
	log << endreq;
      }
    }

    // keep track of all the pocas to hit tiles
    AcdRecon::PocaDataMap upwardPocas;
    AcdRecon::PocaDataMap downwardPocas;
    
    // calculate all the distances to the hit tiles at once
    sc = hitDistances(upwardExtend,acdHits,upwardPocas);
    if (sc.isFailure()) return sc;

    sc = hitDistances(downwardExtend,acdHits,downwardPocas);
    if (sc.isFailure()) return sc;

    // filter the lists for further procsessing
    AcdRecon::PocaDataPtrMap upPocasCut;
    AcdRecon::PocaDataPtrMap downPocasCut;
    
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->filter(upwardPocas,upPocasCut);
      if (sc.isFailure()) return sc;
      sc = m_pocaTool->filter(downwardPocas,downPocasCut);
      if (sc.isFailure()) return sc;
    }

    if ( log.level() <= MSG::DEBUG ) {
      log << MSG::DEBUG << "AcdReconAlg::vertexDistances() poca calculations finished." << std::endl << endreq;
    }

    // extrapolate the track upwards
    sc = extrapolateVertex(upwardExtend, upPocasCut, upwardExit, pocaSet, exitPoints);
    if (sc.isFailure()) return sc;

    // extrapolate the track downwards
    sc = extrapolateVertex(downwardExtend, downPocasCut, downwardExit, pocaSet, exitPoints);
    if (sc.isFailure()) return sc;

    if ( log.level() <= MSG::DEBUG ) {
      log << MSG::DEBUG << "AcdReconAlg::vertexDistances() finished" << std::endl << endreq;
    }

    return sc;
	
};


StatusCode AcdReconAlg::mcDistances(const Event::AcdHitCol& acdHits, 
				    Event::AcdPocaSet& pocaSet,
				    Event::AcdTkrPointCol& exitPoints) {

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
      log << MSG::DEBUG << "AcdReconAlg::mcDistances()" << endreq;
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
    extend.m_point = mcPart->initialPosition();
    extend.m_dir   = mcPart->initialFourMomentum().vect().unit();
    extend.m_energy = mcPart->initialFourMomentum().e();
    extend.m_index = -2;
    extend.m_upward = (extend.m_dir.z() > 0);
    Point point(extend.m_point.x(),extend.m_point.y(),extend.m_point.z());
    Vector dir(extend.m_dir.x(),extend.m_dir.y(),extend.m_dir.z());   

    // get the LAT exit points
    if ( ! AcdRecon::entersLat(extend,s_acdVolume,enter) ) {
      log << MSG::INFO << "AcdRecon::entersLat() failed on MC track- we'll bravely carry on" << endreq;
      return StatusCode::SUCCESS;
    } else {
      if ( log.level() <= MSG::DEBUG ) {
	log << MSG::DEBUG << "  Acd MC entry().  ";
	writeExitPoint(log.stream(),extend,enter);
	log << endreq;
      }
    }

    // keep track of all the pocas to hit tiles
    AcdRecon::PocaDataMap pocas;
    
    // calculate all the distances to the hit tiles at once
    sc = hitDistances(extend,acdHits,pocas);
    if (sc.isFailure()) return sc;

    // filter the lists for further procsessing
    AcdRecon::PocaDataPtrMap pocasCut;
    
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->filter(pocas,pocasCut);
      if (sc.isFailure()) return sc;
    }

    // extrapolate the track upwards
    sc = extrapolateVertex(extend, pocasCut, enter, pocaSet, exitPoints);
    if (sc.isFailure()) return sc;

    log << MSG::DEBUG << "AcdReconAlg::mcDistances() finished" << endreq;

    return sc;	
}

StatusCode AcdReconAlg::hitDistances(const AcdRecon::TrackData& aTrack, 
                                     const Event::AcdHitCol& acdHits, 
				     AcdRecon::PocaDataMap& pocaMap) {
  /// get the all the distances to hit tiles for track in one direction

  StatusCode sc = StatusCode::SUCCESS;

  MsgStream   log( msgSvc(), name() );
  log << MSG::DEBUG << "AcdReconAlg::hitDistances for track " << aTrack.m_index 
      << " going " << (aTrack.m_upward ? "up" : "down" ) <<  " with " << acdHits.size() << " hits." << endreq;

  for (Event::AcdHitCol::const_iterator acdHitIt = acdHits.begin(); acdHitIt != acdHits.end(); acdHitIt++) {
    idents::AcdId acdId = (*acdHitIt)->getAcdId();
    // get the data object to store all the computations
    AcdRecon::PocaData& pocaData = pocaMap[acdId];

    if (acdId.na()) { 
      log << MSG::ERROR << "Skipping NA hit" << endreq;
      continue;
    } else if (acdId.tile()) {      
      const AcdTileDim* tileDim = m_geomMap->getTile(acdId,*m_acdGeoSvc);
      sc = tileDim->statusCode();
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() 
            << endreq;
	return sc;
      }
      if ( m_pocaTool != 0 ) {
	sc = m_pocaTool->tileDistances(*tileDim,aTrack,pocaData);
      } else {
	log << MSG::ERROR << "No Poca Tool" << endreq;
	return StatusCode::FAILURE;
      }
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "Failed to get hit distances for a tile" 
            << acdId.id() << endreq;
	return sc;
      }
    } else if ( acdId.ribbon() ) {
      const AcdRibbonDim* ribbonDim = m_geomMap->getRibbon(acdId,*m_acdGeoSvc);
      sc = ribbonDim->statusCode();
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "Failed to get geom for a ribbon " << acdId.id() 
            << endreq;
	return sc;
      }
      if ( m_pocaTool != 0 ) {
	sc = m_pocaTool->ribbonDistances(*ribbonDim,aTrack,pocaData);
      } else {
	return StatusCode::FAILURE;
      }
      if ( sc.isFailure() ) {
	log << MSG::ERROR << "Failed to get hit distances for a ribbon" 
            << acdId.id() << endreq;
	return sc;
      }
    } else {
      log << MSG::ERROR << "Neither NA, nor tile, nor ribbon" << endreq;
      return StatusCode::FAILURE;
    }
  }
  return sc;
}



StatusCode AcdReconAlg::tileActDist(const AcdRecon::PocaDataMap& pocaMap,
				    std::vector<double> &row_values, 
				    double &return_dist, 
                                    idents::AcdId& maxActDistId){
    // Purpose and Method:  
    // Inputs: 
    // Outputs: 
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // loop on all the pocas we computed
    for ( AcdRecon::PocaDataMap::const_iterator it = pocaMap.begin(); 
          it != pocaMap.end(); it++ ) {
     
        const idents::AcdId& acdId = it->first;
        if (acdId.ribbon()) continue;
        if (acdId.na()) continue;
	const AcdRecon::PocaData& data = it->second;

	// use the original definition of active dist, 2D inside plane, 
        // 3D outside plane
	double dist = data.m_active3D > 0. ? data.m_active2D : data.m_active3D;
	
	if ( dist > return_dist ) {
	  return_dist = dist;
	  maxActDistId = acdId;
	}
	
        // Pick up the max. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && dist > row_values[0]) row_values[0] = dist;
        if (acdId.side()) {
            unsigned int k = acdId.row()+1;
            if( k >= row_values.size()){
                log << MSG::WARNING << "rejecting bad ACD id, row = " 
		    << k-1 << ' ' << row_values.size() << endreq;
            } else
	      if ( data.m_active3D > row_values[k]) row_values[k] = dist;
        }
    }
	
    return sc;
}
		
StatusCode AcdReconAlg::hitRibbonDist(const AcdRecon::PocaDataMap& pocaMap,
				      double &return_dist, idents::AcdId& maxActDistId) {
    // Purpose and Method: 
    // 

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // loop on all the pocas we computed
    for ( AcdRecon::PocaDataMap::const_iterator it = pocaMap.begin(); it != pocaMap.end(); it++ ) {
     
        const idents::AcdId& acdId = it->first;
        if (acdId.tile()) continue;
        if (acdId.na()) continue;
	const AcdRecon::PocaData& data = it->second;
	
	if ( data.m_active2D > return_dist ) {
	  return_dist = data.m_active2D;
	  maxActDistId = acdId;
	}

    }  // end loop over AcdDigis
    return sc;
}

StatusCode AcdReconAlg::extrapolateTrack(const Event::TkrTrack& aTrack,
					 const AcdRecon::TrackData& trackData,
					 const AcdRecon::PocaDataPtrMap& pocaDataMap,
					 const AcdRecon::ExitData& isectData,
					 Event::AcdPocaSet& pocaSet,
					 Event::AcdTkrIntersectionCol& acdIntersections,
					 Event::AcdTkrGapPocaCol& gapPocas,
					 Event::AcdTkrPointCol& points) {

  MsgStream   log( msgSvc(), name() );
  StatusCode sc = StatusCode::SUCCESS;  

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlg::extrapolateTrack(" 
	<< trackData.m_index << ',' << (trackData.m_upward ? "up" : "down" ) << ')' << endreq;
  }

  // first figure out how far to extrapolate track
  double maxArcLength(0.);

  // figure out which direction we are going in
  bool forward = isectData.m_arcLength > 0;
  
  // take the furtherest POCA
  AcdRecon::PocaDataPtrMap::const_iterator itr;
  for ( itr = pocaDataMap.begin(); itr != pocaDataMap.end(); itr++ ) {
    const AcdRecon::PocaData& pocaData = *(itr->second);
    if ( (forward && pocaData.m_arcLength > maxArcLength) || 
	 ((!forward) && pocaData.m_arcLength < maxArcLength ) ) {
      maxArcLength = pocaData.m_arcLength;
    }
  }
  
  // compare that to direct calculation
  double arcToIsect = isectData.m_arcLength;
  arcToIsect += forward ? 100. : -100;  
  if ( (forward && arcToIsect > maxArcLength) || 
       ((!forward) &&  arcToIsect < maxArcLength) ) {
    maxArcLength = arcToIsect;
  }

  // protect against negative arcLengths
  if ( forward && maxArcLength < 0 ) {
    log << MSG::ERROR << "Negative Arclength to upper intersection " << maxArcLength << endreq;
    return sc;
  }
  if ( !forward && maxArcLength > 0 ) {
    log << MSG::ERROR << "Positive Arclength to lower intersection " << maxArcLength << endreq;
    return sc;
  }

  // run the propagator out to the right arclength
  if ( !forward ) maxArcLength *= -1.;

  try {
    AcdRecon::startPropagator(*m_G4PropTool,aTrack,trackData,maxArcLength);
  } catch (...) {
    log << MSG::ERROR << "Caught exception starting propagator on track " << trackData.m_index 
	<< ".  No intersection or POCA's will be calculated for that track." <<   endreq;
    // Don't crash, just continue.
    return sc;
  }  

  // build all the intersections
  if ( m_intersectionTool != 0 ) {
    try {
      sc = m_intersectionTool->makeIntersections(*m_G4PropTool,trackData,isectData,pocaDataMap,m_hitMap,*m_geomMap,
						 acdIntersections,gapPocas);
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
  
  if ( log.level() <= MSG::DEBUG ) {    
    gapPocas.writeOut(log);
  }

  // build all the pocas
  AcdTkrParams paramsAtArcLength;
  Event::TkrTrackParams next_params;

  for ( itr = pocaDataMap.begin(); itr != pocaDataMap.end(); itr++ ) {
    AcdRecon::PocaData& pocaData = *(itr->second);
    float pocaArcLength = forward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    try {
      Vector voca(pocaData.m_voca.x(),pocaData.m_voca.y(),pocaData.m_voca.z());
      AcdRecon::propagateToArcLength(*m_G4PropTool,trackData,pocaArcLength,next_params,paramsAtArcLength);
      AcdRecon::projectErrorToPocaVector(paramsAtArcLength,voca,pocaData.m_active3DErr_proj);
    } catch (...) {
       log << MSG::ERROR << "Caught exception using propagator make POCAs with track " << trackData.m_index 
	<< ".  No ACD POCAs calculated for that track." <<   endreq;
      // Don't crash, just continue.
      return sc;
    }
    Event::AcdTkrHitPoca* aPoca(0);
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->makePoca(trackData,pocaData,aPoca);
    }
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "AcdPocaTool::makePoca failed" << endreq;
      return sc;
    }
    
    if ( aPoca != 0 ) {      
      if ( log.level() <= MSG::DEBUG ) {    
	aPoca->writeOut(log);
      }
      pocaSet.insert(aPoca);
    }
  }

  // build the TrkPoint
  AcdRecon::propagateToArcLength(*m_G4PropTool,trackData,isectData.m_arcLength,next_params,paramsAtArcLength);
  Event::AcdTkrPoint* exitPoint(0);
  if ( m_intersectionTool != 0 ) {
    sc = m_intersectionTool->makeTkrPoint(trackData,isectData,next_params,exitPoint);
    if ( sc.isFailure() ){
      log << MSG::ERROR << "AcdTkrIntersectionTool::makeTkrPoint failed" << endreq;
      return sc;
    }
    if ( exitPoint != 0 ) {
      if ( log.level() <= MSG::DEBUG ) {    
	exitPoint->writeOut(log);
      }
      points.push_back(exitPoint);      
    }
  }

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlg::extrapolateTrack() finished" << std::endl << endreq;
  }

  return sc;

}

StatusCode AcdReconAlg::extrapolateVertex(const AcdRecon::TrackData& trackData,
					  const AcdRecon::PocaDataPtrMap& pocaDataMap,
					  const AcdRecon::ExitData& isectData,
					  Event::AcdPocaSet& pocaSet,
					  Event::AcdTkrPointCol& points) {

  MsgStream   log( msgSvc(), name() );
  StatusCode sc = StatusCode::SUCCESS;  
 
  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlg::extrapolateVertex(" 
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
    AcdRecon::PocaData& pocaData = *(itr->second);
    //float pocaArcLength = forward ? pocaData.m_arcLength : -1* pocaData.m_arcLength;
    //Event::TkrTrackParams next_params = m_G4PropTool->getTrackParams(pocaArcLength,startEnergy,true);
    //AcdRecon::projectErrorAtPoca(trackData,next_params,pocaData.m_poca,pocaData.m_voca,pocaData.m_active3DErr);
    Event::AcdTkrHitPoca* aPoca(0);
    if ( m_pocaTool != 0 ) {
      sc = m_pocaTool->makePoca(trackData,pocaData,aPoca);
    }
    if ( sc.isFailure() ) return sc;
    if ( aPoca != 0 ) {
      aPoca->writeOut(log);
      pocaSet.insert(aPoca);
    }
  }

  // build the TrkPoint
  //Event::TkrTrackParams next_params = m_G4PropTool->getTrackParams(isectData.m_arcLength,startEnergy,true);
  Event::TkrTrackParams next_params;
  Event::AcdTkrPoint* exitPoint(0);
  if ( m_intersectionTool != 0 ) {
    sc = m_intersectionTool->makeTkrPoint(trackData,isectData,next_params,exitPoint);
    if ( sc.isFailure() ) return sc;
    if ( exitPoint != 0 ) {
      exitPoint->writeOut(log);
      points.push_back(exitPoint);      
    }
  }

  if ( log.level() <= MSG::DEBUG ) {
    log << MSG::DEBUG << "AcdReconAlg::extrapolateTrack() finished" << std::endl << endreq;
  }

  return sc;

}

StatusCode AcdReconAlg::calcCornerDoca(const HepPoint3D &x0, const HepVector3D &t0, double &return_dist) {

    return_dist = maxDoca;

    // Use this flag to determine whether to apply sign or not.. if we never
    // find a corner DOCA where the intersection is within the limits of the 
    // side, the return distance should remain at maxDoca
    bool foundOne = false;

    // Form a Ray using the input track
    Point trackPos(x0.x(), x0.y(), x0.z());
    Vector trackDir(t0.x(), t0.y(), t0.z());
    Ray track(trackPos, trackDir);

    // iterate over all corner gaps
    unsigned int iCorner;
    for (iCorner=0; iCorner<4; iCorner++) {
        const Ray gapRay = m_acdGeoSvc->getCornerGapRay(iCorner);
        // Compute DOCA between the track and gap ray 
        AcdUtil::RayDoca doca = AcdUtil::RayDoca(track, gapRay);

        // Check if x,y,z along corner gap ray falls within limits of LAT.
        // where the top is defined as the top ACD tiles and bottom is the
        // base of the bottom row of ACD side tiles
        double length_2_intersect = doca.arcLenRay2();
        if (length_2_intersect > 0 && length_2_intersect < gapRay.getArcLength()) {

            foundOne = true;
            double test_dist = doca.docaRay1Ray2();
            return_dist = (return_dist > test_dist) ? test_dist : return_dist;
        }
    }

    // If no DOCAs were found, just return, and skip the sign calculation
    if(!foundOne) return StatusCode::SUCCESS;

    // Now we have DOCA to the corners
    // Next compute sign based on (Tkr1X0*Tkr1YDir - Tkr1Y0*Tkr1XDir)
    double sign = trackPos.x()*trackDir.y() - trackPos.y()*trackDir.x();
    if (sign < 0) return_dist = -return_dist;
 
    return StatusCode::SUCCESS;
  
}

StatusCode AcdReconAlg::doBacksplash(const Event::AcdDigiCol& /* digiCol */, Event::AcdSplashVarsCol& acdSplashVars) {

  StatusCode sc = StatusCode::SUCCESS;
  MsgStream   log( msgSvc(), name() );
  
  // Retrieve the information on fit tracks
  SmartDataPtr<Event::TkrTrackCol> tracksTds(eventSvc(), 
					     EventModel::TkrRecon::TkrTrackCol);
  
  if (!tracksTds) {
    log << MSG::DEBUG << "No reconstructed tracks found on the TDS" 
	<< endreq;
    return StatusCode::SUCCESS;
  }
  
  AcdRecon::TrackData downwardExtend;
  AcdRecon::SplashData splashData;

  int iTrack(-1);
  Event::TkrTrackColPtr trkPtr = tracksTds->begin();
  while(trkPtr != tracksTds->end()) {

    const Event::TkrTrack* trackTds  = *trkPtr++;       // The TDS track
    iTrack++;

    // only look at downward extension
    const unsigned int lastHitIdx = trackTds->getNumHits() - 1;
    const Event::TkrTrackHit* lastHit = (*trackTds)[lastHitIdx];
    downwardExtend.m_point = lastHit->getPoint(Event::TkrTrackHit::SMOOTHED);
    downwardExtend.m_dir   = lastHit->getDirection(Event::TkrTrackHit::SMOOTHED);
    downwardExtend.m_energy = lastHit->getEnergy();
    downwardExtend.m_index = iTrack;
    downwardExtend.m_upward = false;

    splashData.m_trackIndex = iTrack;
    AcdRecon::entersCal(downwardExtend,0.,splashData.m_calEntryPoint,splashData.m_calEntryVector,splashData.m_region);

    // only keep stuff that hits the CAL
    if ( splashData.m_region != 0 ) continue;

    static std::list<idents::AcdId> acdList;
    if ( acdList.size() == 0 ) {
      acdList.push_back( idents::AcdId(0,0,0,0) );
      acdList.push_back( idents::AcdId(0,1,0,0) );
      acdList.push_back( idents::AcdId(0,1,1,0) );
      acdList.push_back( idents::AcdId(0,1,2,0) );
      acdList.push_back( idents::AcdId(0,1,3,0) );      
    }

    for ( std::list<idents::AcdId>::const_iterator acdIt = acdList.begin(); acdIt != acdList.end(); acdIt++ ) {
    //for (Event::AcdDigiCol::const_iterator acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
      // get the id
      //idents::AcdId acdId = (*acdDigiIt)->getId();
      const idents::AcdId& acdId = *acdIt;
      // get the data object to store all the computations
      if (acdId.na()) continue;
      if (acdId.ribbon()) continue;
      if (acdId.tile()) {      
	const AcdTileDim* tileDim = m_geomMap->getTile(acdId,*m_acdGeoSvc);
	sc = tileDim->statusCode();
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() 
	      << endreq;
	  return sc;
	}
	splashData.resetTileData(2000.);
	AcdRecon::splashVariables(*tileDim,splashData.m_calEntryPoint,splashData.m_calEntryVector,
				  splashData.m_tileSolidAngle, splashData.m_weightedTrackAngle, splashData.m_weightedPathlength );

	Event::AcdSplashVars* newSplash = 
	  new Event::AcdSplashVars(acdId,iTrack,
				   splashData.m_calEntryPoint,splashData.m_calEntryVector,
				   splashData.m_tileSolidAngle , splashData.m_weightedTrackAngle , splashData.m_weightedPathlength);

	acdSplashVars.push_back(newSplash);
      }
    }	
  }
  return sc;
}

void AcdReconAlg::writeExitPoint(std::ostream& os, 
				 const AcdRecon::TrackData& trackData, const AcdRecon::ExitData& isectData) {
  os << "Tk: " 
     << trackData.m_index << " s= " << isectData.m_arcLength << "; P=" << isectData.m_x << "; face=" << isectData.m_face;
}
