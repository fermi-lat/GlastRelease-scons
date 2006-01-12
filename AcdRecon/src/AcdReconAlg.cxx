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

#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/AcdRecon/AcdTkrIntersection.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "AcdUtil/AcdDetectorList.h"

#include "geometry/Ray.h"
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "./RayDoca.h"

#include <algorithm>
#include <cstdio>
#include <stdlib.h>

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
  
    // Determine the rays for corner gaps once and for all
    m_calcCornerDoca = true;
    sc = m_acdGeoSvc->findCornerGaps();
    if (sc.isFailure()) {
        MsgStream log(msgSvc(), name());
        log << MSG::WARNING << "Could not construct corner gap rays,"
            << " will not calculate AcdCornerDoca" << endreq;
        m_calcCornerDoca = false;
        sc = StatusCode::SUCCESS;
    }

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
	
    SmartDataPtr<Event::AcdDigiCol> acdDigiCol(eventSvc(), EventModel::Digi::AcdDigiCol);
    if (!acdDigiCol) {
        log << MSG::INFO << "No AcdDigiCol found on the TDS" << endreq;
        return sc;
    }
		
    // run the reconstruction
    sc = reconstruct(acdDigiCol);
	
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

    m_energyCol.clear();
    m_idCol.clear();
    m_energyRibbonCol.clear();
    m_idRibbonCol.clear();
    m_act_dist = -maxDoca;
    m_act_dist3D = -maxDoca;  // for new active distance
    m_ribbon_act_dist = -maxDoca;

    m_cornerDoca = maxDoca;

    m_hitMap.clear();

    idents::AcdId resetId;
    resetId.na(1);
    m_minDocaId = resetId;
    m_ribbon_act_dist_id = resetId;
    m_maxActDistId = resetId;
    m_maxActDist3DId = resetId;

}


StatusCode AcdReconAlg::reconstruct (const Event::AcdDigiCol& digiCol) {
    // Purpose and Method:  Actually performs the ACD reconstruction.
    //        Counts the number of hit tiles and determines the total energy deposited
    //        in the ACD.
    // Inputs:  digiCol is a pointer to the TDS ACD detector data.
    // Outputs:  Gaudi StatusCode:  returns StatusCode::FAILURE if an error occurs
    // TDS Output:  EventModel::AcdRecon
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    // reset all member variables to their defaults
    clear();
	
	
    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        
        idents::AcdId id = (*acdDigiIt)->getId();

	// caculate the hitMask and stick it in the map
	unsigned char hitMask = 0;
	
	hitMask |= (*acdDigiIt)->getAcceptMapBit(Event::AcdDigi::A) ? 1 : 0;
	hitMask |= (*acdDigiIt)->getAcceptMapBit(Event::AcdDigi::B) ? 2 : 0;
	hitMask |= (*acdDigiIt)->getVeto(Event::AcdDigi::A) ? 4 : 0;
	hitMask |= (*acdDigiIt)->getVeto(Event::AcdDigi::B) ? 8 : 0;
	hitMask |= (*acdDigiIt)->getCno(Event::AcdDigi::A) ? 16 : 0;
	hitMask |= (*acdDigiIt)->getCno(Event::AcdDigi::B) ? 32 : 0;	
	m_hitMap[id] = hitMask;

        if (id.tile()) {
        // toss out hits below threshold -- OLD
        //if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
#if 0 //THB: for analysis, as opposed to a hardware veto, we want to see *all* tiles with signals
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 
#endif
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
	
    log << MSG::DEBUG;
    if ( log.isActive()) 
        log.stream() << "num Tiles = " << m_tileCount 
            << " total energy = " << m_totEnergy
            << "  numRibbons = " << m_ribbonCount
            << " total Ribbon energy = " << m_totRibbonEnergy;
    log << endreq;

    Event::AcdPocaSet acdPocaSet;
    sc = trackDistances(digiCol,acdPocaSet);
    if (sc.isFailure()) return sc;

    static Event::AcdPocaMap acdPocaMap;
    static Event::AcdTkrPocaCol acdPocas;
    for (  Event::AcdPocaSet::iterator itrPoca = acdPocaSet.begin(); itrPoca != acdPocaSet.end();
	   itrPoca++ ) {
      Event::AcdTkrPoca* sortPoca = const_cast<Event::AcdTkrPoca*>(*itrPoca);
      if ( sortPoca == 0 ) continue;
      acdPocaMap.add(*sortPoca);
      acdPocas.push_back(sortPoca);
    }

    log << MSG::DEBUG;
    if (log.isActive()) log.stream() << "DOCA: " << m_doca << " "
        << "ActDist: " << m_act_dist;
    log << endreq;
	
    
    static Event::AcdTkrIntersectionCol acdIntersections;

    // get pointer to the tracker vertex collection
    SmartDataPtr<Event::TkrTrackCol> trackCol(eventSvc(),EventModel::TkrRecon::TkrTrackCol);
    // if reconstructed tracker data doesn't exist - put the debugging message    
    if (trackCol==0) {
      log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
    } else if (m_intersectionTool != 0) {
      sc = m_intersectionTool->findIntersections(trackCol,&acdIntersections,m_hitMap);
      if ( sc.isFailure() ) {
          log << MSG::WARNING << "AcdIntersectionTool Failed - we'll bravely carry on" << endreq;
          sc = StatusCode::SUCCESS;
      }
    }

    static Event::AcdHitCol acdHits;
    if (m_hitTool != 0) {
      sc = m_hitTool->makeAcdHits(&digiCol,&acdHits);
      if ( sc.isFailure() ) {
	log << MSG::WARNING << "AcdHitTool Failed - we'll bravely carry on" << endreq;
	sc = StatusCode::SUCCESS;
      }
    }

    SmartDataPtr<Event::AcdRecon> checkAcdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (checkAcdRecTds) {
        log << MSG::DEBUG;
        if (log.isActive()) log.stream() << "AcdRecon data already on TDS!";
        log << endreq;
        checkAcdRecTds->clear();
        checkAcdRecTds->init(m_totEnergy, m_totRibbonEnergy, 
			     m_tileCount, m_ribbonCount, 
			     m_gammaDoca, m_doca, m_minDocaId,  
			     m_act_dist, m_maxActDistId, 
			     m_rowDocaCol, m_rowActDistCol, m_idCol, m_energyCol,
			     m_act_dist3D, m_maxActDist3DId, m_rowActDist3DCol,
			     acdIntersections, acdPocas, acdHits, m_ribbon_act_dist, 
                             m_ribbon_act_dist_id, m_cornerDoca);
	// ownership handed to TDS, clear local copies
	acdIntersections.clear();
	acdPocas.clear();
	acdPocaMap.clear();
	acdHits.clear();	
    } else {
        // create the TDS location for the AcdRecon
        Event::AcdRecon *acdRecon = new Event::AcdRecon(m_totEnergy, m_totRibbonEnergy, 
							m_tileCount, m_ribbonCount, 
							m_gammaDoca, m_doca, m_minDocaId, 
							m_act_dist, m_maxActDistId, 
							m_rowDocaCol, m_rowActDistCol, m_idCol, m_energyCol, 
							m_ribbon_act_dist, m_ribbon_act_dist_id, 
							acdIntersections, acdPocas, acdHits,
							m_act_dist3D, m_maxActDist3DId, m_rowActDist3DCol, m_cornerDoca);
	// ownership handed to TDS, clear local copies	
	acdIntersections.clear();
	acdPocas.clear();
	acdPocaMap.clear();
	acdHits.clear();
        sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, acdRecon);
        if (sc.isFailure()) {
            log << "Failed to register AcdRecon" << endreq;
            return StatusCode::FAILURE;
        }
    }

    return sc;
}


StatusCode AcdReconAlg::trackDistances(const Event::AcdDigiCol& digiCol, Event::AcdPocaSet& pocaSet) {
    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities.  Updates the
    // local data members m_doca, m_rowDocaCol, m_act_dist, m_rowActDistCol
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCol
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrTrackCol);
	
    if (!tracksTds) {
        log << MSG::DEBUG << "No reconstructed tracks found on the TDS" << endreq;
        return StatusCode::SUCCESS;
    }
	
    bool firstPassDone = false;
    double testDoca, test_dist, ribDist, newActDist;
    Event::TkrTrackColPtr trkPtr = tracksTds->begin();
    idents::AcdId tempId;

    int iTrack(-1);

    while(trkPtr != tracksTds->end())
    {
        const Event::TkrTrack* trackTds  = *trkPtr++;       // The TDS track
	iTrack++;

        sc = doca(digiCol, *trackTds, m_rowDocaCol, testDoca, tempId);
        if (sc.isFailure()) return sc;
        if(testDoca < m_doca) {
	  m_minDocaId = tempId;
	  m_doca = testDoca;
	}


        // Original Active Distance Calc
        sc = hitTileDist(digiCol, *trackTds, m_rowActDistCol, test_dist, tempId);
        if (sc.isFailure()) return sc;
        if(test_dist > m_act_dist) {
	  m_act_dist = test_dist;
	  m_maxActDistId = tempId;
	}

        // New Active Distance Calc
        sc = tileActDist(digiCol, *trackTds, iTrack, m_rowActDist3DCol, newActDist, tempId, pocaSet);
        if (sc.isFailure()) return sc;
        if(newActDist > m_act_dist3D) {
	  m_act_dist3D = newActDist;
	  m_maxActDist3DId = tempId;
	}

	// "Active Distance" from ribbons
        sc = hitRibbonDist(digiCol, *trackTds, iTrack, ribDist, tempId, pocaSet);
        if (ribDist > m_ribbon_act_dist) {
	  m_ribbon_act_dist = ribDist;
	  m_ribbon_act_dist_id = tempId;
	}
        if (sc.isFailure()) return sc;

        if ((!firstPassDone) && (m_calcCornerDoca)) {
            // take the first track, since it is the best track ignore the rest
            calcCornerDoca(trackTds->getInitialPosition(), 
                          -(trackTds->getInitialDirection()), m_cornerDoca);
            // First track in the list, is the reconstructed gamma
            firstPassDone = true;
            log << MSG::DEBUG << "AcdCornerDoca = " << m_cornerDoca << endreq;
        }
    }
	
    return sc;
	
}


StatusCode AcdReconAlg::doca(const Event::AcdDigiCol& digiCol, 
			     const Event::TkrTrack& aTrack,
                             std::vector<double> &doca_values, double &minDoca, idents::AcdId& minDocaId) {
    // Purpose and Method:  This method looks for close-by hits to the ACD tiles
    //        Calculates the minimum distance between the track and the center
    //        of all ACD tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: doca_values represent the DOCA for regions of the ACD: top, s0, s1, s2,...
    //          returns minimum DOCA value
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    minDoca = maxDoca;
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    AcdPocaTool::PocaData data;

    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        if (acdId.ribbon()) continue;
	
        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        AcdTileDim tileDim(acdId,volId,*m_glastDetSvc);
        sc = tileDim.statusCode();
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() << endreq;
	  return sc;
	}

        // toss out hits below threshold -- OLD
        // if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue;
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
#if 0 //THB: for analysis, as opposed to a hardware veto, we want to see *all* tiles with signals
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 
#endif		

	sc = m_pocaTool->doca(tileDim,aTrack,data);
	if ( data.m_region != 0 && data.m_dist < minDoca ) {
	  minDoca = data.m_dist;
	  minDocaId = acdId;
	}
		
        // Pick up the min. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && data.m_dist < doca_values[0]) doca_values[0] = data.m_dist;
        if (acdId.side()) {
            unsigned int k = acdId.row()+1;
            if( k >= doca_values.size()){
	        log << MSG::WARNING << "rejecting bad ACD id, row = " << k-1 << endreq;
            }else
	      if (data.m_dist < doca_values[k]) doca_values[k] = data.m_dist;
        }
		
    }
	
    return StatusCode::SUCCESS;
}


StatusCode AcdReconAlg::hitTileDist(const Event::AcdDigiCol& digiCol, 
				    const Event::TkrTrack& aTrack,
				    std::vector<double> &row_values, 
				    double &return_dist, idents::AcdId& maxActDistId) {
    // Purpose and Method:  Bill Atwood's edge DOCA algorithm called active distance
    // Determines minimum distance between a track and the edges of ACD
    // tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: Returns minimum distance
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // this is an active dist caculation, we want to pick larger numbers
    return_dist = -maxDoca;
	
    AcdPocaTool::PocaData data;

    // iterate over all digis
    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        if (acdId.ribbon()) continue;

	idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
	AcdTileDim tileDim(acdId,volId,*m_glastDetSvc);
	sc = tileDim.statusCode();
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() << endreq;
	  return sc;
	}

        // toss out hits below threshold -- OLD
        // if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
#if 0 //THB: for analysis, as opposed to a hardware veto, we want to see *all* tiles with signals
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 
#endif

	sc = m_pocaTool->hitTileDist(tileDim,aTrack,data);
	if ( data.m_region != 0 && data.m_dist > return_dist ) {
	  return_dist = data.m_dist;
	  maxActDistId = acdId;
	}
	     
        // Pick up the max. distance from each type of tile
        // i.e. top, and each type of side row tile	
        if (acdId.top() && data.m_region != 0 && data.m_dist > row_values[0]) row_values[0] = data.m_dist;
        if (acdId.side()) {
            unsigned int k = acdId.row()+1;
            if( k >= row_values.size()){
                log << MSG::WARNING << "rejecting bad ACD id, row = " 
                                    << k-1 << endreq;
            } else
	      if (data.m_region != 0 && data.m_dist > row_values[k]) row_values[k] = data.m_dist;
        }
		
    }
	
    return StatusCode::SUCCESS;    
}

StatusCode AcdReconAlg::tileActDist(const Event::AcdDigiCol& digiCol, 
				    const Event::TkrTrack& aTrack, int iTrack,
				    std::vector<double> &row_values, 
				    double &return_dist, idents::AcdId& maxActDistId,
				    Event::AcdPocaSet& pocaSet) {
    // Purpose and Method:  Bill Atwood's edge DOCA algorithm called active distance
    // Determines minimum distance between a track and the edges of ACD
    // tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: Returns minimum distance
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    // this is an active dist caculation, we want to pick larger numbers
    return_dist = -maxDoca;

    AcdPocaTool::PocaData data;
	
    // iterate over all digis
    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        if (acdId.ribbon()) continue;

        // toss out hits below threshold -- OLD
        // if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
#if 0 //THB: for analysis, as opposed to a hardware veto, we want to see *all* tiles with signals
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 
#endif

        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
	AcdTileDim tileDim(acdId,volId,*m_glastDetSvc);
	sc = tileDim.statusCode();
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a tile " << acdId.id() << endreq;
	  return sc;
	}

	sc = m_pocaTool->tileActiveDist(tileDim,aTrack,data);
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get tile active distance " << acdId.id() << endreq;
	  return sc;
	}
	
	if ( data.m_region != 0 && data.m_dist > return_dist ) {
	  return_dist = data.m_dist;
	  maxActDistId = acdId;
	}
	
	Event::AcdTkrPoca* aPoca(0);

	sc = m_pocaTool->makePoca(aTrack,iTrack,data,acdId,*m_G4PropTool,aPoca);
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to make a AcdTkrPoca object" << acdId.id() << endreq;
	  return sc;
	}
	if (aPoca != 0) {
	  pocaSet.insert(aPoca);
	}
	
        // Pick up the max. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && data.m_region != 0 && data.m_dist > row_values[0]) row_values[0] = data.m_dist;
        if (acdId.side()) {
            unsigned int k = acdId.row()+1;
            if( k >= row_values.size()){
                log << MSG::WARNING << "rejecting bad ACD id, row = " 
                                    << k-1 << endreq;
            } else
	      if (data.m_region != 0 && data.m_dist > row_values[k]) row_values[k] = data.m_dist;
        }
    }
	
    return StatusCode::SUCCESS;    
}
		
bool AcdReconAlg::withinTileEdge(const Ray& edge, const HepPoint3D& pos) {
    // Purpose and Method:  Determine if a point occurs within the limits of
    //                      a Ray.  This is done by determining the distance
    //                      between the point and the two "endpoints".  If
    //                      either is greater than the arcLength, return false.
	// Note: This function is not used - its simpler to ask if the arc-length 
	//       along the edge is within the edge...  
 
    HepPoint3D start = edge.position(0);
    double distFromStartToPos = start.distance(pos);
    if (distFromStartToPos > edge.getArcLength()) return false;

    HepPoint3D end = edge.position(edge.getArcLength());
    double distFromEndToPos = end.distance(pos);
    if (distFromEndToPos > edge.getArcLength()) return false;
    return true;
}

StatusCode AcdReconAlg::hitRibbonDist(const Event::AcdDigiCol& digiCol, const Event::TkrTrack& aTrack, int iTrack,
				      double &return_dist, idents::AcdId& maxActDistId, Event::AcdPocaSet& pocaSet) {
    //Purpose and Method:  Calculate ActiveDistance for Ribbons
    // To simplify the situation - we assume the ribbons are merely lines 
    // Since the ribbons are segmented in the geometry - we treat the side segments and 
    // top segments separately, taking the minimum distance over all hit segments

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // this is an active dist caculation, we want to pick larger numbers
    return_dist = -maxDoca;

    AcdPocaTool::PocaData data;

    // iterate over all digis and search for ribbons
    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        // if a tile - skip we want ribbons
        if (acdId.tile()) continue;
#if 0 //THB: for analysis, as opposed to a hardware veto, we want to see *all* ribbons with signals
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 
#endif

	idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
	AcdRibbonDim ribbonDim(acdId,volId,*m_glastDetSvc);
	sc = ribbonDim.statusCode();
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get geom for a ribbon " << acdId.id() << endreq;
	  return sc;
	}

	sc = m_pocaTool->hitRibbonDist(ribbonDim,aTrack,data);	
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to get hitRibbonDist " << acdId.id() << endreq;
	  return sc;
	}

	Event::AcdTkrPoca* aPoca(0);
	sc = m_pocaTool->makePoca(aTrack,iTrack,data,acdId,*m_G4PropTool,aPoca);
	if ( sc.isFailure() ) {
	  log << MSG::ERROR << "Failed to make a AcdTkrPoca object" << acdId.id() << endreq;
	  return sc;
	}
	if (aPoca != 0) {
	  pocaSet.insert(aPoca);
	}

	if ( data.m_region != 0 && data.m_dist > return_dist ) {
	  return_dist = data.m_dist;
	  maxActDistId = acdId;
	}

    }  // end loop over AcdDigis
    return StatusCode::SUCCESS;    
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
        RayDoca doca = RayDoca(track, gapRay);

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
