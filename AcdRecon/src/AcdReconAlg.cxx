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

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"

#include "Event/Recon/TkrRecon/TkrFitTrackBase.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "xml/IFile.h"

#include <algorithm>
#include <cstdio>
#include <stdlib.h>

double AcdReconAlg::s_vetoThresholdMeV;

unsigned int AcdReconAlg::s_numSideRows;

static double maxDoca = 2000.0;

// Define the factory for this algorithm
static const AlgFactory<AcdReconAlg>  Factory;
const IAlgFactory& AcdReconAlgFactory = Factory;

// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor.
AcdReconAlg::AcdReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
	
}

StatusCode AcdReconAlg::initialize ( ) {
    StatusCode sc = StatusCode::SUCCESS;
	
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
	
    m_glastDetSvc = 0;
    sc = service("GlastDetSvc", m_glastDetSvc, true);
    if (sc.isSuccess() ) {
        sc = m_glastDetSvc->queryInterface(IID_IGlastDetSvc, (void**)&m_glastDetSvc);
    }
	
    if( sc.isFailure() ) {
        log << MSG::ERROR << "AcdReconAlg failed to get the GlastDetSvc" << endreq;
        return sc;
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
	
    m_acdDigiCol = acdDigiCol;
	
    // run the reconstruction
    reconstruct(m_acdDigiCol);

    m_acdDigiCol = 0;
	
    return sc;
}


StatusCode AcdReconAlg::finalize() {    
    clear();
    m_acdDigiCol = 0;
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
    m_acdRecon = 0;
    m_totEnergy = 0.0;
    m_tileCount = 0;
    m_gammaDoca = maxDoca;
    m_doca = maxDoca;
    m_rowDocaCol.clear();
    // one for each side, plus one for the top
    m_rowDocaCol.resize(s_numSideRows+1, maxDoca);  
    m_rowActDistCol.clear();
    // one for each side, plus one for the top
    m_rowActDistCol.resize(s_numSideRows+1, maxDoca);
    m_energyCol.clear();
	m_idCol.clear();
    m_act_dist = -maxDoca;
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
        
        // toss out hits below threshold
        if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        
        m_tileCount++;
        double tileEnergy = (*acdDigiIt)->getEnergy();
        m_totEnergy += tileEnergy;
        idents::AcdId id = (*acdDigiIt)->getId();
        
        // Temporarily populate reconstructed energy collection with digi energy
        m_idCol.push_back(id);
        m_energyCol.push_back(tileEnergy);
    }
	
    log << MSG::DEBUG << "num Tiles = " << m_tileCount << endreq;
    log << MSG::DEBUG << "total energy = " << m_totEnergy << endreq;
	
    trackDistances();
	
    log << MSG::DEBUG << "DOCA: " << m_doca << " "
        << "ActDist: " << m_act_dist << endreq;
	
    SmartDataPtr<Event::AcdRecon> checkAcdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (checkAcdRecTds) {
        log << MSG::DEBUG << "AcdRecon data already on TDS!" << endreq;
        checkAcdRecTds->clear();
        checkAcdRecTds->init(m_totEnergy, m_tileCount, m_gammaDoca, m_doca, 
            m_act_dist, m_minDocaId, m_rowDocaCol, m_rowActDistCol, m_idCol, m_energyCol);
    } else {
        // create the TDS location for the AcdRecon
        m_acdRecon = new Event::AcdRecon(m_totEnergy, m_tileCount, m_gammaDoca, m_doca, 
            m_act_dist, m_minDocaId, m_rowDocaCol, m_rowActDistCol, m_idCol, m_energyCol);
        
        sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, m_acdRecon);
        if (sc.isFailure()) {
            log << "Failed to register AcdRecon" << endreq;
            return StatusCode::FAILURE;
        }
    }
    
    return sc;
}


StatusCode AcdReconAlg::trackDistances() {
    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities.  Updates the
    // local data members m_doca, m_rowDocaCol, m_act_dist, m_rowActDistCol
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCol
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrFitTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrFitTrackCol);
	
    if (!tracksTds) {
        log << MSG::INFO << "No reconstructed tracks found on the TDS" << endreq;
        return StatusCode::SUCCESS;
    }
	
    Event::TkrFitColPtr trkPtr = tracksTds->begin();
    while(trkPtr != tracksTds->end())
    {
        const Event::TkrFitTrackBase* trackTds  = *trkPtr++;       // The TDS track
        double testDoca = doca(trackTds->getPosition(), trackTds->getDirection(), m_rowDocaCol);
        if(testDoca < m_doca) m_doca = testDoca;
        double test_dist= hitTileDist(trackTds->getPosition(), -(trackTds->getDirection()), m_rowActDistCol);
        if(test_dist > m_act_dist) m_act_dist = test_dist;
    }
	
    return sc;
	
}

double AcdReconAlg::doca(const HepPoint3D &x0, const HepVector3D &t0, 
						 std::vector<double> &doca_values) {
    // Purpose and Method:  This method looks for close-by hits to the ACD tiles
    //        Calculates the minimum distance between the track and the center
    //        of all ACD tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: doca_values represent the DOCA for regions of the ACD: top, s0, s1, s2,...
    //          returns minimum DOCA value
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    double minDoca = maxDoca;
    double dist;
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = m_acdDigiCol.begin(); acdDigiIt != m_acdDigiCol.end(); acdDigiIt++) {
		
		// toss out hits below threshold
		if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue;
		
        idents::AcdId acdId = (*acdDigiIt)->getId();
        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        std::string str;
        std::vector<double> dim;
        sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::DEBUG << "Failed to retrieve Shape by Id" << endreq;
            return sc;
        }
        HepTransform3D transform;
        sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::DEBUG << "Failed to get transformation" << endreq;
            return sc;
        }
		
        HepPoint3D center(0., 0., 0.);
        HepPoint3D acdCenter = transform * center;
		
        HepVector3D dX = acdCenter - x0;
		
        double prod = dX * t0;
        dist = sqrt(dX.mag2() - prod*prod);
        if(dist < minDoca){
            minDoca = dist;
            m_minDocaId = acdId;
        }
		
        // Pick up the min. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && dist < doca_values[0]) doca_values[0] = dist;
        if (acdId.side()) {
			unsigned int k = acdId.row()+1;
			if( k >= doca_values.size()){
				log << MSG::WARNING << "rejecting bad ACD id, row = " << k-1 << endreq;
			}else
				if (dist < doca_values[k]) doca_values[k] = dist;
        }
		
    }
	
    return minDoca;
}

double AcdReconAlg::hitTileDist(const HepPoint3D &x0, const HepVector3D &t0, 
                                std::vector<double> &row_values) {
    // Purpose and Method:  Bill Atwood's edge DOCA algorithm called active distance
    //       Determines minimum distance between a track and the edges of ACD
    //       tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: Returns minimum distance
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    double return_dist = -200.;
	
    // iterate over all tiles
    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = m_acdDigiCol.begin(); acdDigiIt != m_acdDigiCol.end(); acdDigiIt++) {
        
        // toss out hits below threshold
		if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        idents::AcdId acdId = (*acdDigiIt)->getId();
        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        std::string str;
        std::vector<double> dim;
        sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::DEBUG << "Failed to retrieve Shape by Id" << endreq;
            return sc;
        }
        HepTransform3D transform;
        sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::DEBUG << "Failed to get trasnformation" << endreq;
            return sc;
        }
		
        HepPoint3D center(0., 0., 0.);
        HepPoint3D xT = transform * center;
        
        int iFace = acdId.face();
		
        double dX = dim[0];
        double dY = dim[1];
        double dZ = dim[2];
        
        // Figure out where in the plane of this face the trajectory hits
        double arc_dist = -1.; 
        if(iFace == 0) {// Top Tile. 
            arc_dist = (xT.z()-x0.z())/t0.z();	                
        }
        else if(iFace == 1 || iFace == 3) {// X Side Tile 
            arc_dist = (xT.x()-x0.x())/t0.x();
        }
        else if(iFace == 2 || iFace == 4) {// Y Side Tile
            arc_dist = (xT.y()-x0.y())/t0.y();
        }
        // If arc_dist is negative... had to go backwards to hit plane... 
        if(arc_dist < 0.) continue;
        
        HepPoint3D x_isec = x0 + arc_dist*t0;
        
        HepVector3D local_x0 = xT - x_isec; 
        double test_dist;
        if(iFace == 0) {// Top Tile
            double dist_x = dX/2. - fabs(local_x0.x());
            double dist_y = dY/2. - fabs(local_x0.y());	                
            test_dist = (dist_x < dist_y) ? dist_x : dist_y;
            if(test_dist > return_dist) return_dist = test_dist;
        }
        else if(iFace == 1 || iFace == 3) {// X Side Tile
            double dist_z = dZ/2. - fabs(local_x0.z());
            double dist_y = dY/2. - fabs(local_x0.y());	                
            test_dist = (dist_z < dist_y) ? dist_z : dist_y;
            if(test_dist > return_dist) return_dist = test_dist;
        }
        else if(iFace == 2 || iFace == 4) {// Y Side Tile
            double dist_z = dZ/2. - fabs(local_x0.z());
            double dist_x = dY/2. - fabs(local_x0.x());	                
            test_dist = (dist_z < dist_x) ? dist_z : dist_x;
            if(test_dist > return_dist) return_dist = test_dist;
        }
		
        // Pick up the min. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && test_dist < row_values[0]) row_values[0] = test_dist;
        if (acdId.side()) {
			unsigned int k = acdId.row()+1;
			if( k >= row_values.size()){
				log << MSG::WARNING << "rejecting bad ACD id, row = " << k-1 << endreq;
			}else
				if (test_dist < row_values[k]) row_values[k] = test_dist;
        }
		
    }
	
    return return_dist;    
}
