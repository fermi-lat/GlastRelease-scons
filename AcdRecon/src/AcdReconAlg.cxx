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

#include "Event/Recon/TkrRecon/TkrTrack.h"

#include "CLHEP/Geometry/Transform3D.h"

#include <algorithm>
#include <cstdio>
#include <stdlib.h>

double AcdReconAlg::s_vetoThresholdMeV;

unsigned int AcdReconAlg::s_numSideRows;

static double maxDoca = 2000.0;



// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor
AcdReconAlg::AcdReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
	
}

// ____________________________________________________________________________

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

// ____________________________________________________________________________

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

// ____________________________________________________________________________

StatusCode AcdReconAlg::finalize() {    
    clear();
    return StatusCode::SUCCESS;
}

// ____________________________________________________________________________

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

// ____________________________________________________________________________

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
    //m_rowActDistCol.resize(s_numSideRows+1, maxDoca); WRONG - Active Distance > 0 is a HIT!!!!!!
    m_rowActDistCol.resize(s_numSideRows+1, -maxDoca);
	m_energyCol.clear();
	m_idCol.clear();
    m_energyRibbonCol.clear();
    m_idRibbonCol.clear();
    m_act_dist = -maxDoca;
    m_ribbon_act_dist = -maxDoca;
}

// ____________________________________________________________________________

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

        if (id.tile()) {
        // toss out hits below threshold -- OLD
        //if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 

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
	
    sc = trackDistances(digiCol);
    if (sc.isFailure()) return sc;
	
    log << MSG::DEBUG;
    if (log.isActive()) log.stream() << "DOCA: " << m_doca << " "
        << "ActDist: " << m_act_dist;
    log << endreq;
	
    SmartDataPtr<Event::AcdRecon> checkAcdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (checkAcdRecTds) {
        log << MSG::DEBUG;
        if (log.isActive()) log.stream() << "AcdRecon data already on TDS!";
        log << endreq;
        checkAcdRecTds->clear();
        checkAcdRecTds->init(m_totEnergy, m_tileCount, m_gammaDoca, m_doca, 
            m_minDocaId, m_act_dist, m_maxActDistId, m_rowDocaCol, 
            m_rowActDistCol, m_idCol, m_energyCol,
            m_ribbon_act_dist, m_ribbon_act_dist_id);
    } else {
        // create the TDS location for the AcdRecon
        Event::AcdRecon *acdRecon = new Event::AcdRecon(m_totEnergy, 
                                   m_tileCount, m_gammaDoca, m_doca, 
                                   m_minDocaId, m_act_dist, m_maxActDistId, 
                           m_rowDocaCol, m_rowActDistCol, m_idCol, m_energyCol, 
                           m_ribbon_act_dist, m_ribbon_act_dist_id);
        
        sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, acdRecon);
        if (sc.isFailure()) {
            log << "Failed to register AcdRecon" << endreq;
            return StatusCode::FAILURE;
        }
    }
    
    return sc;
}

// ____________________________________________________________________________

StatusCode AcdReconAlg::trackDistances(const Event::AcdDigiCol& digiCol) {
    // Purpose and Method:  Retrieves the TkrFitTrackCol from the TDS and 
    //  calculates the DOCA and Active Distance quantities.  Updates the
    // local data members m_doca, m_rowDocaCol, m_act_dist, m_rowActDistCol
    // TDS Input: EventModel::TkrRecon::TkrFitTrackCol
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrTrackCol);
	
    if (!tracksTds) {
        log << MSG::INFO << "No reconstructed tracks found on the TDS" << endreq;
        return StatusCode::SUCCESS;
    }
	
    double testDoca, test_dist, ribDist;
    Event::TkrTrackColPtr trkPtr = tracksTds->begin();
    while(trkPtr != tracksTds->end())
    {
        const Event::TkrTrack* trackTds  = *trkPtr++;       // The TDS track
        sc = doca(digiCol, trackTds->getInitialPosition(), trackTds->getInitialDirection(), 
            m_rowDocaCol, testDoca);
        if (sc.isFailure()) return sc;
        if(testDoca < m_doca) m_doca = testDoca;
        sc = hitTileDist(digiCol, trackTds->getInitialPosition(), 
            -(trackTds->getInitialDirection()), m_rowActDistCol, test_dist);
        if (sc.isFailure()) return sc;
        if(test_dist > m_act_dist) m_act_dist = test_dist;

        sc = hitRibbonDist(digiCol, trackTds->getInitialPosition(), -(trackTds->getInitialDirection()),
            ribDist);
        if (ribDist > m_ribbon_act_dist) m_ribbon_act_dist = ribDist;
        if (sc.isFailure()) return sc;
    }
	
    return sc;
	
}

// ____________________________________________________________________________

StatusCode AcdReconAlg::doca(const Event::AcdDigiCol& digiCol, const HepPoint3D &x0, const HepVector3D &t0, 
						 std::vector<double> &doca_values, double &minDoca) {
    // Purpose and Method:  This method looks for close-by hits to the ACD tiles
    //        Calculates the minimum distance between the track and the center
    //        of all ACD tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: doca_values represent the DOCA for regions of the ACD: top, s0, s1, s2,...
    //          returns minimum DOCA value
    // Dependencies: None
    // Restrictions and Caveats:  None
	
    minDoca = maxDoca;
    double dist;
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
         idents::AcdId acdId = (*acdDigiIt)->getId();
         if (acdId.ribbon()) continue;

        // toss out hits below threshold -- OLD
        // if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue;
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 
		
        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        std::string str;
        std::vector<double> dim;
        sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
            return sc;
        }
        HepTransform3D transform;
        sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::WARNING << "Failed to get transformation" << endreq;
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
	
    return StatusCode::SUCCESS;
}

// ____________________________________________________________________________

StatusCode AcdReconAlg::hitTileDist(const Event::AcdDigiCol& digiCol, const HepPoint3D &x0, const HepVector3D &t0, 
                                std::vector<double> &row_values, double &return_dist) {
    // Purpose and Method:  Bill Atwood's edge DOCA algorithm called active distance
    // Determines minimum distance between a track and the edges of ACD
    // tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: Returns minimum distance
	
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
	
    return_dist = -200.;
	
    // iterate over all tiles
    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        if (acdId.ribbon()) continue;

        // toss out hits below threshold -- OLD
        // if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 
        // Use Veto Discrim instead
        // Skip this ACD detector if neither PMT has veto discrim set
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 


        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();
        std::string str;
        std::vector<double> dim;
        sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
            return sc;
        }
        HepTransform3D transform;
        sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::WARNING << "Failed to get trasnformation" << endreq;
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
        
        HepVector3D local_x0 = x_isec - xT;
        double test_dist;
        if(iFace == 0) {// Top Tile
            double dist_x = dX/2. - fabs(local_x0.x());
            double dist_y = dY/2. - fabs(local_x0.y());	 
            // Choose which is furthest away from edge (edge @ 0.)
            test_dist = (dist_x < dist_y) ? dist_x : dist_y;
            // Choose closest to tile center
            if(test_dist > return_dist) { 
                return_dist = test_dist;
                m_maxActDistId = acdId;
            }
        }
        else if(iFace == 1 || iFace == 3) {// X Side Tile
            double dist_z = dZ/2. - fabs(local_x0.z());
            double dist_y = dX/2. - fabs(local_x0.y());	
            test_dist = (dist_z < dist_y) ? dist_z : dist_y;
            if(test_dist > return_dist) { 
                return_dist = test_dist;
                m_maxActDistId = acdId;          
            }
        }
        else if(iFace == 2 || iFace == 4) {// Y Side Tile
            double dist_z = dZ/2. - fabs(local_x0.z());
            double dist_x = dX/2. - fabs(local_x0.x());
            test_dist = (dist_z < dist_x) ? dist_z : dist_x;
            if(test_dist > return_dist) {
                return_dist = test_dist;
                m_maxActDistId = acdId;
            }
        }
		
        // Pick up the max. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && test_dist > row_values[0]) row_values[0] = test_dist;
        if (acdId.side()) {
            unsigned int k = acdId.row()+1;
            if( k >= row_values.size()){
                log << MSG::WARNING << "rejecting bad ACD id, row = " << k-1 << endreq;
            }else
                if (test_dist > row_values[k]) row_values[k] = test_dist;
        }
		
    }
	
    return StatusCode::SUCCESS;    
}


// ____________________________________________________________________________

StatusCode AcdReconAlg::hitRibbonDist(const Event::AcdDigiCol& digiCol, const HepPoint3D &x0, const HepVector3D &t0, double &return_dist) {
    //Purpose and Method:  Calculate ActiveDistance for Ribbons
    // To simplify the situation - we assume the ribbons are merely lines 
    // Since the ribbons are segmented in the geometry - we treat the side segments and 
    // top segments separately, taking the minimum distance over all hit segments

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    return_dist = -maxDoca;
    int ribbonX = 5, ribbonY = 6;

    // iterate over all digis and search for ribbons
    Event::AcdDigiCol::const_iterator acdDigiIt;
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        idents::AcdId acdId = (*acdDigiIt)->getId();
        // if a tile - skip we want ribbons
        if (acdId.tile()) continue;
        if ( (!(*acdDigiIt)->getVeto(Event::AcdDigi::A)) && (!(*acdDigiIt)->getVeto(Event::AcdDigi::B)) ) continue; 

        // Need to reconstruct the required volume ids to retrieve the geometry information
        // For now we brute force it.. and construct what we need
        int topSegment;
        int sideFace[2];
        // check orientation to determine which segment number to retrieve for top
        if (acdId.ribbonOrientation() == ribbonX) {
            topSegment = 2;
            // ribbons that are along x-axis on the top go down faces 1,3
            sideFace[0] = 1; sideFace[1] = 3;
        } else {
            topSegment = 1;
            // ribbons that are along the y-axis on the top go down faces 2,4
            sideFace[0] = 2; sideFace[1] = 4;
        }

        // toss out hits below threshold
        //if ((*acdDigiIt)->getEnergy() < s_vetoThresholdMeV) continue; 

        // Get the volumeId for this acdId
        idents::VolumeIdentifier volId = (*acdDigiIt)->getVolId();

        // Now loop over 3 segments for this ribbon
        // We want to grab the 2 side segments and retrieve their dimensions
        // Also want to grap one of the top segments and for the orientation that is cut
        // into 5 segments - we want to take the middle segment and use that to construct
        // a line covering the whole top - we'll use that line to calculate active distance
        int isegment;
        for (isegment = 0; isegment < 3; isegment++) {
            // Manually create our volId
            idents::VolumeIdentifier segmentVolId;
            segmentVolId.append(volId[0]);
            short face;
            if (isegment == 1) {
                face = 0;
            } else {
                face = sideFace[isegment/2];
            }
            segmentVolId.append(face);
            segmentVolId.append(volId[2]);
            segmentVolId.append(volId[3]);
            segmentVolId.append(volId[4]);
            if (isegment != 1) { // side ribbon segment
                segmentVolId.append(0); //put back
            } else { // top ribbon segment
                segmentVolId.append(topSegment);
            }

            // Variables for storing ribbon data
            HepPoint3D center;
            double x1=0.0, y1=0.0, z1=0.0;
            double x2=0.0, y2=0.0, z2=0.0;
            double ribbonHalfWidth;

            // in this case, we need to extract the dimensions from 3 other top segments
            // to extend an imaginary ribbon across the whole top of the instrument
            if (acdId.ribbonOrientation() == ribbonX && isegment == 1) {
                int iseg;
                for(iseg = 1; iseg <= 3; iseg++) {
                    idents::VolumeIdentifier volId1;
                    volId1.append(volId[0]); volId1.append(0); volId1.append(volId[2]); 
                    volId1.append(volId[3]); volId1.append(volId[4]);
                    if (iseg == 3) {
                        // grab the last segment
                        volId1.append(5);
                    } else {
                        volId1.append(iseg);
                    }

                    std::vector<double> dim1;
                    sc = getDetectorDimensions(volId1, dim1, center);
                    if (sc.isFailure()) continue;
                    // pick up the beginning from the first segment
                    if (iseg == 1) x1 = center.x() - dim1[0]/2.;
                    if (iseg == topSegment){
                        // pick up the other 2 dimensions from a ribbon in the middle
                        y1 = center.y(); 
                        z1 = center.z();
                        y2 = center.y();
                        z2 = center.z();
                        ribbonHalfWidth = dim1[1]/2.;
                    } else if (iseg == 3) {
                        // Pick up the ending point from the last segment
                        x2 = center.x() + dim1[0]/2.;
                    }
                }
            } else if (acdId.ribbonOrientation() == ribbonY && isegment == 1) {
                std::vector<double> dim;
                sc = getDetectorDimensions(segmentVolId, dim, center);
                if (sc.isFailure()) continue;
                x1 = center.x();
                y1 = center.y() - dim[1]/2.;
                z1 = center.z();
                x2 = center.x();
                y2 = center.y() + dim[1]/2.;
                z2 = center.z();
                ribbonHalfWidth = dim[0]/2.;
            } else { // side ribbons - which are in one segment
                std::vector<double> dim;
                sc = getDetectorDimensions(segmentVolId, dim, center);
                if (sc.isFailure()) continue;
                x1 = center.x();
                y1 = center.y();
                z1 = center.z() - dim[2]/2.;
                x2 = center.x();
                y2 = center.y();
                z2 = center.z() + dim[2]/2.;
                // Determine the half width of the ribbon
                if ((face == 1) || (face == 3)) {
                    ribbonHalfWidth = dim[0]/2.;
                } else {
                    ribbonHalfWidth = dim[0]/2.;
                }
            }

            // After all that, we have the beginning and ending points for a line segment
            // that defines a ribbon
            HepPoint3D ribbonBeginPos(x1, y1, z1);
            HepPoint3D ribbonEndPos(x2, y2, z2);

            int iFace = face; 

            // Figure out where in the plane of this face the trajectory hits
            double arc_dist = -1.; 
            if(iFace == 0) {// Top ribbon
                arc_dist = (z1-x0.z())/t0.z();	                
            } else if(iFace == 1 || iFace == 3) {// X Side ribbon
                arc_dist = (x1-x0.x())/t0.x();
            } else if(iFace == 2 || iFace == 4) {// Y Side ribbon
                arc_dist = (y1-x0.y())/t0.y();
            }

            if (arc_dist < 0) continue;

            // position of hit in the plane of the ribbon
            HepPoint3D x_isec = x0 + arc_dist*t0;

            // Form vector between the beginning of the ribbon and the point where the
            // track intersects the plane of the ribbon
            HepVector3D delta = x_isec - ribbonBeginPos;
            // Form a vector for the ribbon
            HepVector3D ribbonVec = ribbonEndPos - ribbonBeginPos;
            double prod = delta * ribbonVec.unit();
            // check that the projection of the point to the ribbon occurs within the
            // length of the ribbon segment
            if ((prod < 0) || (prod > ribbonVec.mag())) continue;
            double dist = sqrt(delta.mag2() - prod*prod);
            // Make this an Active Distance calculation
            dist = ribbonHalfWidth - dist;
            if(dist > return_dist){
                return_dist = dist;
                m_ribbon_act_dist_id = acdId;
            } // endif
        } // end loop over ribbon segments
    }  // end loop over AcdDigis
    return StatusCode::SUCCESS;    
}

// ____________________________________________________________________________

StatusCode AcdReconAlg::getDetectorDimensions(const idents::VolumeIdentifier &volId, std::vector<double> &dims, HepPoint3D &xT ) {
    MsgStream   log( msgSvc(), name() );
    std::string str;
    StatusCode sc = StatusCode::SUCCESS;
    sc = m_glastDetSvc->getShapeByID(volId, &str, &dims);
    if ( sc.isFailure() ) {
        log << MSG::WARNING << "Failed to retrieve Shape by Id" << endreq;
        return sc;
    }
    HepTransform3D transform;
    sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
    if (sc.isFailure() ) {
        log << MSG::WARNING << "Failed to get trasnformation" << endreq;
        return sc;
    }
    HepPoint3D center(0., 0., 0.);
    xT = transform * center;
    return sc;
}
