// File and Version Information:
//      $Header$
//
// Description:
//      AcdReconAlg is a Gaudi algorithm which performs the ACD reconstruction.
//          
// Author(s):
//      Heather Kelly			

#include "AcdReconAlg.h"

// Gaudi specific include files
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/StatusCode.h"


// TkrRecon classes
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrPatCandCol.h"


#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/Box.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "xml/IFile.h"

#include "idents/AcdId.h"

#include <algorithm>
#include <cstdio>
#include <stdlib.h>

double AcdReconAlg::s_thresholdEnergy;

unsigned int AcdReconAlg::s_numSideRows;

static float maxDoca = 200.0;

// Define the factory for this algorithm
static const AlgFactory<AcdReconAlg>  Factory;
const IAlgFactory& AcdReconAlgFactory = Factory;

// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor.
AcdReconAlg::AcdReconAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {

}

StatusCode AcdReconAlg::initialize ( )
{

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
    // Inputs:  None
    // Outputs:  Gaudi StatusCode
    // Dependencies: None
    // Restrictions and Caveats:  None
    
    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    log << MSG::DEBUG << "execute" << endreq;

    SmartDataPtr<Event::AcdDigiCol> acdDigiCol(eventSvc(), EventModel::Digi::AcdDigiCol);
    if (!acdDigiCol) return sc;

    m_acdDigiCol = acdDigiCol;


    // run the reconstruction
    reconstruct(m_acdDigiCol);

    // We are now ready to end the routine.
    // Do not delete any allocated memory that has been registered with the TDS - that memory now
    // belongs to the TDS and will be cleaned up automatically at the end of the event.

    return sc;
}


StatusCode AcdReconAlg::finalize() {
    
    MsgStream log(msgSvc(), name());
    
    return StatusCode::SUCCESS;
}


void AcdReconAlg::getParameters ()
{
    MsgStream   log( msgSvc(), name() );
    StatusCode sc;

    sc = m_glastDetSvc->getNumericConstByName("threshold", &s_thresholdEnergy);

    double temp;
    sc = m_glastDetSvc->getNumericConstByName("numSideRows", &temp);
    s_numSideRows = temp;        
}

void AcdReconAlg::clear() {
    m_acdRecon = 0;
    m_totEnergy = 0.0;
    m_tileCount = 0;
    m_gammaDoca = maxDoca;
    m_doca = maxDoca;
    m_rowDocaCol.clear();
    m_rowDocaCol.resize(s_numSideRows+1, maxDoca);  // one for each side, plus one for the top
    m_energyCol.clear();
    m_act_dist = -200.0;
}

StatusCode AcdReconAlg::reconstruct (const Event::AcdDigiCol& digiCol)
{
    // Purpose and Method:  Actually performs the ACD reconstruction.
    //        Counts the number of hit tiles and determines the total energy deposited
    //        in the ACD.
    // Inputs:  v is a pointer to the TDS ACD detector data.
    // Outputs:  Gaudi StatusCode:  returns StatusCode::FAILURE if an error occurs
    // Dependencies: None
    // Restrictions and Caveats:  None

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    // reset all member variables to their defaults
    clear();

    // create the TDS location for the AcdRecon
    m_acdRecon = new Event::AcdRecon;
    sc = eventSvc()->registerObject(EventModel::AcdRecon::Event, m_acdRecon);
    if (sc.isFailure()) {
        log << "Failed to register AcdRecon" << endreq;
        return StatusCode::FAILURE;
    }

    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        
        if ((*acdDigiIt)->getEnergy() < s_thresholdEnergy) continue; // toss out hits below threshold

        m_tileCount++;
        double tileEnergy = (*acdDigiIt)->getEnergy();
        m_totEnergy += tileEnergy;
        idents::AcdId id = (*acdDigiIt)->getId();
    }

    log << MSG::DEBUG << "num Tiles = " << m_tileCount << endreq;
    log << MSG::DEBUG << "total energy = " << m_totEnergy << endreq;

    acdDoca();

    log << MSG::DEBUG << "DOCA: " << m_doca << " "
        << "ActDist: " << m_act_dist << endreq;

    m_acdRecon->initialize(m_totEnergy, m_tileCount, m_gammaDoca, m_doca, 
        m_act_dist, m_minDocaId, m_rowDocaCol, m_energyCol);

    return sc;
}


StatusCode AcdReconAlg::acdDoca() {

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrFitTrackCol> tracksTds(eventSvc(), EventModel::TkrRecon::TkrFitTrackCol);

    if (!tracksTds) return StatusCode::SUCCESS;

    Event::TkrFitColPtr trkPtr = tracksTds->begin();
    while(trkPtr != tracksTds->end())
    {
        Event::TkrFitTrack* trackTds  = *trkPtr++;       // The TDS track
        float testDoca = doca(trackTds->getPosition(), trackTds->getDirection(), m_rowDocaCol);
        if(testDoca < m_doca) m_doca = testDoca;
        float test_dist= hitTileDist(trackTds->getPosition(), -(trackTds->getDirection()));
        if(test_dist > m_act_dist) m_act_dist = test_dist;

    }


    return sc;

}

double AcdReconAlg::doca (const Point &x0, const Vector &t0, std::vector<double> &doca_values)
{
    // Purpose and Method:  This section looks for close-by hits to the ACD tiles
    //        Calculates the minimum distance between the track and the center
    //        of all ACD tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: doca_values represent the DOCA for regions of the ACD: top, s0, s1, s2,...
    //          returns minimum DOCA value
    // Dependencies: None
    // Restrictions and Caveats:  None

    float minDoca = maxDoca;
    float dist;
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = m_acdDigiCol.begin(); acdDigiIt != m_acdDigiCol.end(); acdDigiIt++) {
        
        if ((*acdDigiIt)->getEnergy() < s_thresholdEnergy) continue; // toss out hits below threshold

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
        HepPoint3D acdCenter = transform * center;

        Vector dX = acdCenter - x0;

        float prod = dX * t0;
        dist = sqrt(dX.mag2() - prod*prod);
        if(dist < minDoca){
            minDoca = dist;
            m_minDocaId = acdId;
        }

        // Pick up the min. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (acdId.top() && dist < doca_values[0]) doca_values[0] = dist;
        if (acdId.side()) {
            if (dist < doca_values[acdId.row()+1]) doca_values[acdId.row()+1] = dist;
        }

        return minDoca;


    }

    return 0;
}

double AcdReconAlg::hitTileDist(const Point &x0, const Vector &t0)
{
    // Purpose and Method:  Bill Atwood's new edge DOCA algorithm
    //       Determines minimum distance between a track and the edges of ACD
    //       tiles above veto threshold.
    // Inputs:  (x0, t0) defines a track
    // Outputs: Returns minimum distance
    // Dependencies: None
    // Restrictions and Caveats:  None
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );

    double return_dist = -200.;
        
    // iterate over all tiles
    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = m_acdDigiCol.begin(); acdDigiIt != m_acdDigiCol.end(); acdDigiIt++) {
        
        if ((*acdDigiIt)->getEnergy() < s_thresholdEnergy) continue; // toss out hits below threshold
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

        float dX = dim[0];
        float dY = dim[1];
        float dZ = dim[2];
        
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
        
        Point x_isec = x0 + arc_dist*t0;
        
        Vector local_x0 = xT - x_isec; 
        if(iFace == 0) {// Top Tile
            double dist_x = dX/2. - fabs(local_x0.x());
            double dist_y = dY/2. - fabs(local_x0.y());	                
            double test_dist = (dist_x < dist_y) ? dist_x : dist_y;
            if(test_dist > return_dist) return_dist = test_dist;
        }
        else if(iFace == 1 || iFace == 3) {// X Side Tile
            double dist_z = dZ/2. - fabs(local_x0.z());
            double dist_y = dY/2. - fabs(local_x0.y());	                
            double test_dist = (dist_z < dist_y) ? dist_z : dist_y;
            if(test_dist > return_dist) return_dist = test_dist;
        }
        else if(iFace == 2 || iFace == 4) {// Y Side Tile
            double dist_z = dZ/2. - fabs(local_x0.z());
            double dist_x = dY/2. - fabs(local_x0.x());	                
            double test_dist = (dist_z < dist_x) ? dist_z : dist_x;
            if(test_dist > return_dist) return_dist = test_dist;
        }
    }
    return return_dist;
    
}


