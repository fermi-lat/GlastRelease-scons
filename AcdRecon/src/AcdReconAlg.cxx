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
#include "GaudiKernel/ObjectVector.h"


// TkrRecon classes
//#include "TkrRecon/GaudiAlg/TkrReconAlg.h"
//#include "TkrRecon/Track/SiRecObjs.h"
//#include "TkrRecon/Track/GFgamma.h"


#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "geometry/Box.h"

#include "xml/IFile.h"

#include "idents/AcdId.h"

#include <algorithm>
#include <cstdio>
#include <stdlib.h>

double AcdReconAlg::s_threshold_energy;

unsigned int AcdReconAlg::s_numSideRows;

static float maxDOCA = 200.0;

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

    SmartDataPtr<Event::AcdDigiCol> acdDigiColTds(eventSvc(), EventModel::Digi::AcdDigiCol);
    if (!acdDigiColTds) return sc;

    // run the reconstruction
    reconstruct(acdDigiColTds);

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

    sc = m_glastDetSvc->getNumericConstByName("threshold", &s_threshold_energy);

    double temp;
    sc = m_glastDetSvc->getNumericConstByName("numSideRows", &temp);
    s_numSideRows = temp;        
}

void AcdReconAlg::clear() {
    m_totEnergy = 0.0;
    m_tileCount = 0;
    m_gammaDOCA = maxDOCA;
    m_DOCA = maxDOCA;
    m_rowDOCA_vec.clear();
    m_rowDOCA_vec.resize(s_numSideRows+1, maxDOCA);  // one for each side, plus one for the top
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

    Event::AcdDigiCol::const_iterator acdDigiIt;
    
    for (acdDigiIt = digiCol.begin(); acdDigiIt != digiCol.end(); acdDigiIt++) {
        
        if ((*acdDigiIt)->getEnergy() < s_threshold_energy) continue; // toss out hits below threshold

        m_tileCount++;
        double tileEnergy = (*acdDigiIt)->getEnergy();
        m_totEnergy += tileEnergy;
        idents::AcdId id = (*acdDigiIt)->getId();
    }

    log << MSG::DEBUG << "num Tiles = " << m_tileCount << endreq;
    log << MSG::DEBUG << "total energy = " << m_totEnergy << endreq;

    acdTileDoca();

    return sc;
}


StatusCode AcdReconAlg::acdTileDoca() {

    StatusCode sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
/*
    // Retrieve the TKR Reconstruction data
    SmartDataPtr<SiRecObjs> tkrRecData(eventSvc(),"/Event/TkrRecon/SiRecObjs");
    if (tkrRecData == 0) {
        log << MSG::DEBUG << "No TKR Reconstruction available " << endreq;
        return sc;
    }

    std::vector<const GFtrack*> xtracks;
    std::vector<const GFtrack*> ytracks;

    // Now retrieve all of the tracks

    int nparticles = tkrRecData->numParticles();
    if (nparticles > 0) {

        log << MSG::DEBUG << "Number of particles " << nparticles << endreq;

        unsigned int iParticle;
        for (iParticle = 0; iParticle < nparticles; iParticle++) {
            GFparticle *particle = tkrRecData->Particle(iParticle);
            xtracks.push_back(particle->getXGFtrack());
            ytracks.push_back(particle->getYGFtrack());
        }
    } else {
        log << MSG::DEBUG << "No reconstructed particles " << endreq;
    }

    int ngammas = tkrRecData->numGammas();
    log << MSG::DEBUG << "number of gammas = " << ngammas << endreq;

    if (ngammas > 0) {
        unsigned int iGamma;
        // Create a temporary vector to store the DOCA's for the 
        // top and side tiles - we don't want to use the reconstructed
        // gamma direction & vector for these DOCAs...
        std::vector<double> temprowDOCA_vec;
        temprowDOCA_vec.resize(numSideRows+1, maxDOCA); 
        for (iGamma = 0; iGamma < ngammas; iGamma++) {
            GFgamma *gamma = tkrRecData->Gamma(iGamma);
            Point gammaVertex = gamma->vertex();
            Vector gammaDirection = gamma->baseDirection();
            double tempDOCA = DOCA(gammaVertex, gammaDirection, m_rowDOCA_vec);
            if (tempDOCA < m_gammaDOCA) m_gammaDOCA = tempDOCA;
            
            // Store the tracks associated with the gamma
            if (!gamma->getPair(TkrCluster::X)->empty()) 
                xtracks.push_back(gamma->getPair(TkrCluster::X));
            if (!gamma->getBest(TkrCluster::X)->empty()) 
                xtracks.push_back(gamma->getBest(TkrCluster::X));
            if (!gamma->getPair(TkrCluster::Y)->empty())
                ytracks.push_back(gamma->getPair(TkrCluster::Y)); 
            if (!gamma->getBest(TkrCluster::X)->empty())
                ytracks.push_back(gamma->getBest(TkrCluster::Y));
        }
    } else {
        log << MSG::DEBUG << "No reconstructed gammas " << endreq;
    }

    // Now we should have a list of x and y tracks
    // iterate through the list of tracks, and combine to create all possible
    // pairings of X and Y tracks.  The min DOCA is then reported as ACD_DOCA in
    // the ntuple

    for (std::vector<const GFtrack*>::const_iterator xtrk = xtracks.begin(); xtrk != xtracks.end(); xtrk++) {
        for(std::vector<const GFtrack*>::const_iterator ytrk = ytracks.begin(); ytrk != ytracks.end(); ytrk++) {

            const GFtrack* pTrkX = *xtrk;
            const GFtrack* pTrkY = *ytrk;
            
            
            //Vector dir = Vector((*xtrk)->slope(), (*ytrk)->slope(), 1.).unit();
	        //Point x0((*xtrk)->positionAtZ(0), (*ytrk)->positionAtZ(0), 0.);
            Vector xDir   = pTrkX->direction();
            double slopeX = xDir.x();
            Vector yDir   = pTrkY->direction();
            double slopeY = yDir.y();

            Vector dir = Vector(slopeX, slopeY, 1.).unit();
	        Point x0(((*xtrk)->positionAtZ(0)).x(), ((*ytrk)->positionAtZ(0)).y(), 0.);

            float testDOCA = DOCA(x0, dir, m_rowDOCA_vec);
	    if(testDOCA < m_DOCA) m_DOCA = testDOCA;
            
            float test_dist= hitTileDist(x0, dir);
            if(test_dist > m_act_dist) m_act_dist = test_dist;

	}
    }
   
    log << MSG::DEBUG << "ACD_DOCA = " << m_DOCA << endreq;
    log << MSG::DEBUG << "ACD_Act_dist = " << m_act_dist << endreq;
*/
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

    float minDOCA = maxDOCA;
    float dist;
/*
    // iterate over all tiles
    for( IVetoData::const_iterator it= m_AcdData->begin(); it != m_AcdData->end(); ++it) {

        float tile_energy = (*it).energy();
	if(tile_energy < threshold_energy) continue;

        Vector dX = (*it).position() - x0;

        float prod = dX * t0;
        dist = sqrt(dX.mag2() - prod*prod);
        if(dist < minDOCA){
            minDOCA = dist;
            m_hit_tile = *it;
        }
        
        idents::AcdId type = (*it).type();

        // Pick up the min. distance from each type of tile
        // i.e. top, and each type of side row tile
        if (type.top() && dist < doca_values[0]) doca_values[0] = dist;
        if (type.side()) {
            if (dist < doca_values[type.row()+1]) doca_values[type.row()+1] = dist;
        }
    }

    return minDOCA;
    */

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
/*
    double return_dist = -200.;
    
    if(!m_AcdData) return return_dist;
    
    // iterate over all tiles
    for( IVetoData::const_iterator it= m_AcdData->begin(); it != m_AcdData->end(); ++it) {
        
        float eT = (*it).energy();
        if(eT < threshold_energy) continue;
        
        Point xT = (*it).position();
        idents::AcdId tileType((*it).type());
        int iFace = tileType.face();

        float dX = m_tileParams->length((*it).type());
        float dY = m_tileParams->width((*it).type());
        float dZ = m_tileParams->height((*it).type());
        
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
    */

    return 0;
}


