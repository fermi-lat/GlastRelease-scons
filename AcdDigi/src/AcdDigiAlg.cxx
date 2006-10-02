#define AcdDigi_AcdDigiAlg_CXX

// File and Version Information:
// $Header$
// Description:
// Implementation of the latest digitization algorithm for the ACD where
// the Monte Carlo hit information is assumed to be stored in McPositionHits.

#include "AcdDigiAlg.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Geometry/Transform3D.h"

// for min and floor functions
#include <algorithm>
#include <cmath>

#include <utility>
#include <map>


// Define the factory for this algorithm
static const AlgFactory<AcdDigiAlg>  Factory;
const IAlgFactory& AcdDigiAlgFactory = Factory;


// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor.
AcdDigiAlg::AcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    
    declareProperty ("xmlFile", m_xmlFile="$(ACDDIGIROOT)/xml/acdDigi.xml");

    declareProperty("autoCalibrate", m_auto_calibrate=true);

    declareProperty("applyPoisson", m_apply_poisson=true);

    declareProperty("applyGaussianNoise", m_apply_noise=true);

    declareProperty("edgeEffect", m_edge_effect=true);
}


StatusCode AcdDigiAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    StatusCode  sc = StatusCode::SUCCESS;
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    getParameters();

    declareProperty("lowThreshold", m_low_threshold_mips = m_low_threshold_mips_xml);
    // Call setProperies again to grab the low threshold if it is available
    // in the JO file - we needed the name of the input xml file first
    setProperties();
    log << MSG::DEBUG << "Set LowThreshold = " << m_low_threshold_mips << endreq; 

    // read in the parameters from our input XML file
    util.getParameters(m_xmlFile);

    m_glastDetSvc = 0;
    sc = service("GlastDetSvc", m_glastDetSvc, true);
    if (sc.isSuccess() ) {
        sc = m_glastDetSvc->queryInterface(IID_IGlastDetSvc, (void**)&m_glastDetSvc);
    }
   
    if( sc.isFailure() ) {
        log << MSG::ERROR << "AcdDigiAlg failed to get the GlastDetSvc" << endreq;
        return sc;
    }

    // get the list of layers, to be used to add noise to otherwise empty layers
    m_tiles.setPrefix(m_glastDetSvc->getIDPrefix());
    
    // Find all the ACD detectors in our geometry
    m_glastDetSvc->accept(m_tiles);
    if (m_tiles.size() > 0) 
        log << MSG::INFO << "Located  " << m_tiles.size() 
            << " ACD volumes, ids from " << m_tiles.front().name() 
            << " to " << m_tiles.back().name() << endreq;

    
    for(AcdTileList::const_iterator it=m_tiles.begin(); it!=m_tiles.end(); 
                                        ++it) {
        std::string str;
        std::vector<double> dim;
        idents::VolumeIdentifier volId = *it;
        idents::AcdId tileId(volId);
        int iFace = volId[1];

        // Keep a count of volumes associated with each AcdId
        // This will help us determine bent tiles later during edgeEffects
        if (m_acdId_volCount.find(tileId) != m_acdId_volCount.end()) {
            m_acdId_volCount[tileId] += 1;
        } else {
           m_acdId_volCount[tileId] = 1;
        }

        // retrieve the dimensions of this volume from the GlastDetSvc
        sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
        if ( sc.isFailure() ) {
            log << MSG::DEBUG << "Failed to retrieve Shape by Id" << endreq;
            return sc;
        }
        HepGeom::Transform3D transform;
        sc = m_glastDetSvc->getTransform3DByID(volId, &transform);
        if (sc.isFailure() ) {
            log << MSG::WARNING << "Failed to get transformation" << endreq;
            return sc;
         }

         HepPoint3D center(0., 0., 0.);
         HepPoint3D acdCenter = transform * center;

        log << MSG::DEBUG << "VolId " << volId.name() << " AcdId " 
            << tileId.id() << endreq;
        log << MSG::DEBUG << "Dimensions:  " << dim[0] << " " << dim[1] << " " 
            << dim[2] << endreq;
        log << MSG::DEBUG << "Center: " << acdCenter.x() << " " << acdCenter.y()
            << " " << acdCenter.z() << endreq;

    }

    log << MSG::INFO << "Located  " << m_acdId_volCount.size() 
        << " ACD detectors" << endreq;

    return StatusCode::SUCCESS;
}


StatusCode AcdDigiAlg::execute() {
    // Purpose and Method:  Using the McPositionHits that hit the ACD tiles
    //   construct the AcdDigi collection.
    // TDS Input:  EventModel::MC::McPositionHitCol
    // TDS Output: EventModel::Digi::AcdDigiCol

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    using namespace Event;

    clear();
    
    SmartDataPtr<Event::McPositionHitCol> allhits(eventSvc(),EventModel::MC::McPositionHitCol );

    if (!allhits) {
        log << MSG::INFO << "No McPositionHits were found in the TDS" << endreq;
		if (!m_apply_noise) return sc;
		// Otherwise, Do not return - there may be noise hits
    }

    //Take care of insuring that data area has been created
    DataObject* pNode = 0;
    sc = eventSvc()->retrieveObject( EventModel::Digi::Event , pNode);
    
    if (sc.isFailure()) {
        sc = eventSvc()->registerObject(EventModel::Digi::Event, 
                                        new Event::DigiEvent);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register " 
                << EventModel::Digi::Event << endreq;
            return sc;
        }
    }

    std::map<idents::AcdId, double> energyIdMap;
    
    // Create the new AcdDigi collection for the TDS
    Event::AcdDigiCol* digiCol = new Event::AcdDigiCol;
    
    // loop over hits, skip if hit is not in ACD
    // Accumulate the deposited energies, applying edge effects if requested
    if (allhits) {
        for (Event::McPositionHitVector::const_iterator hit = allhits->begin();
              hit != allhits->end(); hit++) {
        
            idents::VolumeIdentifier volId = 
                                ((idents::VolumeIdentifier)(*hit)->volumeID());
            // Check to see if this is an ACD volume
            if(volId[0] != 1 ) continue; 

            idents::AcdId id(volId);
            double energy;
            // No edge effects for ribbons, tiles only
            energy = (m_edge_effect && id.tile()) ? edgeEffect(*hit) : 
                                                    (*hit)->depositedEnergy();
            log << MSG::DEBUG << "AcdId " << id.id() << " E: " 
                << energy << endreq;
            if (energyIdMap.find(id) != energyIdMap.end()) {
                energyIdMap[id] += energy;
            } else {
                energyIdMap[id] = energy;
            }
         } 
    }

    // Add noise to all tiles if requested
    if (m_apply_noise) addNoise();

    // Now loop over the map of AcdId and their corresponding energies deposited
    std::map<idents::AcdId, double>::const_iterator acdIt;
    for (acdIt = energyIdMap.begin(); acdIt != energyIdMap.end(); acdIt++) {

        idents::AcdId id = acdIt->first; //(volId);

        double energyMevDeposited = acdIt->second;
        m_energyDepMap[id] = energyMevDeposited;

        log << MSG::DEBUG << "tile id found: " << id.id() 
            << ", energy deposited: "<< energyMevDeposited<< " MeV" << endreq;
        
        double mips = util.convertMevToMips(energyMevDeposited);
                
        // apply the signal to the two PMT's
        double pmtA_mips = mips;
        double pmtB_mips = mips;
        
        // Number of photoelectrons for each PMT, A and B
        unsigned int pmtA_pe, pmtB_pe;
        util.convertMipsToPhotoElectrons(id, pmtA_mips, pmtA_pe, 
                                         pmtB_mips, pmtB_pe);
        
        // Apply Poisson fluctuations to the number of pe's for each PMT
        if (m_apply_poisson) {
            pmtA_pe = util.shootPoisson(pmtA_pe);
            pmtB_pe = util.shootPoisson(pmtB_pe);
        }
        
        util.convertPhotoElectronsToMips(id, pmtA_pe, pmtA_mips, 
                                             pmtB_pe, pmtB_mips);
        
        double pmtA_mipsToFullScale, pmtB_mipsToFullScale;
        
        // If in auto calibrate mode, determine the conversion factor from MIPs
        // to PHA now, at runtime
        if (m_auto_calibrate) {
            util.calcMipsToFullScale(id, pmtA_mips, pmtA_pe, 
                pmtA_mipsToFullScale, pmtB_mips, pmtB_pe, pmtB_mipsToFullScale);
        } else {
            util.applyGains(id, pmtA_mipsToFullScale, pmtB_mipsToFullScale);
        }
        
        m_pmtA_toFullScaleMap[id] = pmtA_mipsToFullScale;
        m_pmtB_toFullScaleMap[id] = pmtB_mipsToFullScale;

        double pmtA_observedMips_pha = pmtA_mips;
        double pmtA_observedMips_veto = pmtA_mips;
        double pmtA_observedMips_cno = pmtA_mips;
        
        double pmtB_observedMips_pha = pmtB_mips;
        double pmtB_observedMips_veto = pmtB_mips;
        double pmtB_observedMips_cno = pmtB_mips;
        
        // If Gaussian noise is requested it, apply it here.
        // Noise is applied separately for pha, veto and CNO discriminators
        if (m_apply_noise) { 
            pmtA_observedMips_pha += m_pmtA_phaMipsMap[id];
            if (pmtA_observedMips_pha < 0.) pmtA_observedMips_pha = 0.;
            pmtA_observedMips_veto += m_pmtA_vetoMipsMap[id];
            if (pmtA_observedMips_veto < 0.) pmtA_observedMips_veto = 0.;
            pmtA_observedMips_cno += m_pmtA_cnoMipsMap[id];
            if (pmtA_observedMips_cno < 0.) pmtA_observedMips_cno = 0.;
            
            pmtB_observedMips_pha += m_pmtB_phaMipsMap[id];
            if (pmtB_observedMips_pha < 0.) pmtB_observedMips_pha = 0.;
            pmtB_observedMips_veto += m_pmtB_vetoMipsMap[id];
            if (pmtB_observedMips_veto < 0.) pmtB_observedMips_veto = 0.;
            pmtB_observedMips_cno += m_pmtB_cnoMipsMap[id];
            if (pmtB_observedMips_cno < 0.) pmtB_observedMips_cno = 0.;
        }
        
        // Save pha, veto, and CNO MIPs for use later
        m_pmtA_phaMipsMap[id] = pmtA_observedMips_pha;
        m_pmtB_phaMipsMap[id] = pmtB_observedMips_pha;

        m_pmtA_vetoMipsMap[id] = pmtA_observedMips_veto;
        m_pmtB_vetoMipsMap[id] = pmtB_observedMips_veto;

        m_pmtA_cnoMipsMap[id] = pmtA_observedMips_cno;
        m_pmtB_cnoMipsMap[id] = pmtB_observedMips_cno;
        
    } // end loop over MC hits


    std::map<idents::AcdId, bool> doneMap;

    // Now fill the TDS with AcdDigis
    // Loop over all tiles in the geometry
    for(AcdTileList::const_iterator it=m_tiles.begin(); it!=m_tiles.end(); ++it){
        idents::VolumeIdentifier volId = *it;
        idents::AcdId tileId(volId);

        // Check to see if we have already processed this AcdId - 
        // necessary since we also have ribbons in the mix and curved tiles
        // have muliple pieces (volumes)
        if (doneMap.find(tileId) != doneMap.end()) continue;
        
        // First check to see if this tile Id has a PHA value associated with 
        // PMT A. If not, we do not need add an AcdDigi in the TDS for this tile
        if (m_pmtA_phaMipsMap.find(tileId) == m_pmtA_phaMipsMap.end()) continue;

        // Next check that the PHA values from both PMTs combined results is a i
        // value above low threshold
        bool lowThresh = true, vetoThresh = true, cnoThresh = true;
        if ((m_pmtA_phaMipsMap[tileId] < m_low_threshold_mips) && 
            (m_pmtB_phaMipsMap[tileId] < m_low_threshold_mips)) 
                 lowThresh = false;
        if ((m_pmtA_vetoMipsMap[tileId] < m_veto_threshold_mips) && 
            (m_pmtB_vetoMipsMap[tileId] < m_veto_threshold_mips)) 
                vetoThresh = false;
        if ((m_pmtA_cnoMipsMap[tileId] < m_high_threshold_mips) && 
            (m_pmtB_cnoMipsMap[tileId] < m_high_threshold_mips)) 
                cnoThresh = false;

        // If neither PMT is above any threshold - skip this one
        if (!lowThresh && !vetoThresh && !cnoThresh) continue;

        // Initialize discriminators
        bool lowArr[2] = { false, false };
        bool vetoArr[2] = { false, false };
        bool highArr[2] = { false, false };
	Event::AcdDigi::Range rangeArr[2] = { Event::AcdDigi::LOW, Event::AcdDigi::LOW };
        
        // Set the discriminators if above threshold
        if (m_pmtA_phaMipsMap[tileId] > m_low_threshold_mips) lowArr[0] = true;
        if (m_pmtB_phaMipsMap[tileId] > m_low_threshold_mips) lowArr[1] = true;

        if (m_pmtA_vetoMipsMap[tileId] > m_veto_threshold_mips) 
            vetoArr[0] = true;
        if (m_pmtB_vetoMipsMap[tileId] > m_veto_threshold_mips) 
            vetoArr[1] = true;
        
        if (m_pmtA_cnoMipsMap[tileId] > m_high_threshold_mips) 
            highArr[0] = true;
        if (m_pmtB_cnoMipsMap[tileId] > m_high_threshold_mips) 
            highArr[1] = true;

        // Now convert MIPs into PHA values for each PMT

        unsigned short pmtA_pha = util.convertMipsToPha(m_pmtA_phaMipsMap[tileId], m_pmtA_toFullScaleMap[tileId],
							rangeArr[0]);
        unsigned short pmtB_pha = util.convertMipsToPha(m_pmtB_phaMipsMap[tileId], m_pmtB_toFullScaleMap[tileId],
							rangeArr[1]);

        unsigned short phaArr[2] = { pmtA_pha, pmtB_pha };

        doneMap[tileId] = true;	

	AcdDigi* aDigi = new AcdDigi(tileId, volId,
				     m_energyDepMap[tileId], phaArr, 
				     vetoArr, lowArr, highArr);

    log << MSG::DEBUG << "making ranges " << rangeArr[0] << ' ' << rangeArr[1] << endreq;
	aDigi->setRanges(rangeArr);
        digiCol->push_back( aDigi );

    }

    clear();

    // Put the AcdDigi collection on the TDS
    return eventSvc()->registerObject(EventModel::Digi::AcdDigiCol, digiCol);    
}


StatusCode AcdDigiAlg::finalize() {
    MsgStream   log( msgSvc(), name() );
    log << MSG::DEBUG;
    if (log.isActive()) util.dumpMeanPePerPmt(log.stream());
    log << endreq;
    return StatusCode::SUCCESS;
}

void AcdDigiAlg::clear() {

    m_energyDepMap.clear();
    m_pmtA_toFullScaleMap.clear();
    m_pmtA_phaMipsMap.clear();
    m_pmtA_vetoMipsMap.clear();
    m_pmtA_cnoMipsMap.clear();
    m_pmtB_toFullScaleMap.clear();
    m_pmtB_phaMipsMap.clear();
    m_pmtB_vetoMipsMap.clear();
    m_pmtB_cnoMipsMap.clear();

    m_acdId_volCount.clear();

}


void AcdDigiAlg::getParameters() {    
    // Purpose and Method:  Read in the parameters from the XML file using IFile.
    //   Perform some spot checking to see if the parameters exist in the input 
    //   XML file - if not provide an informational message.

    MsgStream log(msgSvc(), name());
 
    xmlBase::IFile xmlFilePtr(m_xmlFile.c_str());
    
    // Perform some spot checking to see if our XML file contains our constants
    if (!xmlFilePtr.contains("thresholds", "low_threshold_mips")) {
        log << MSG::INFO << "XML file " << m_xmlFile
            << " does not contain low_threshold_mips using default" << endreq;
    }

    m_low_threshold_mips_xml = xmlFilePtr.getDouble(
                                   "thresholds", "low_threshold_mips", 0.1);
    m_veto_threshold_mips = xmlFilePtr.getDouble(
                                   "thresholds", "veto_threshold_mips", 0.3);
    m_high_threshold_mips = xmlFilePtr.getDouble(
                                   "thresholds", "high_threshold_mips", 20.0);
    
    if (!xmlFilePtr.contains("global_constants", "mean_pe_per_mip")) {
        log << MSG::INFO << "XML file " << m_xmlFile
            << "does not contain mean_pe_per_mip using default" << endreq;
    }

    m_mean_pe_per_mip = xmlFilePtr.getInt(
                            "global_constants", "mean_pe_per_mip", 18);
    
    m_noise_std_dev_pha = xmlFilePtr.getDouble(
                              "global_constants", "noise_std_dev_pha", 0.02);
    m_noise_std_dev_veto = xmlFilePtr.getDouble(
                               "global_constants", "noise_std_dev_veto", 0.02);
    m_noise_std_dev_cno = xmlFilePtr.getDouble(
                              "global_constants", "noise_std_dev_cno", 0.02);
    
    m_full_scale = xmlFilePtr.getInt("global_constants", "full_scale", 4095);
    
    m_mips_full_scale = xmlFilePtr.getDouble(
                            "global_constants", "mips_full_scale", 20.0);
    
    m_mev_per_mip = xmlFilePtr.getDouble(
                        "global_constants", "mev_per_mip", 1.9);

    m_max_edge_dist = xmlFilePtr.getDouble(
                          "edge_effects", "max_edge_dist", 20.0);

    m_edge_slope = xmlFilePtr.getDouble("edge_effects", "edge_slope", 0.01);

    m_edge_intercept = xmlFilePtr.getDouble(
                           "edge_effects", "edge_intercept", 0.8);

    return;
}


void AcdDigiAlg::addNoise()  {
    // Purpose and Method:  add noise to PMTs
    // Inputs:  AcdDigiAlg::m_tiles
    // Outputs: additional hits in the ACD subsystem
    // TDS Inputs:  None 
    // TDS Outputs:  None
    // Dependencies: None
    // Restrictions and Caveats:  None
        
    std::map<idents::AcdId, bool> doneMap;

    // loop over list of possible tile ids
    for(AcdTileList::const_iterator it=m_tiles.begin(); it!=m_tiles.end(); ++it){
        idents::VolumeIdentifier volId = *it;
        idents::AcdId tileId(volId);
        
        // Check to see if we have processed this AcdId already
        // due to presence of ribbons and multiple volumes for curved tiles
        if (doneMap.find(tileId) != doneMap.end()) continue;
        
        
        m_pmtA_phaMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_pha);
        m_pmtA_vetoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_veto);
        m_pmtA_cnoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_cno);
        
        m_pmtB_phaMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_pha);
        m_pmtB_vetoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_veto);
        m_pmtB_cnoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_cno);

        // Number of photoelectrons for each PMT, A and B
        unsigned int pmtA_pe, pmtB_pe;
        util.convertMipsToPhotoElectrons(tileId, m_pmtA_phaMipsMap[tileId], 
                                pmtA_pe, m_pmtB_phaMipsMap[tileId], pmtB_pe);
        // If in auto calibrate mode, determine the conversion factor from MIPs
        // to PHA now, at runtime
        double pmtA_mipsToFullScale, pmtB_mipsToFullScale;
        if (m_auto_calibrate) {
            util.calcMipsToFullScale(tileId, m_pmtA_phaMipsMap[tileId], 
                pmtA_pe, pmtA_mipsToFullScale, m_pmtB_phaMipsMap[tileId], 
                pmtB_pe, pmtB_mipsToFullScale);
        } else {
            util.applyGains(tileId, pmtA_mipsToFullScale, pmtB_mipsToFullScale);
        }
        m_pmtA_toFullScaleMap[tileId] = pmtA_mipsToFullScale;
        m_pmtB_toFullScaleMap[tileId] = pmtB_mipsToFullScale;

        doneMap[tileId] = true;

    }

    return;
}


double AcdDigiAlg::edgeEffect(const Event::McPositionHit *hit)  {
    // Purpose and Method:  Compute a modified energy deposit depending upon the
    //   distance the hit is from the center of the tile.
    // Input: Event::McPositionHit pointer
    // Output: The energy associated with this hit

    MsgStream log(msgSvc(), name());
    StatusCode  sc = StatusCode::SUCCESS;

    std::string str;
    std::vector<double> dim;

    idents::VolumeIdentifier volId = hit->volumeID();
    idents::AcdId tileId(volId);
    int iFace = volId[1];

    // retrieve the dimensions of this volume from the GlastDetSvc
    sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
    if ( sc.isFailure() ) {
        log << MSG::WARNING << "Failed to retrieve Shape by Id " 
            << volId.name() << " will not apply edgeEffects for this volume" 
            << endreq;
        return hit->depositedEnergy();
    }
    
    // In local coordinates the box should be centered at (0,0,0)
    const HepPoint3D local_x0 = hit->entryPoint();

    double dX = dim[0];
    double dY = dim[1];
    double dZ = dim[2];

    double dist;
    
    // Calculate distance from edges based on the type of tile and its
    // location in the geometry.
    if(iFace == 0) { // Top Tile
        double dist_x = dX/2. - fabs(local_x0.x());
        double dist_y = dY/2. - fabs(local_x0.y());	                
        dist = (dist_x < dist_y) ? dist_x : dist_y;

        // Multi-volume tiles only occur for top bent tiles
        // which currently are either on the +/- Y sides of the instrument
        if (m_acdId_volCount[tileId] > 1) {
            // First figure out if we're in the main part of the tile or 
            // the bent part
            if (volId[5] == 0) { // main part of tile
                // Retrieve dimensions for the bent part of the tile
                // we'll add the Z-length to the Y edge distance if appropriate
                idents::VolumeIdentifier bentVolId = tileId.volId(true);
                std::string strBent;
                std::vector<double> dimBent;
                sc = m_glastDetSvc->getShapeByID(bentVolId, &strBent, &dimBent);
                if ( sc.isFailure() ) {
                    log << MSG::WARNING << "Failed to retrieve Shape by Id " 
                        << bentVolId.name()
                        << " when applying edgeEffects to bent tiles "
                        << " will ignore bent portion for edgeEffects."
                        << endreq;
                }
                // Row 0 and along edge that is "covered" by bent tile
                if ( (tileId.row() == 0) && (local_x0.y() < 0) )  
                    dist_y += dimBent[2];
                // Row 4 and along edge that is "covered" by bent tile
                else if ( (tileId.row() == 4) && (local_x0.y() > 0) )  
                    dist_y += dimBent[2];

                // Recalculate the minimum distance from the edges
                dist = (dist_x < dist_y) ? dist_x : dist_y;
            

            } else if (volId[5] == 1) { // bent part of tile that is vertical
                // In this case, we add the local Z to the half-length because
                // if we're located in the upper part (z > 0) we need to 
                // measure the distance to the lower edge, if we're in (z < 0)
                // we measure distance to lower edge, and we'll end up 
                // subtracting as we should.
                double dist_z = dZ/2. + local_x0.z();
                // Using Z and X coordinates to check edges, recalc dist
                dist = (dist_x < dist_z) ? dist_x : dist_z;
            } else {
                log << MSG::WARNING << "Unexpected volId for multi-volume "
                    << " ACD tile detector " << volId.name() 
                    << " applying edge effects as if this detector is not bent"
                    << endreq;
            }

        } // end multi-volume handling

    } else if(iFace == 1 || iFace == 3) { // X Side Tile
        double dist_z = dZ/2. - fabs(local_x0.z());
        double dist_y = dX/2. - fabs(local_x0.x()); // these faces are rotated
        dist = (dist_z < dist_y) ? dist_z : dist_y;

    } else if(iFace == 2 || iFace == 4) { // Y Side Tile
        double dist_z = dZ/2. - fabs(local_x0.z());
        double dist_x = dX/2. - fabs(local_x0.x());	                
        dist = (dist_z < dist_x) ? dist_z : dist_x;
    }
    

  
    // Apply edge correction if within m_max_edge_dist (mm) of the edge
    if (dist < m_max_edge_dist) {
        return ( 
           (m_edge_slope*dist + m_edge_intercept) * hit->depositedEnergy() );   
    } else {
        return hit->depositedEnergy();
    }
}

