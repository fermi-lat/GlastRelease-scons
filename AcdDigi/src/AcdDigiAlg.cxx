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
#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"

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
    
    m_glastDetSvc->accept(m_tiles);
    log << MSG::INFO << "will add noise to "<< m_tiles.size() << " ACD tiles, ids from "
        << m_tiles.front().name() << " to " << m_tiles.back().name() << endreq;

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
        return sc;
    }

    //Take care of insuring that data area has been created
    DataObject* pNode = 0;
    sc = eventSvc()->retrieveObject( EventModel::Digi::Event , pNode);
    
    if (sc.isFailure()) {
        sc = eventSvc()->registerObject(EventModel::Digi::Event ,new DataObject);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register " << EventModel::Digi::Event << endreq;
            return sc;
        }
    }

   // std::map<idents::VolumeIdentifier, double> energyVolIdMap;
    std::map<idents::AcdId, double> energyIdMap;
    
    // Create the new AcdDigi collection for the TDS
    Event::AcdDigiCol* digiCol = new Event::AcdDigiCol;
    
    // loop over hits, skip if hit is not in ACD
    // Accumulate the deposited energies, applying edge effects if requested
    for (Event::McPositionHitVector::const_iterator hit = allhits->begin(); hit != allhits->end(); hit++) {
        
        idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)(*hit)->volumeID());
        // Check to see if this is an ACD volume
        if(volId[0] != 1 ) continue; 

        idents::AcdId id(volId);
        double energy;
        // No edge effects for ribbons, tiles only
        energy = (m_edge_effect && id.tile()) ? edgeEffect(*hit) : (*hit)->depositedEnergy();
        if (energyIdMap.find(id) != energyIdMap.end()) {
            energyIdMap[id] += energy;
        } else {
            energyIdMap[id] = energy;
        }
        
    }

    if (m_apply_noise) addNoise();

    // Now loop over the map of AcdId and their corresponding energies deposited
    std::map<idents::AcdId, double>::const_iterator acdIt;
    for (acdIt = energyIdMap.begin(); acdIt != energyIdMap.end(); acdIt++) {

        //idents::VolumeIdentifier volId = acdIt->first;
        idents::AcdId id = acdIt->first; //(volId);

        double energyMevDeposited = acdIt->second;
        m_energyDepMap[id] = energyMevDeposited;

        log << MSG::DEBUG << "tile id found: " << id.id() 
            << ", energy deposited: "<< energyMevDeposited<< " MeV" << endreq;
        
        double mips = util.convertMevToMips(energyMevDeposited);
                
        // split the signal evenly between the 2 PMTs (A and B)
        double pmtA_mips = mips * 0.5;
        double pmtB_mips = mips - pmtA_mips;
        
        // Number of photoelectrons for each PMT, A and B
        unsigned int pmtA_pe, pmtB_pe;
        util.convertMipsToPhotoElectrons(id, pmtA_mips, pmtA_pe, pmtB_mips, pmtB_pe);
        
        // Apply Poisson fluctuations to the number of pe's for each PMT
        if (m_apply_poisson) {
            pmtA_pe = util.shootPoisson(pmtA_pe);
            pmtB_pe = util.shootPoisson(pmtB_pe);
        }
        
        util.convertPhotoElectronsToMips(id, pmtA_pe, pmtA_mips, pmtB_pe, pmtB_mips);
        
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
        m_pmtA_cnoMipsMap[id] = pmtB_observedMips_cno;

        // check that we're above low threshold, otherwise skip
        //if ((pmtA_observedMips_pha + pmtB_observedMips_pha) < m_low_threshold_mips) continue;


        // Now convert MIPs into PHA values for each PMT
        //unsigned short pmtA_pha = util.convertMipsToPha(pmtA_observedMips_pha, pmtA_mipsToFullScale);
        //unsigned short pmtB_pha = util.convertMipsToPha(pmtB_observedMips_pha, pmtB_mipsToFullScale);
        
        // Initialize discriminators
        //bool lowArr[2] = { true, true };
        //bool vetoArr[2] = { false, false };
        //bool highArr[2] = { false, false };
        
        // Set the discriminators if above threshold
        //if (pmtA_observedMips_veto > m_veto_threshold_mips) vetoArr[0] = true;
        //if (pmtB_observedMips_veto > m_veto_threshold_mips) vetoArr[1] = true;
        
      //  if (pmtA_observedMips_cno > m_high_threshold_mips) highArr[0] = true;
        //if (pmtB_observedMips_cno > m_high_threshold_mips) highArr[1] = true;

        //unsigned short phaArr[2] = { pmtA_pha, pmtB_pha };
        
        log << MSG::DEBUG << "AcdId: " << id.id() <<endreq;
        

        // Add this AcdDigi to the TDS collection
        //digiCol->push_back(
        //    new AcdDigi(
        //    id, volId,
        //    energyMevDeposited, phaArr, 
        //    vetoArr, lowArr, highArr) );   
        
    } // end loop over MC hits


    std::map<idents::AcdId, bool> doneMap;

    // Now fill the TDS with AcdDigis
    for(AcdTileList::const_iterator it=m_tiles.begin(); it!=m_tiles.end(); ++it){
        idents::VolumeIdentifier volId = *it;
        idents::AcdId tileId(volId);

        // Check to see if we have already processed this AcdId - necessary since we also have
        // ribbons in the mix.
        if (doneMap.find(tileId) != doneMap.end()) continue;
        
        // First check to see if this tile Id has a PHA value associated with PMT A
        // If not, we do not need add an AcdDigi in the TDS for this tile
        if (m_pmtA_phaMipsMap.find(tileId) == m_pmtA_phaMipsMap.end()) continue;

        // Next check that the PHA values from both PMTs combined results is a value above low 
        // threshold
        if ((m_pmtA_phaMipsMap[tileId] + m_pmtB_phaMipsMap[tileId]) < m_low_threshold_mips) continue;

        // Initialize discriminators
        bool lowArr[2] = { true, true };
        bool vetoArr[2] = { false, false };
        bool highArr[2] = { false, false };
        
        // Set the discriminators if above threshold
        if (m_pmtA_vetoMipsMap[tileId] > m_veto_threshold_mips) vetoArr[0] = true;
        if (m_pmtB_vetoMipsMap[tileId] > m_veto_threshold_mips) vetoArr[1] = true;
        
        if (m_pmtA_cnoMipsMap[tileId] > m_high_threshold_mips) highArr[0] = true;
        if (m_pmtB_cnoMipsMap[tileId] > m_high_threshold_mips) highArr[1] = true;

        // Now convert MIPs into PHA values for each PMT
        unsigned short pmtA_pha = util.convertMipsToPha(m_pmtA_phaMipsMap[tileId], m_pmtA_toFullScaleMap[tileId]);
        unsigned short pmtB_pha = util.convertMipsToPha(m_pmtB_phaMipsMap[tileId], m_pmtB_toFullScaleMap[tileId]);

        unsigned short phaArr[2] = { pmtA_pha, pmtB_pha };

        doneMap[tileId] = true;

        digiCol->push_back(
            new AcdDigi(
            tileId, volId,
            m_energyDepMap[tileId], phaArr, 
            vetoArr, lowArr, highArr) );   
    }

    clear();

    // Put the AcdDigi collection on the TDS
    return eventSvc()->registerObject(EventModel::Digi::AcdDigiCol, digiCol);    
}


StatusCode AcdDigiAlg::finalize() {
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

}


void AcdDigiAlg::getParameters() {    
    // Purpose and Method:  Read in the parameters from the XML file using IFile.
    //   Perform some spot checking to see if the parameters exist in the input 
    //   XML file - if not provide an informational message.

    MsgStream log(msgSvc(), name());
 
    xml::IFile xmlFilePtr(m_xmlFile.c_str());
    
        // Perform some spot checking to see if our XML file contains our constants
    if (!xmlFilePtr.contains("thresholds", "low_threshold_mips")) {
        log << MSG::INFO << "XML file " << m_xmlFile
            << " does not contain low_threshold_mips using default" << endreq;
    }

    m_low_threshold_mips = xmlFilePtr.getDouble("thresholds", "low_threshold_mips", 0.1);
    m_veto_threshold_mips = xmlFilePtr.getDouble("thresholds", "veto_threshold_mips", 0.3);
    m_high_threshold_mips = xmlFilePtr.getDouble("thresholds", "high_threshold_mips", 10.5);
    
    if (!xmlFilePtr.contains("global_constants", "mean_pe_per_mip")) {
        log << MSG::INFO << "XML file " << m_xmlFile
            << "does not contain mean_pe_per_mip using default" << endreq;
    }

    m_mean_pe_per_mip = xmlFilePtr.getInt("global_constants", "mean_pe_per_mip", 18);
    
    m_noise_std_dev_pha = xmlFilePtr.getDouble("global_constants", "noise_std_dev_pha", 0.02);
    m_noise_std_dev_veto = xmlFilePtr.getDouble("global_constants", "noise_std_dev_veto", 0.02);
    m_noise_std_dev_veto = xmlFilePtr.getDouble("global_constants", "noise_std_dev_cno", 0.02);
    
    m_full_scale = xmlFilePtr.getInt("global_constants", "full_scale", 4095);
    
    m_mips_full_scale = xmlFilePtr.getDouble("global_constants", "mips_full_scale", 20.0);
    
    m_mev_per_mip = xmlFilePtr.getDouble("global_constants", "mev_per_mip", 1.9);

    m_max_edge_dist = xmlFilePtr.getDouble("edge_effects", "max_edge_dist", 20.0);

    m_edge_slope = xmlFilePtr.getDouble("edge_effects", "edge_slope", 0.01);

    m_edge_intercept = xmlFilePtr.getDouble("edge_effects", "edge_intercept", 0.8);

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
        // due to presence of ribbons
        if (doneMap.find(tileId) != doneMap.end()) continue;
        
        
        m_pmtA_phaMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_pha);
        m_pmtA_vetoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_veto);
        m_pmtA_cnoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_cno);
        
        m_pmtB_phaMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_pha);
        m_pmtB_vetoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_veto);
        m_pmtB_cnoMipsMap[tileId] = util.shootGaussian(m_noise_std_dev_cno);

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
    int iFace = volId[1];

    // retrieve the dimensions of this volume from the GlastDetSvc
    sc = m_glastDetSvc->getShapeByID(volId, &str, &dim);
    if ( sc.isFailure() ) {
        log << MSG::DEBUG << "Failed to retrieve Shape by Id" << endreq;
        return sc;
    }
    
    // In local coordinates the box should be centered at (0,0,0)
    const HepPoint3D local_x0 = hit->entryPoint();

    double dX = dim[0];
    double dY = dim[1];
    double dZ = dim[2];
    
    double dist;
    
    if(iFace == 0) { // Top Tile
        double dist_x = dX/2. - fabs(local_x0.x());
        double dist_y = dY/2. - fabs(local_x0.y());	                
        dist = (dist_x < dist_y) ? dist_x : dist_y;
    }
    else if(iFace == 1 || iFace == 3) { // X Side Tile
        double dist_z = dZ/2. - fabs(local_x0.z());
        double dist_y = dY/2. - fabs(local_x0.y());	                
        dist = (dist_z < dist_y) ? dist_z : dist_y;
    }
    else if(iFace == 2 || iFace == 4) { // Y Side Tile
        double dist_z = dZ/2. - fabs(local_x0.z());
        double dist_x = dY/2. - fabs(local_x0.x());	                
        dist = (dist_z < dist_x) ? dist_z : dist_x;
    }
    
    // Apply edge correction if within m_max_edge_dist (mm) of the edge
    if (dist < m_max_edge_dist) {
        return ( (m_edge_slope*dist + m_edge_intercept) * hit->depositedEnergy() );   
    } else {
        return hit->depositedEnergy();
    }
}