#define AcdDigi_AcdDigiAlg_CPP 

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

    declareProperty("applyGaussianNoise", m_apply_noise=false);

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

    std::map<idents::VolumeIdentifier, double> energyVolIdMap;
    
    // Create the new AcdDigi collection for the TDS
    Event::AcdDigiCol* digiCol = new Event::AcdDigiCol;
    
    // loop over hits, skip if hit is not in ACD
    // Accumulate the deposited energies, applying edge effects if requested
    for (Event::McPositionHitVector::const_iterator hit = allhits->begin(); hit != allhits->end(); hit++) {
        
        idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)(*hit)->volumeID());
        // Check to see if this is an ACD volume
        if(volId[0] != 1 ) continue; 

        double energy = (m_edge_effect) ? edgeEffect(*hit) : (*hit)->depositedEnergy();

        if (energyVolIdMap.find(volId) != energyVolIdMap.end()) {
            energyVolIdMap[volId] += energy;
        } else {
            energyVolIdMap[volId] = energy;
        }
    }


    // Now loop over the map of ACD VolId and their corresponding energies deposited
    std::map<idents::VolumeIdentifier, double>::const_iterator acdIt;
    for (acdIt = energyVolIdMap.begin(); acdIt != energyVolIdMap.end(); acdIt++) {

        idents::VolumeIdentifier volId = acdIt->first;
        idents::AcdId id(volId);
        
        double energyMevDeposited = acdIt->second;
        log << MSG::DEBUG << "tile volId found: " << volId.name() 
            << ", energy deposited: "<< energyMevDeposited<< " MeV" << endreq;
        
        float mips = util.convertMevToMips(energyMevDeposited);
        
        /// check that we're above low threshold, otherwise skip
        if (mips < m_low_threshold_mips) continue;
        
        // split the signal evenly between the 2 PMTs (A and B)
        float pmtA_mips = mips * 0.5;
        float pmtB_mips = mips - pmtA_mips;
        
        // Number of photoelectrons for each PMT, A and B
        unsigned int pmtA_pe, pmtB_pe;
        util.convertMipsToPhotoElectrons(id, pmtA_mips, pmtA_pe, pmtB_mips, pmtB_pe);
        
        // Apply Poisson fluctuations to the number of pe's for each PMT
        if (m_apply_poisson) {
            pmtA_pe = util.shootPoisson(pmtA_pe);
            pmtB_pe = util.shootPoisson(pmtB_pe);
        }
        
        util.convertPhotoElectronsToMips(id, pmtA_pe, pmtA_mips, pmtB_pe, pmtB_mips);
        
        float pmtA_mipsToFullScale, pmtB_mipsToFullScale;
        
        // If in auto calibrate mode, determine the conversion factor from MIPs
        // to PHA now, at runtime
        if (m_auto_calibrate) {
            util.calcMipsToFullScale(id, pmtA_mips, pmtA_pe, 
                pmtA_mipsToFullScale, pmtB_mips, pmtB_pe, pmtB_mipsToFullScale);
        } else {
            util.applyGains(id, pmtA_mipsToFullScale, pmtB_mipsToFullScale);
        }
        
        float pmtA_observedMips_pha = pmtA_mips;
        float pmtA_observedMips_veto = pmtA_mips;
        float pmtA_observedMips_cno = pmtA_mips;
        
        float pmtB_observedMips_pha = pmtB_mips;
        float pmtB_observedMips_veto = pmtB_mips;
        float pmtB_observedMips_cno = pmtB_mips;
        
        // If Gaussian noise is requested it, apply it here.
        // Noise is applied separately for pha, veto and CNO discriminators
        if (m_apply_noise) {
            pmtA_observedMips_pha += util.shootGaussian(m_noise_std_dev_pha);
            if (pmtA_observedMips_pha < 0.) pmtA_observedMips_pha = 0.;
            pmtA_observedMips_veto += util.shootGaussian(m_noise_std_dev_veto);
            if (pmtA_observedMips_veto < 0.) pmtA_observedMips_veto = 0.;
            pmtA_observedMips_cno += util.shootGaussian(m_noise_std_dev_cno);
            if (pmtA_observedMips_cno < 0.) pmtA_observedMips_cno = 0.;
            
            pmtB_observedMips_pha += util.shootGaussian(m_noise_std_dev_pha);
            if (pmtB_observedMips_pha < 0.) pmtB_observedMips_pha = 0.;
            pmtB_observedMips_veto += util.shootGaussian(m_noise_std_dev_veto);
            if (pmtB_observedMips_veto < 0.) pmtB_observedMips_veto = 0.;
            pmtB_observedMips_cno += util.shootGaussian(m_noise_std_dev_cno);
            if (pmtB_observedMips_cno < 0.) pmtB_observedMips_cno = 0.;
        }
        
        // Now convert MIPs into PHA values for each PMT
        unsigned short pmtA_pha = util.convertMipsToPha(pmtA_observedMips_pha, pmtA_mipsToFullScale);
        unsigned short pmtB_pha = util.convertMipsToPha(pmtB_observedMips_pha, pmtB_mipsToFullScale);
        
        // Initialize discriminators
        bool lowArr[2] = { true, true };
        bool vetoArr[2] = { false, false };
        bool highArr[2] = { false, false };
        
        // Set the discriminators if above threshold
        if (pmtA_observedMips_veto > m_veto_threshold_mips) vetoArr[0] = true;
        if (pmtB_observedMips_veto > m_veto_threshold_mips) vetoArr[1] = true;
        
        if (pmtA_observedMips_cno > m_high_threshold_mips) highArr[0] = true;
        if (pmtB_observedMips_cno > m_high_threshold_mips) highArr[1] = true;
        
        unsigned short phaArr[2] = { pmtA_pha, pmtB_pha };
        
        log << MSG::DEBUG << "AcdId: " << id.id() <<endreq;
        
        // Add this AcdDigi to the TDS collection
        digiCol->push_back(
            new AcdDigi(
            id, volId,
            energyMevDeposited, phaArr, 
            vetoArr, lowArr, highArr) );   
        
    } // end loop over MC hits
    
    // Put the AcdDigi collection on the TDS
    return eventSvc()->registerObject(EventModel::Digi::AcdDigiCol, digiCol);    
}


StatusCode AcdDigiAlg::finalize() {
    return StatusCode::SUCCESS;
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
    
    // Apply edge correction if within 20 mm of the edge
    if (dist < 20.0) {
        return ( (0.1*dist + 0.8) * hit->depositedEnergy() );   
    } else {
        return hit->depositedEnergy();
    }
}