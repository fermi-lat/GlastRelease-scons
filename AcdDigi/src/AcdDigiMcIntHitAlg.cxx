#define AcdDigi_AcdDigiMcIntHitAlg_CPP 

#include "AcdDigiMcIntHitAlg.h"

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/TopLevel/Event.h"
#include "Event/Digi/AcdDigi.h"
#include "Event/MonteCarlo/McIntegratingHit.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGauss.h"

// for min and floor functions
#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

// Define the factory for this algorithm
static const AlgFactory<AcdDigiMcIntHitAlg>  Factory;
const IAlgFactory& AcdDigiMcIntHitAlgFactory = Factory;


// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor.
AcdDigiMcIntHitAlg::AcdDigiMcIntHitAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    
    declareProperty ("xmlFile", m_xmlFile="$(ACDDIGIROOT)/xml/acdDigi.xml");
}


StatusCode AcdDigiMcIntHitAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    log << MSG::INFO << "initialize" << endreq;
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
    
    getParameters();

    util.getParameters(m_xmlFile);
    
    return StatusCode::SUCCESS;
}


StatusCode AcdDigiMcIntHitAlg::execute() {
    // Purpose and Method:  Using the McIntegratingHits that hit the ACD tiles
    //   construct the AcdDigi collection.
    // TDS Input:  EventModel::MC::McIntegratingHitCol
    // TDS Output: EventModel::Digi::AcdDigiCol

    StatusCode  sc = StatusCode::SUCCESS;
    MsgStream   log( msgSvc(), name() );
    
    using namespace Event;
    
    SmartDataPtr<Event::McIntegratingHitCol> allhits(eventSvc(),EventModel::MC::McIntegratingHitCol );
    
    if (!allhits) {
        log << MSG::INFO << "No McIntegratingHits for this event" << endreq;
        return sc;
    }

    //Take care of insuring that data area has been created
    DataObject* pNode = 0;
    sc = eventSvc()->retrieveObject( EventModel::Digi::Event , pNode);
    
    if (sc.isFailure()) {
        sc = eventSvc()->registerObject(EventModel::Digi::Event, new DataObject);
        if( sc.isFailure() ) {
            log << MSG::ERROR << "could not register " << EventModel::Digi::Event << endreq;
            return sc;
        }
    }
    
    // Create the new AcdDigi collection for the TDS
    Event::AcdDigiCol* digiCol = new Event::AcdDigiCol;
    
    // loop over hits, skip if hit is not in ACD
    // if below low threshold, skip
    for (Event::McIntegratingHitVector::const_iterator it = allhits->begin(); it != allhits->end(); it++) {
        
        idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)(*it)->volumeID());
        // Check to see if this is an ACD volume
        if(volId[0] != 1 ) continue; 
        
        // unpack from volume id
        idents::AcdId id(volId);
        
        double energyMevDeposited = (*it)->totalEnergy();
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
        
        // Apply Gaussian noise if requested
        // where the noise is applied separately for pha, veto and CNO discrim
        if (m_apply_noise) {
            pmtA_observedMips_pha += util.shootGaussian(m_noise_std_dev_pha);
            pmtA_observedMips_veto += util.shootGaussian(m_noise_std_dev_veto);
            pmtA_observedMips_cno += util.shootGaussian(m_noise_std_dev_cno);
            
            pmtB_observedMips_pha += util.shootGaussian(m_noise_std_dev_pha);
            pmtB_observedMips_veto += util.shootGaussian(m_noise_std_dev_veto);
            pmtB_observedMips_cno += util.shootGaussian(m_noise_std_dev_cno);
        }
        
        // Now convert MIPs into PHA values for each PMT
        unsigned short pmtA_pha = util.convertMipsToPha(pmtA_observedMips_pha, pmtA_mipsToFullScale);
        unsigned short pmtB_pha = util.convertMipsToPha(pmtB_observedMips_pha, pmtB_mipsToFullScale);
        
        // initialize the discriminators
        bool lowArr[2] = { true, true };
        bool vetoArr[2] = { false, false };
        bool highArr[2] = { false, false };
        
        // Set the discriminators if above threshold
        if (pmtA_observedMips_veto > m_veto_threshold_mips) vetoArr[0] = true;
        if (pmtB_observedMips_veto > m_veto_threshold_mips) vetoArr[1] = true;
        
        if (pmtA_observedMips_cno > m_high_threshold_mips) highArr[0] = true;
        if (pmtB_observedMips_cno > m_high_threshold_mips) highArr[1] = true;
        
        unsigned short phaArr[2] = { pmtA_pha, pmtB_pha };
        
        log << MSG::DEBUG << "AcdId: " << id.id() << endreq;
        
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


StatusCode AcdDigiMcIntHitAlg::finalize() {    
    return StatusCode::SUCCESS;
}


void AcdDigiMcIntHitAlg::getParameters() {    
    // Purpose and Method:  Read in the parameters from the XML file using IFile
    //   Performs some checks to see if the XML file contains some of the expected
    //   parameters, if not, an information message is provided as output.
    
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

    m_auto_calibrate = xmlFilePtr.getBool("processing", "auto_calibrate", true);
    
    m_apply_poisson = xmlFilePtr.getBool("processing", "apply_poisson", true);

    
    return;
}


