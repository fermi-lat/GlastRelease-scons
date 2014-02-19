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
#include "OverlayEvent/GemOverlay.h"
#include "OverlayEvent/OverlayEventModel.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "AcdUtil/IAcdGeometrySvc.h"
#include "AcdUtil/IAcdCalibSvc.h"
#include "AcdUtil/IAcdFailureModeSvc.h"

#include "CLHEP/Random/RandExponential.h"

// Define the factory for this algorithm
//static const AlgFactory<AcdDigiAlg>  Factory;
//const IAlgFactory& AcdDigiAlgFactory = Factory;
DECLARE_ALGORITHM_FACTORY(AcdDigiAlg);


// Algorithm parameters which can be set at run time must be declared.
// This should be done in the constructor.
AcdDigiAlg::AcdDigiAlg(const std::string& name, ISvcLocator* pSvcLocator) :
Algorithm(name, pSvcLocator) {
    
    declareProperty("xmlFile", m_xmlFileName="$(ACDDIGIXMLPATH)/acdDigi.xml");
    declareProperty("AcdSimCalibSvc",    m_calibSvcName = "AcdSimCalibSvc");
    declareProperty("applyPoisson", m_apply_poisson=true);
    declareProperty("applyGaussianNoise", m_apply_noise=true);
    declareProperty("applyCoherentNoise", m_apply_coherent_noise=true);
    declareProperty("coherentNoiseInOverlay", m_coherent_noise_in_overlay=true);
    declareProperty("edgeEffect", m_edge_effect=true);
    declareProperty("lightYeildRatio", m_lightYeildRatio=1.0);
    //Change made by D. Green to incorporate zero suppression as a JO
    declareProperty("phaZeroThreshold", m_phaZeroThreshold=25.0);
}


StatusCode AcdDigiAlg::initialize() {
    
    MsgStream log(msgSvc(), name());
    StatusCode  sc = StatusCode::SUCCESS;
    
    AcdUtil::IAcdCalibSvc* calibSvc(0);
    sc = service(m_calibSvcName, calibSvc, true);
    if (sc.isSuccess() ) 
    {
        sc = calibSvc->queryInterface(IID_IAcdCalibSvc, (void**)&calibSvc);
    }

    if ( !sc.isSuccess() ) 
    {
        log << MSG::WARNING << "Could not get CalibDataSvc " << m_calibSvcName << endreq;
        return sc;
    } else {
        log << MSG::INFO << "Got CalibDataSvc " << m_calibSvcName << endreq;
    }
      
    sc = toolSvc()->retrieveTool("AcdDigiRandom", m_randTool);
    if (sc.isFailure() )
      log << MSG::WARNING << "Unable to register AcdDigiRandom" << endreq;

    IGlastDetSvc* glastDetSvc(0);
    sc = service("GlastDetSvc", glastDetSvc, true);
    if (sc.isSuccess() ) 
    {
        sc = glastDetSvc->queryInterface(IID_IGlastDetSvc, (void**)&glastDetSvc);
    }
   
    if ( !sc.isSuccess() ) 
    {
        log << MSG::WARNING << "AcdDigiAlg failed to get the GlastDetSvc" << endreq;
        return sc;
    } else {
        log << MSG::INFO << "Got GlastDetSvc" << endreq;
    }

    // get the list of layers, to be used to add noise to otherwise empty layers
    m_tiles.setPrefix(glastDetSvc->getIDPrefix());
    
    // Find all the ACD detectors in our geometry
    glastDetSvc->accept(m_tiles);
    if (m_tiles.size() > 0) 
        log << MSG::INFO << "Located  " << m_tiles.size() 
            << " ACD volumes, ids from " << m_tiles.front().name() 
            << " to " << m_tiles.back().name() << endreq;

    IAcdGeometrySvc* acdGeomSvc(0);
    sc = service("AcdGeometrySvc", acdGeomSvc, true);
    if (sc.isSuccess() ) 
    {
        sc = acdGeomSvc->queryInterface(IID_IAcdGeometrySvc, (void**)&acdGeomSvc);
    }
   
    if ( !sc.isSuccess() ) 
    {
        log << MSG::WARNING << "AcdDigiAlg failed to get the AcdGeometrySvc" << endreq;
        return sc;
    } else {
        log << MSG::INFO << "Got AcdGeometrySvc" << endreq;
    }

    sc =service("AcdFailureModeSvc", m_acdFailureSvc);
    if (sc.isSuccess() ) 
        log << MSG::INFO << "AcdFailureModeSvc loaded" << endreq;
    else {
        log << MSG::INFO << "AcdFailureModeSvc not found, skipping" << endreq;
        m_acdFailureSvc = 0;
    }

    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();

    
    // read in the parameters from our input XML file
    sc = m_util.initialize(*calibSvc,*acdGeomSvc,m_xmlFileName);

    if ( !sc.isSuccess() ) 
    {
        log << MSG::WARNING << "Failed to intiliaze AcdDigiUtil" << endreq;
        return sc;
    } else {
        log << MSG::INFO << "Intiliazed AcdDigiUtil" << endreq;
    }


    if ( m_apply_coherent_noise ) {
      // Convention for multiple input overlay files is that there will be separate OverlayDataSvc's with 
      // names appended by "_xx", for example OverlayDataSvc_1 for the second input file. 
      // In order to ensure the data read in goes into a unique section of the TDS we need to modify the 
      // base root path, which we do by examining the name of the service
      std::string dataSvcName = "OverlayDataSvc";
      int         subPos      = name().rfind("_");
      std::string nameEnding  = subPos > 0 ? name().substr(subPos, name().length() - subPos) : "";
      
      if (nameEnding != "") dataSvcName += nameEnding;
      
      IService* dataSvc = 0;
      sc = service(dataSvcName, dataSvc);
      if (sc.isFailure() ) {
          log << MSG::ERROR << "  can't get OverlayDataSvc " << endreq;
	  return sc;
      }

      // Caste back to the "correct" pointer
      m_dataSvc = dynamic_cast<DataSvc*>(dataSvc);      
    } else {
      m_dataSvc = 0;
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

    clear();
    
    SmartDataPtr<Event::McPositionHitCol> allhits(eventSvc(),EventModel::MC::McPositionHitCol );

    if (!allhits) 
    {
        log << MSG::INFO << "No McPositionHits were found in the TDS" << endreq;
        if (!m_apply_noise) return sc;
        // Otherwise, Do not return - there may be noise hits
    }

    //Take care of insuring that data area has been created
    DataObject* pNode = 0;
    sc = eventSvc()->retrieveObject( EventModel::Digi::Event , pNode);
    
    if (sc.isFailure()) 
    {
        sc = eventSvc()->registerObject(EventModel::Digi::Event, 
                                        new Event::DigiEvent);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "could not register " 
                << EventModel::Digi::Event << endreq;
            return sc;
        }
    }

    if ( m_apply_coherent_noise ) {
      SmartDataPtr<Event::GemOverlay> gemOverlay(m_dataSvc, m_dataSvc->rootName() + OverlayEventModel::Overlay::GemOverlay);
      if( !gemOverlay) {
	log << MSG::ERROR << "could not find gemOverlay object." << endreq;
	return sc;
      }


      m_gemDeltaEventTime = gemOverlay->getDeltaEventTime();
    } else {
      m_gemDeltaEventTime = 5000;
    }


    // Create the new AcdDigi collection for the TDS
    Event::AcdDigiCol* digiCol = new Event::AcdDigiCol;
    
    // loop over hits, skip if hit is not in ACD
    // Accumulate the deposited energies, applying edge effects if requested
    if (allhits) 
    {
        sc = fillEnergyAndPeMaps(*allhits,m_energyDepMap, m_peMap, m_statusMap);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Failed to fill energy deposit and p.e. maps" << endreq;
            return sc;
        }  
      
        sc = convertPeToMips(m_peMap, m_statusMap, m_mipsMap);
        if( sc.isFailure() ) 
        {
            log << MSG::ERROR << "Failed to convert energy to mips" << endreq;
            return sc;
        }  
    }
    
    sc = makeDigis(m_mipsMap, m_statusMap, *digiCol);

    clear();

    // Put the AcdDigi collection on the TDS
    return eventSvc()->registerObject(EventModel::Digi::AcdDigiCol, digiCol);    
}


StatusCode AcdDigiAlg::fillEnergyAndPeMaps( const Event::McPositionHitCol& mcHits,
                        std::map<idents::AcdId, double>& energyIdMap,
                        std::map<idents::AcdId, std::pair<double, double> >& peMap,
                        std::map<idents::AcdId, unsigned int>& statusMap) 
{  
    MsgStream log( msgSvc(), name() );
      
    for (Event::McPositionHitVector::const_iterator hit = mcHits.begin();
         hit != mcHits.end(); hit++) 
    {
        idents::VolumeIdentifier volId = ((idents::VolumeIdentifier)(*hit)->volumeID());

        // Check to see if this is an ACD volume
        if(volId[0] != 1 ) continue; 
    
        idents::AcdId id(volId);
        // No edge effects for ribbons, tiles only

        // energy from GEANT
        double energy = (*hit)->depositedEnergy();
        // mean # of p.e. at each PMT
        double pe_pmtA_mean(0.);
        double pe_pmtB_mean(0.);

        if ( id.tile() ) 
        {
            // For tiles we check for edge effects
            StatusCode sc = m_util.photoElectronsFromEnergy_tile((*hit),
                                                                  m_edge_effect,
                                                                  log,
                                                                  pe_pmtA_mean,
                                                                  pe_pmtB_mean);
            if ( sc.isFailure() ) 
            {
                log << MSG::ERROR << "Couldn't get p.e. for tile " << volId.name() << ' ' << energy << endreq;
                return sc;
            }
        } else if ( id.ribbon() ) {
            // For ribbons don't bother with edge effects, but do care about attenuation effects
            StatusCode sc = m_util.photoElectronsFromEnergy_ribbon((*hit),log,pe_pmtA_mean,pe_pmtB_mean);  
            if ( sc.isFailure() ) 
            {
                log << MSG::ERROR << "Couldn't get p.e. for ribbon " << volId.name() << ' ' << energy << endreq;
                return sc;
            }
        } else {
            log << MSG::ERROR << "Deposited energy in ACD from neither tile nor ribbon " << volId.name() << ' ' << energy << endreq;
            return StatusCode::FAILURE;
        }
    
        // Adjust the light yield
        pe_pmtA_mean *= m_lightYeildRatio;
        pe_pmtB_mean *= m_lightYeildRatio;
    
        // Add the deposited energy to the map
        if (energyIdMap.find(id) != energyIdMap.end()) 
        {
            energyIdMap[id] += energy;
        } else {
            energyIdMap[id] = energy;
        }   
    
        // Add the p.e. expectation values to the map
        if ( peMap.find(id) != peMap.end() ) 
        {
            peMap[id].first  += pe_pmtA_mean;
            peMap[id].second += pe_pmtB_mean;
        } else {
            peMap[id] = std::make_pair(pe_pmtA_mean,pe_pmtB_mean);
        }

        // Set the status bits
        unsigned int statusWord = 0;
        
        // If overlay bit not set (and in this code) then must be MC
        if (!((*hit)->getPackedFlags() & Event::McPositionHit::overlayHit)) 
             statusWord = Event::AcdDigi::DIGI_MC;

        // If overlay bit is set, copy entire packed word (note that top bits have same def digi<->McPositionHit
        else statusWord = (*hit)->getPackedFlags();

        if (statusMap.find(id) != statusMap.end())
        {
            statusMap[id] |= statusWord;
        } else {
            statusMap[id]  = statusWord;
        }

    } // finished loop over MC hits

    return StatusCode::SUCCESS;
}


StatusCode AcdDigiAlg::convertPeToMips( const std::map<idents::AcdId, std::pair<double, double> >& peMap,
                                        const std::map<idents::AcdId, unsigned int>&               statusMap,
                                              std::map<idents::AcdId, std::pair<double, double> >& mipsMap) 
{
    MsgStream log( msgSvc(), name() );

    for ( std::map<idents::AcdId, std::pair<double, double> >::const_iterator itr = peMap.begin(); 
          itr != peMap.end(); itr++ ) 
    {
        // get the expectation values
        double pe_pmtA_mean = itr->second.first;
        double pe_pmtB_mean = itr->second.second;

        // initialize correct responses
        double pe_pmtA = pe_pmtA_mean;
        double pe_pmtB = pe_pmtB_mean;

        // Retrieve the status bits 
        unsigned int status = 0;
        
        if (statusMap.find(itr->first) != statusMap.end()) status = statusMap.find(itr->first)->second;

        // If not a pure overlay hit then do the dynode simulation
        if ((status & 0xF0000000) != (unsigned int)Event::AcdDigi::DIGI_OVERLAY)
        {
            // Throw the Poisson stats on the expected # of PE.
            pe_pmtA = AcdDigiUtil::simulateDynodeChain(pe_pmtA_mean);
            pe_pmtB = AcdDigiUtil::simulateDynodeChain(pe_pmtB_mean);

            // Adjust the light yield
            pe_pmtA /= m_lightYeildRatio;
            pe_pmtB /= m_lightYeildRatio;
        }

        double mipA(0.); 
        double mipB(0.);
        // Convert from "measured" PE back to mips

        StatusCode sc  = m_util.mipEquivalentLightYeild(itr->first,pe_pmtA,pe_pmtB,log,mipA,mipB);
        if ( sc.isFailure() ) return sc;

        mipsMap[itr->first] = std::make_pair(mipA,mipB);
    
        //uncomment to make it easier to test ribbons
        //if ( itr->first.tile() ) continue;
        log << MSG::DEBUG << "PMT: " << itr->first.id() 
        << "  E: " << m_energyDepMap[itr->first]
        << "  pe_mean: " << pe_pmtA_mean << ',' << pe_pmtB_mean
        << "  pe: " << pe_pmtA << ',' << pe_pmtB
        << "  mips: " << mipA << ',' << mipB << endreq;  
    }

    return StatusCode::SUCCESS;
}

StatusCode AcdDigiAlg::makeDigis(const std::map<idents::AcdId, std::pair<double,double> >& mipsMap,
                                 const std::map<idents::AcdId, unsigned int>& statusMap,
                                 Event::AcdDigiCol& digiCol) 
{
    MsgStream log( msgSvc(), name() );

    std::set<idents::AcdId> doneMap;

    // Now fill the TDS with AcdDigis
    // Loop over all tiles in the geometry
    for(AcdTileList::const_iterator it=m_tiles.begin(); it!=m_tiles.end(); ++it)
    {
        idents::VolumeIdentifier volId = *it;
        idents::AcdId acdId(volId);
 
        // Check to see if AcdFailureModeSvc is loaded
        if (m_acdFailureSvc) {
            // Check if this AcdId is in the list of killed detectors
            if (m_acdFailureSvc->matchAcdId(acdId))
                continue; // Skip this ACD detector
        }

        // Check to see if we have already processed this AcdId - 
        // necessary since we also have ribbons in the mix and curved tiles
        // have muliple pieces (volumes)
        if (doneMap.find(acdId) != doneMap.end()) continue;
        doneMap.insert(acdId);

        double mipsPmt[2] = {0.,0.};

        std::map<idents::AcdId, std::pair<double,double> >::const_iterator itrFind = mipsMap.find(acdId);
        if ( itrFind != mipsMap.end() ) 
        {
            mipsPmt[0] = itrFind->second.first;
            mipsPmt[1] = itrFind->second.second;;      
        } 
	
        unsigned short phaArr[2] = { 0, 0 };
        Event::AcdDigi::Range rangeArr[2] = { Event::AcdDigi::LOW, Event::AcdDigi::LOW };
        bool phaThreshArr[2] = { false, false };
        bool vetoArr[2] = { false, false };
        bool highArr[2] = { false, false };
        bool makeDigi(false);

        // Recover the status word for this digi
        unsigned int statusWord = Event::AcdDigi::DIGI_MC;
        if (statusMap.find(acdId) != statusMap.end()) statusWord = statusMap.find(acdId)->second;

        StatusCode sc = m_util.phaCounts(acdId,mipsPmt,m_apply_noise,log,rangeArr,phaArr);
        if ( sc.isFailure() ) return sc;
    
        if ( m_apply_coherent_noise ) 
        {  
       	    // Make sure we aren't double counting coherent noise from overlays 
	    bool noiseFromOverlay = m_coherent_noise_in_overlay && ( statusWord & (Event::AcdDigi::DIGI_OVERLAY) != 0 );
   	    if  ( ! noiseFromOverlay ) 
	    {     
                sc = m_util.applyCoherentNoiseToPha(acdId,(unsigned)m_gemDeltaEventTime,log,phaArr);
		if ( sc.isFailure() ) return sc;  
	    }
        }

        //Change made by D. Green to incorporate zero suppression as a JO
        sc = m_util.checkThresholds(acdId,mipsPmt,phaArr,rangeArr,m_apply_noise,log,makeDigi,phaThreshArr,vetoArr,highArr,m_phaZeroThreshold);
        if ( sc.isFailure() ) return sc;


        // Now do special handling for pure Overlay events
        // Pure overlay?
        if ((statusWord & (Event::AcdDigi::DIGI_OVERLAY | Event::AcdDigi::DIGI_MC)) == (unsigned int) Event::AcdDigi::DIGI_OVERLAY)
        {
            // Use the Overlay event's bits to attempt to reproduce ghost event behavior
            // NOTE: bit definions in OverlayEvent::AcdOverlay - need to define in better location... (upgrade!)
            phaThreshArr[Event::AcdDigi::A] = (statusWord & 0x00000001) != 0;
            vetoArr[Event::AcdDigi::A]      = (statusWord & 0x00000002) != 0;
            highArr[Event::AcdDigi::A]      = (statusWord & 0x00000008) != 0;
            phaThreshArr[Event::AcdDigi::B] = (statusWord & 0x00010000) != 0;
            vetoArr[Event::AcdDigi::B]      = (statusWord & 0x00020000) != 0;
            highArr[Event::AcdDigi::B]      = (statusWord & 0x00080000) != 0;
        }
        else
        {   
            // Ok, kill PHA below ZS threshold
            phaArr[0] = phaThreshArr[0] ? phaArr[0] : 0;
            phaArr[1] = phaThreshArr[1] ? phaArr[1] : 0;
        }

        bool gemBit = vetoArr[0] || vetoArr[1];
        bool ninjaBit = gemBit && ( ! ( phaThreshArr[0] || phaThreshArr[1] ) );               

        if ( makeDigi ) 
        {
            // Create the digi 
            Event::AcdDigi* aDigi = new Event::AcdDigi(acdId, 
                                                       volId,
                                                       m_energyDepMap[acdId], 
                                                       phaArr, 
                                                       vetoArr, 
                                                       phaThreshArr, 
                                                       highArr);

            // Fill in the rest of the stuff
            aDigi->setRanges(rangeArr);
            aDigi->initGem(ninjaBit,gemBit);
            aDigi->setStatus(statusWord);
            digiCol.push_back( aDigi );

            //uncomment to make it easier to test ribbons
            //if ( acdId.tile() ) continue;
            log << MSG::DEBUG << "PMT: " << acdId.id()
                << "  E: " << m_energyDepMap[acdId]
                << "  mips: " << mipsPmt[0] << ',' <<  mipsPmt[1]
                << "  PHA: " << phaArr[0] << ',' << phaArr[1]
                << "  Range: " << (rangeArr[0] == Event::AcdDigi::LOW ? 'L' : 'H') << ',' 
                << (rangeArr[1] == Event::AcdDigi::LOW ? 'L' : 'H')
                << "  Veto: " << (vetoArr[0] ? 'T' : 'F') << ',' << (vetoArr[1] ? 'T' : 'F')
                << "  CNO: " << (highArr[0] ? 'T' : 'F') << ',' << (highArr[1] ? 'T' : 'F')
                << endreq;  
        }
    }
    return StatusCode::SUCCESS;
}


StatusCode AcdDigiAlg::finalize() {
  toolSvc()->releaseTool(m_randTool);
  return StatusCode::SUCCESS;
}

void AcdDigiAlg::clear() {

    m_energyDepMap.clear();
    m_peMap.clear();
    m_mipsMap.clear();
    m_statusMap.clear();

}


