// $Header$
#include "IRFConverter.h"


#include "instrument/Scintillator.h"
#include "instrument/SiDetector.h"
#include "instrument/CsIDetector.h"
#include "instrument/MCTruth.h"


#include "GlastEvent/Raw/TdCsIData.h"


//! constructor - create all object containers
IRFConverter::IRFConverter() {
    MCACDcontainer = new MCACDHitVector;
    MCCalorimeterHitContainer = new MCCalorimeterHitVector;
    MCSiLayerContainer = new MCSiLayerVector;
    
    // Added because not sure as of yet how to make the
    // converter without using object vectors
    allcsiData =  new TdCsIData(9);
}

//! destructor - delete all object containers
IRFConverter::~IRFConverter() { 
    if (MCACDcontainer)  delete MCACDcontainer;
    if (MCCalorimeterHitContainer) delete MCCalorimeterHitContainer;
    if (MCSiLayerContainer) delete MCSiLayerContainer;
}

//! called due to a GlastDetector::accept(DetectorConverter) call
void IRFConverter::forward (const Scintillator& s) {
    // ACD tile data
    // retrieve data only if this detector has data
    if ( !(s.empty()) ) {
        MCACDHit* tile = new MCACDHit(); 
        tile->setEnergy(s.energy());
        tile->setId(s.id());
        MCACDcontainer->push_back(tile);
    } 
}

//! called due to GlastDetector::accept(), handles Cal data
void IRFConverter::forward ( const CsIDetector& csi) {
    // CAL CsI log data
    if ( !(csi.empty()) ) {
        
        // Do the MCInformation
        MCCalorimeterHit* mcCal = new MCCalorimeterHit();
        mcCal->setEnergy(csi.energy());
        mcCal->setLeftResponse(csi.Lresp());
        mcCal->setRightResponse(csi.Rresp());
        MCCalorimeterHitContainer->push_back(mcCal);
        
        //Do the Raw information
        allcsiData->load(csi);
        
        
    }
    
}

//! not implemented
void IRFConverter::forward ( const MCTruth& mc) {
    // monte carlo truth
    /** Going to increase the functionality of
    so that this converter can handle MCTrack class*/
    
    if ( !(mc.empty()) ) {
        
    /* Put on back burner untill Toby looks at it tomorrow
    
      MCTrack* particle = new MCTrack();
      HepLorentzVector* lv = new HepLorentzVector();
      lv->setX(mc.particle()->getIanCrap());
      lv->setY(mc.particle().m_pos.y());
      lv->setZ(mc.particle().m_pos.z());
      // need to get the set up the different 
      //properties here. Just for 
      
        particle->setFourMomentum(lv);
        MCTrackVector->push_back(particle); */
        
    }
}

//! called due to GlastDetector::accept(), handles TKR strip data
void IRFConverter::forward ( const SiDetector& si) {
    // TKR silicon strip detectors
    if ( !(si.empty()) ) {
        MCSiLayer* siLayer = new MCSiLayer();
        siLayer->setId(si.id());
        siLayer->setMaxEnergy(si.maxELoss());
        // loop over all hit strips for this layer
        for (SiDetector::const_iterator it = si.begin(); it != si.end(); ++it) 
        {
            MCTKRHit* ssd = new MCTKRHit();
            ssd->setEnergy((*it).energy());
            ssd->setId((*it).index());                
            ssd->setNoise((*it).noise());
            siLayer->addHit(ssd);
        }
        
        MCSiLayerContainer->push_back(siLayer);
    }
}

//! provide access to the list of ACD tile data
MCACDHitVector* IRFConverter::getACDTiles() { return MCACDcontainer; }
//! provide access to the list of CAL crystal data
MCCalorimeterHitVector* IRFConverter::getMCCalHits() { return MCCalorimeterHitContainer; }
//! provide access to the TKR strip data
MCSiLayerVector* IRFConverter::getMCTKRHits() { return MCSiLayerContainer; }
//! provide access to the Raw CsI  Data geometry included as well
TdCsIData* IRFConverter::getTdCsIData() { return allcsiData; }

