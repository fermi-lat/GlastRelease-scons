// $Header$
#include "IRFConverter.h"


#include "instrument/Scintillator.h"
#include "instrument/SiDetector.h"
#include "instrument/CsIDetector.h"
#include "instrument/MCTruth.h"
#include "instrument/Tower.h"


#include "src/data/TdCsIData.h"
#include "src/data/TdSiData.h"
#include "src/data/TdVetoData.h"


//! constructor - create all object containers
IRFConverter::IRFConverter() {
    // Added because not sure as of yet how to make the
    // converter without using object vectors
    allcsiData =  new TdCsIData(9);
    allsiData = new TdSiData(18); //TODO: Change hardcode.
    allvetoData = new TdVetoData();

    IrfAcdHitContainer = new IrfAcdHitVector;
    IrfCalHitContainer = new IrfCalHitVector;
    IrfTkrLayerContainer = new IrfTkrLayerVector;
}

//! destructor - delete all object containers
IRFConverter::~IRFConverter() { 
    if (IrfAcdHitContainer) delete IrfAcdHitContainer;
    if (IrfCalHitContainer) delete IrfCalHitContainer;
    if (IrfTkrLayerContainer) delete IrfTkrLayerContainer;
}

//! called due to a GlastDetector::accept(DetectorConverter) call
void IRFConverter::forward (const Scintillator& s) {
    // ACD tile data
    // retrieve data only if this detector has data
    if ( !(s.empty()) ) {
        IrfAcdHit* newTile = new IrfAcdHit();
        newTile->setEnergy(s.energy());
        newTile->setId(s.id());
        IrfAcdHitContainer->push_back(newTile);
    }

    allvetoData->load(s);
}

// just grab the ID to pass onto subsequent visits
void IRFConverter::forward (const Tower& t) {
    m_towerId = t.getId();
}

//! called due to GlastDetector::accept(), handles Cal data
void IRFConverter::forward ( const CsIDetector& csi) {
    // CAL CsI log data
    if ( !(csi.empty()) ) {
        
        IrfCalHit* irfCal = new IrfCalHit();
        irfCal->setEnergy(csi.energy());
        irfCal->setMinusResponse(csi.Lresp());
        irfCal->setPlusResponse(csi.Rresp());
        IrfCalHitContainer->push_back(irfCal);


        //Do the Raw information
        allcsiData->load(csi, m_towerId);
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

        IrfTkrLayer* tkrLayer = new IrfTkrLayer();
        tkrLayer->setId(si.id());
        tkrLayer->setMaxEnergy(si.maxELoss());

        // loop over all hit strips for this layer
        for (SiDetector::const_iterator it = si.begin(); it != si.end(); ++it) 
        {
            IrfTkrHit* newSSD = new IrfTkrHit();
            newSSD->setEnergy((*it).energy());
            newSSD->setId((*it).index());
            newSSD->setNoise((*it).noise());
            tkrLayer->addHit(newSSD);
        }
        
        IrfTkrLayerContainer->push_back(tkrLayer);

        allsiData->load(si,m_towerId);
    }
}

//! provide access to the Raw CsI  Data geometry included as well
TdCsIData* IRFConverter::getTdCsIData() { return allcsiData; }

TdSiData* IRFConverter::getTdSiData() { return allsiData; }

TdVetoData* IRFConverter::getTdVetoData() { return allvetoData;}
//! Provide access to the IRF ACD data
IrfAcdHitVector* IRFConverter::getIrfAcdHits() {return IrfAcdHitContainer; }
//! Provide access to the IRF Cal data
IrfCalHitVector* IRFConverter::getIrfCalHits() {return IrfCalHitContainer; }
//! Provide access to the IRF Tkr data
IrfTkrLayerVector* IRFConverter::getIrfTkrHits() {return IrfTkrLayerContainer; }