// $Header$

#include "GlastDetectorManager.h"
#include <iostream>
#include "instrument/Instrument.h"
#include "instrument/DetectorVisitor.h"
#include "instrument/GlastDetector.h"
#include "instrument/SiDetector.h"
#include "instrument/CsIDetector.h"
#include "instrument/DiodeDetector.h"
#include "instrument/Scintillator.h"
#include "CLHEP/Geometry/Transform3D.h"
#include "geometry/Box.h"
#include "GlastEvent/MonteCarlo/McPositionHit.h"
#include "idents/VolumeIdentifier.h"

#include "McHitsManager.h"

// Geant4 interface
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"


namespace {
    class Visitor : public DetectorVisitor {
    public:
        Visitor(GlastDetectorManager* dm):m_dm(dm){}
        virtual void visit (Glast* )
        {}
        
        
        virtual void visit (SiDetector* si) {
            m_dm->addSiDetector(si);
        }
        
        virtual void visit (CsIDetector* csi){
            m_dm->addCsIDetector(csi);
            
        }
        
        virtual void visit (DiodeDetector* diode){
            m_dm->addDiodeDetector(diode);
        }
        
        virtual void visit (Scintillator* acd ) {
            m_dm->addACDtile(acd);
        }
        // leave these for reference: not needed
        virtual void visit (Tower* ){ }
        virtual void visit (MCTruth* ) {}
        virtual void visit (Calorimeter* cal) { }
        virtual void visit (SiTracker* st) { }
    private:
        GlastDetectorManager* m_dm;
    };
    // get box params (dimensions, transformation) from the touched guy.
    
    void makeBox(G4TouchableHistory* touched)
    {
        G4VPhysicalVolume* pvol = touched->GetVolume(); 
        
        HepTransform3D 
            global(*(touched->GetRotation()), 
            touched->GetTranslation());
        
        
        const G4LogicalVolume* lvol = pvol->GetLogicalVolume();
        const G4VSolid * solid = lvol->GetSolid();
        const G4Box* box = dynamic_cast<const G4Box*>(solid);
        if( box !=0){
            double 
                x = 2*box->GetXHalfLength(), 
                y = 2*box->GetYHalfLength(), 
                z = 2*box->GetZHalfLength();
            
            DisplayManager::instance()->addHitBox(global, x,y,z);
        }
        
    }
    
}
//------------------------------
/**
Base class for conversion objects to convert G4Step info into appropriate calls to
the corresponding GlastDetector objects.
  */
class GenericDet {
public:
    GenericDet(){}
    //! pure virtual to be implemented by subclasses.
    virtual void score(G4Step * aStep)=0;
    //! utillity to extract and convert energy
    double depositedEnergy(G4Step * aStep){
        return aStep->GetTotalEnergyDeposit()*(MeV/GeV);
    }
};
//------------------------------
class TkrPlane : public GenericDet {
public:
    TkrPlane(SiDetector& si): m_si(si){};
    void score(G4Step * aStep){

        Hep3Vector beginpt = aStep->GetPreStepPoint()->GetPosition()*(mm/cm),
                   endpt =aStep->GetPostStepPoint()->GetPosition()*(mm/cm);

        Point p(beginpt.x(),beginpt.y(), beginpt.z()),
              q( endpt.x(),  endpt.y(),   endpt.z() );
        p.transform(m_si.globalToLocal());
        q.transform(m_si.globalToLocal());
        // NOTE conversion to "fC" units, just like before
        m_si.score(p,q, depositedEnergy(aStep)/ 3.4875E-5);
    }
private:
    SiDetector& m_si;
};
//------------------------------
class CalLog : public GenericDet {
public:
    CalLog(CsIDetector& csi): m_csi(csi){}

    void score(G4Step * aStep){
        // get begin and endpoints in cm units.
        Hep3Vector beginpt = aStep->GetPreStepPoint()->GetPosition()*(mm/cm),
                   endpt =aStep->GetPostStepPoint()->GetPosition()*(mm/cm),
                   midpt = 0.5*(beginpt+endpt);

        Point p(midpt.x(),midpt.y(),midpt.z()); // stupid way to make it a Point
        p.transform(m_csi.globalToLocal());

        m_csi.scoreDetector (depositedEnergy(aStep),  p);
    }
private:
    CsIDetector& m_csi;
};
//------------------------------
class Diode : public GenericDet {
public:
    Diode(DiodeDetector& diode):m_diode(diode){}
    void score(G4Step * aStep){
        double deltaE = aStep->GetTotalEnergyDeposit()*(MeV/GeV);
        m_diode.addEnergy(depositedEnergy(aStep)); 
    }
private:
    DiodeDetector& m_diode;
};
//------------------------------
class ACDtile : public GenericDet {
public:
    ACDtile(Scintillator& acd):m_acd(acd){}
    void score(G4Step * aStep){
        m_acd.addEnergy(depositedEnergy(aStep));
    }
private:
    Scintillator& m_acd;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GlastDetectorManager::GlastDetectorManager(DetectorConstruction *det, const detModel::IDmapBuilder& idmap)
:G4VSensitiveDetector("GlastDetectorManager")
{

    // set up iterators over the map if VolumeId's of all detectors.
    m_id_it = idmap.getIdVector()->begin(); // will be incremented
    m_id_end = idmap.getIdVector()->end();  // used to prevent too many: must equal m_id_it at end.

    m_idMap = det->idMap();

    // construct file name for xml file with detector stuff
    std::string detfile("$(INSTRUMENTROOT)/xml/");
    detfile += det->topVolumeName()+".xml";

    std::cout << "Loading glastdetectors from " << detfile << "... ";
    m_instrument = new Instrument;
    m_instrument->initialize("",detfile );
    std::cout << "found " << m_instrument->detector_count() 
        << " GlastDetector objects" << std::endl;
    
    // now use the local visitor to make a list
    m_instrument->rootDetector()->accept(Visitor(this));
    std::cout << "actually found " << m_detMap.size() << " detectors " << std::endl;
    if( m_id_it != m_id_end ) {
        std::cerr << "Warning!!! Did not use all of the ids in the list of VolumeIdentifiers" << std::endl;
    }
    // and tell G4 about us
    G4SDManager::GetSDMpointer()->AddNewDetector( this );

    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GlastDetectorManager::~GlastDetectorManager()
{
    delete m_instrument;
}
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetectorManager::Initialize(G4HCofThisEvent*HCE)
{
    m_energySummary.clear();
    m_detectorList.clear();
    m_detectorEnergy.clear();
    // clear all glast detectors
    GlastDetector::clear_all();

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

G4bool GlastDetectorManager::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
    
    // Energy Deposition & Step Length
    
    G4double edep = aStep->GetTotalEnergyDeposit()/MeV;
    G4double stepl = aStep->GetStepLength()/mm;
    
    if ((edep==0.)) return false;      
    
    // Physical Volume
    
    G4TouchableHistory* theTouchable
        = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
    G4LogicalVolume* logVol = physVol->GetLogicalVolume();
    G4String material = logVol->GetMaterial()->GetName();
    G4String nameVolume = physVol->GetName();
    
    G4ThreeVector InitPos = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector FinPos = aStep->GetPostStepPoint()->GetPosition();
    
    // determine the ID by studying the history, then call appropriate 
    idents::VolumeIdentifier id = constructId(aStep);
    assert(m_detMap[id]);
    m_detMap[id]->score(aStep);


    // Filling of the hits container
    mc::McPositionHit *hit = new mc::McPositionHit;
    
    hit->setDepositedEnergy(edep);
    hit->setVolumeID(id);
    hit->setEntryPoint(InitPos);
    hit->setExitPoint(FinPos);

    // We add the hit to the manager
    McHitsManager* mchits = McHitsManager::getPointer();
    mchits->addHit(hit);
 

    //**** interface to display *************
    
    DisplayManager::instance()->addHit(InitPos, FinPos);
    
    if( m_detectorList[id]==0) {
        makeBox( theTouchable );
        
    }
    ++ m_detectorList[id]; 
    m_detectorEnergy[id] += edep;
    
#if 0
    std::cout << "Hit -> " 
        << std::setw(8) << aStep->GetTrack()->GetDefinition()->GetParticleName()  
        << std::setw(8) << std::setprecision(3) << edep  
        << std::setw(8) << std::setprecision(3) <<stepl 
        << std::setw(12) << nameVolume << " " 
        << id.name() 
        << std::endl;
#endif
    m_energySummary[logVol->GetName()] += edep;
    return true;
    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

idents::VolumeIdentifier GlastDetectorManager::constructId(G4Step * aStep)
{
    using  idents::VolumeIdentifier;
    VolumeIdentifier ret;
    G4TouchableHistory* theTouchable
        = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    for( int i = 0; i<theTouchable->GetHistoryDepth() ; ++i) {
        const G4VPhysicalVolume* physVol = theTouchable->GetVolume(i); 
        if( physVol->GetMother()==0) break;
        VolumeIdentifier id = (*m_idMap)[physVol];
        ret.prepend(id);
    }
    return ret;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetectorManager::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    m_instrument->rootDetector()->generateResponse(true);
#if 0
    std::cout << "Detector list:" << std::endl;
    
    for( DetectorList::const_iterator itd = m_detectorList.begin(); itd!=m_detectorList.end(); ++ itd){
        std::cout 
            << std::setw(20) << (*itd).first.name()
            << std::setw(5)  << (*itd).second
            << std::setw(10) <<  std::setprecision(3)<< m_detectorEnergy[(*itd).first]
            << std::endl;
    }
    std::cout << "Energy deposit summary: " << std::endl;
    for( std::map<std::string, double>::const_iterator it = m_energySummary.begin(); it!= m_energySummary.end(); ++it){
        std::cout 
            << std::setw(20) << (*it).first 
            << std::setw(10) << std::setprecision(3) << (*it).second 
            << std::endl;
    }
#endif
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GlastDetectorManager::process(G4LogicalVolume* lvol)
{
    lvol->SetSensitiveDetector(this);
}
idents::VolumeIdentifier GlastDetectorManager::checkId(idents::VolumeIdentifier::int64 gid, const char * name)
{
    assert(m_id_it != m_id_end);
    idents::VolumeIdentifier id = *m_id_it++;
   // std::cout << "associate VolumeIdentifier " << id.name() << " with " << (int)gid << " "<< name << std::endl;
    return id;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetectorManager::addSiDetector(SiDetector* det)
{

    m_detMap[checkId(det->id(),"Si")]=new TkrPlane(*det);

}
void GlastDetectorManager::addCsIDetector(CsIDetector* det)
{
    m_detMap[checkId(det->id(),"CsI")]=new CalLog(*det);
    
}
void GlastDetectorManager::addDiodeDetector(DiodeDetector* det)
{
    m_detMap[checkId(det->id(),"diode")] = new Diode(*det);
}
void GlastDetectorManager::addACDtile(Scintillator* det)
{
    m_detMap[checkId(det->id(),"ACD")] = new ACDtile(*det);
}
