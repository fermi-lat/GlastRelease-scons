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
            
            DisplayManager::instance()->addBox(global, x,y,z);
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
};
//------------------------------
class TkrPlane : public GenericDet {
public:
    TkrPlane(SiDetector& si): m_si(si){};
    void score(G4Step * aStep){
    }
private:
    SiDetector& m_si;
};
//------------------------------
class CalLog : public GenericDet {
public:
    CalLog(CsIDetector& csi): m_csi(csi){}
    void score(G4Step * aStep){
    }
private:
    CsIDetector& m_csi;
};
//------------------------------
class Diode : public GenericDet {
public:
    Diode(DiodeDetector& diode):m_diode(diode){}
    void score(G4Step * aStep){
    }
private:
    DiodeDetector& m_diode;
};
//------------------------------
class ACDtile : public GenericDet {
public:
    ACDtile(Scintillator& acd):m_acd(acd){}
    void score(G4Step * aStep){
    }
private:
    Scintillator& m_acd;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GlastDetectorManager::GlastDetectorManager(DetectorConstruction *det)
:G4VSensitiveDetector("GlastDetectorManager")
{
    m_idMap = det->idMap();
#if 0
    static std::string detfile("$(INSTRUMENTROOT)/xml/GlastDetector.xml");
    std::cout << "Loading glastdetectors from " << detfile << "... ";
    m_instrument = new Instrument;
    m_instrument->initialize("",detfile );
    std::cout << "found " << m_instrument->detector_count() 
        << " GlastDetector objects" << std::endl;
    
    // now use the local visitor to make a list
    m_instrument->rootDetector()->accept(Visitor(this));
    //std::cout << "actually found " << count << " detectors " << std::endl;
#endif
    // and tell G4 about us
    G4SDManager::GetSDMpointer()->AddNewDetector( this );

    
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetectorManager::Initialize(G4HCofThisEvent*HCE)
{
    m_energySummary.clear();
    m_detectorList.clear();
    m_detectorEnergy.clear();
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

G4bool GlastDetectorManager::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
    
    // Energy Deposition & Step Lenght
    
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
    
    //**** interface to display *************
    
    DisplayManager::instance()->addHit(InitPos, FinPos);
    
    idents::VolumeIdentifier id = constructId(aStep);
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
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GlastDetectorManager::process(G4LogicalVolume* lvol)
{
    lvol->SetSensitiveDetector(this);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GlastDetectorManager::addSiDetector(SiDetector* det)
{
    m_detMap[det->id()]=new TkrPlane(*det);

}
void GlastDetectorManager::addCsIDetector(CsIDetector* det)
{
    m_detMap[det->id()]=new CalLog(*det);
    
}
void GlastDetectorManager::addDiodeDetector(DiodeDetector* det)
{
    m_detMap[det->id()] = new Diode(*det);
}
void GlastDetectorManager::addACDtile(Scintillator* det)
{
    m_detMap[det->id()] = new ACDtile(*det);
}
