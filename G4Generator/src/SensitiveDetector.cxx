//$Header$

#include <map>


#include "detModel/Management/Manager.h"
#include "detModel/Management/IDmapBuilder.h"
#include "detModel/Utilities/PositionedVolume.h"
#include "detModel/Utilities/Color.h"
#include "detModel/Sections/Box.h"
#include "detModel/Materials/MatCollection.h"
#include "detModel/Materials/Material.h"

#include "HepRepXMLWriter.h"

#include "DetectorConstruction.h"
#include "SensitiveDetector.h"

#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4GRSSolid.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

/**
    This class is used by G4 to retrive hits information during the simulation.
    When a hit occurs in a sensitive detector, the method ProcessHits is called.
    The methods Initialize and EndOfEvent are called before and after the 
    processing of a full event. We are not actually using the mechanism of Hit
    collections of G4, since we want to fill the TDS directly. For now
    there is this SensitiveDetector for all the detectors, but we can think about
    having 3 separate sensitive detectors classes for ACD, TKR and CAL
 */
SensitiveDetector::SensitiveDetector(G4String name,
                                   DetectorConstruction* det)
  :G4VSensitiveDetector(name),m_detector(det)
{
  /// Temp HepRep initialization stuff
  m_hepRepXMLWriter = new HepRepXMLWriter(); 
  m_hepRepXMLWriter->open("GlastHepRep.xml");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SensitiveDetector::~SensitiveDetector()
{
  /// Temp HepRep finalization stuff
  m_hepRepXMLWriter->close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
  m_energySummary.clear();
  /// Add a new Hits type for the HepRep
  m_hepRepXMLWriter->addType("Hits");  
  /// Clear the vector of detectors ID
  m_hitID.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool SensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  typedef std::map<G4VPhysicalVolume*, std::string> M;
  M::const_iterator i; 

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
  G4String name = "";
  G4String particle = aStep->GetTrack()->GetDefinition()->GetParticleName();

  do
    { 
      i = m_detector->g4ID.find(physVol);
      if(i != m_detector->g4ID.end()){
	name = i->second + name;
      }
      theTouchable->MoveUpHistory();
      physVol = theTouchable->GetVolume(); 
    } while(physVol->GetName() != "motherVolume");

  
  // Initial Position
  G4ThreeVector InitPos = aStep->GetPreStepPoint()->GetPosition();
  
  // Final Position 
  G4ThreeVector FinPos = aStep->GetPostStepPoint()->GetPosition();

    /** HepRep stuff for hits
	It is temporary, since this stuff does not belong to SensitiveDetector, but to an 
	HepRepSvc probably; it is here now since in SensitiveDetector
	I have all the info on the hits.
    */
  m_hepRepXMLWriter->addAttDef("Energy","Energy released in the hit","Physics","MeV");
  m_hepRepXMLWriter->addAttDef("detID","Name of the detector of the hit","Physics","");
  m_hepRepXMLWriter->addAttDef("particle","Name of the particle","Physics","");

  m_hepRepXMLWriter->addInstance();
  m_hepRepXMLWriter->addAttValue("DrawAs","Point");
  m_hepRepXMLWriter->addAttValue("Energy",edep);
  m_hepRepXMLWriter->addAttValue("detID",name);
  m_hepRepXMLWriter->addAttValue("particle",particle);
  m_hepRepXMLWriter->addPrimitive();
  m_hepRepXMLWriter->addPoint(InitPos.x(),
			    InitPos.y(),
			    InitPos.z());

  /// I put in a vector the ID of detectors; it will be used for detector rep.
  m_hitID.push_back(name);

#if 1
  std::cout << "Hit -> " 
	    << std::setw(8) << particle  
	    << std::setw(8) << std::setprecision(3) << edep  
	    << std::setw(8) << std::setprecision(3) <<stepl 
	    << std::setw(18) << name << " " 
	    << std::setw(12) << nameVolume << " " << material
	    << std::endl;
#endif

  m_energySummary[logVol->GetName()] += edep;
  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{

    std::cout << "Energy deposit summary: " << std::endl;
    for( std::map<std::string, double>::const_iterator it = m_energySummary.begin(); it!= m_energySummary.end(); ++it){
        std::cout 
            << std::setw(20) << (*it).first 
            << std::setw(10) << std::setprecision(3) << (*it).second << std::endl;
    }




    /** HepRep stuff for detector
	I would like to use a hierarchy (subType in HepRep), but it seems
	not working with the actual version of WIRED .. need to ask Joseph
	For now the detectors are simply positioned without any hierarchy
	This can be moved to a new method. In any case it is temporary, 
	since this stuff does not belong to SensitiveDetector, but to an 
	HepRepSvc probably; it is here now since in SensitiveDetector
	I have all the info on the hits.
    */

    /// Retrive the detModel manager and the Gdd
    detModel::Manager* manager = detModel::Manager::getPointer();
    detModel::Gdd* gdd = manager->getGdd();

    /// Get the colors map from detModel
    std::map<std::string,detModel::Color*> colorsMap = gdd->getMaterials()->getMaterialColors();

    /// Set up some vectors and matrix useful
    G4ThreeVector trans;
    G4RotationMatrix rot;
    G4ThreeVector vertex;

    /// Build the idMap with detModel: this can be put in the constructor 
    detModel::IDmapBuilder* mapBuilder = new detModel::IDmapBuilder("");
    std::map<std::string,detModel::PositionedVolume*> idMap;
    manager->startVisitor(mapBuilder);
    idMap = mapBuilder->m_volMap;
    delete mapBuilder;

    /// Start to fill the HepRep file detector part
    m_hepRepXMLWriter->addType("Detector");  
    m_hepRepXMLWriter->addAttDef("ID","Name of the detector","Physics","");
    m_hepRepXMLWriter->addAttDef("mat","Material of the detector","Physics","");

    /// Cycle on the hitID vector elements
    for(unsigned int i = 0; i<m_hitID.size();i++)
      {
	/// This is dangerous .. can we have non box sensitive detector? 
	detModel::Box* box = dynamic_cast<detModel::Box*>(idMap[m_hitID[i]]->getVolume());
	    
	/// Retrive geometric info
	double x = box->getX()/2;
	double y = box->getY()/2;
	double z = box->getZ()/2;
	/// Retrive positioning info
	trans = idMap[m_hitID[i]]->getTranslation();
	rot = (idMap[m_hitID[i]]->getRotation());

	/// Add the HepRep to the file
	/// Only the XY projection, just for test
	/// I would like to use DrawAs Prism, but it seems that it does not work 
	m_hepRepXMLWriter->addInstance();
	m_hepRepXMLWriter->addAttValue("DrawAs","Polygon");
	m_hepRepXMLWriter->addAttValue("ID",m_hitID[i].c_str());
	m_hepRepXMLWriter->addAttValue("mat",box->getMaterial().c_str());
	m_hepRepXMLWriter->addAttValue("LineColor",
				       colorsMap[box->getMaterial()]->getRed(),
				       colorsMap[box->getMaterial()]->getGreen(),
				       colorsMap[box->getMaterial()]->getBlue());

	m_hepRepXMLWriter->addPrimitive();

	vertex = rot*G4ThreeVector(-x,-y,-z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
	vertex = rot*G4ThreeVector(x,-y,-z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
	vertex = rot*G4ThreeVector(x,y,-z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
	vertex = rot*G4ThreeVector(-x,y,-z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());

	m_hepRepXMLWriter->addPrimitive();

	vertex = rot*G4ThreeVector(-x,-y,z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
	vertex = rot*G4ThreeVector(x,-y,z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
	vertex = rot*G4ThreeVector(x,y,z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
	vertex = rot*G4ThreeVector(-x,y,z) + trans;
	m_hepRepXMLWriter->addPoint(vertex.x(),vertex.y(),vertex.z());
      }



}





