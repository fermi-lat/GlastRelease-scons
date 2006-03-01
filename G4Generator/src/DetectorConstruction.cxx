// File and Version Information:
// $Header$
//
// Description: This class hinerits from G4VUserDetectorConstruction and it is
// used to build the detector geometry of the G4 simulation by using two other
// classes, G4Media and G4Geometry. Both are concrete implementation of
// GlastDetSvc abstract intefaces and are used to set up materials table and
// volumes hierarchy respectively
//
// Author(s):
//      R.Giannitrapani

#include "DetectorConstruction.h"

// geant4 include files
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SDManager.hh"

// sensitive detectors managers
#include "PosDetectorManager.h"
#include "IntDetectorManager.h"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

// visitors interfaces
#include "G4Geometry.h"
#include "G4Media.h"

#include "LocalMagneticFieldDes.h"

#include <iomanip>
#include <cassert>

DetectorConstruction::DetectorConstruction(IGlastDetSvc* gsv,
                                           IDataProviderSvc* esv,
                                           std::string geometry_mode,
                                           std::ostream& log,
					   const LocalMagneticFieldDes& 
					   magFieldDes)
                                           :m_gsv(gsv), 
                                           m_geometryMode(geometry_mode),
  m_log(log), m_magFieldDes(magFieldDes), m_magField(0), m_fieldManager(0),
  m_gMagField(0)
{
  // Purpose and Method: the constructor needs the GlastDetSvc and
  // DataProviderSvc; the first one is used to start the visitor mechanism using
  // G4Geometry and G4Media, while the DataProviderSvc is passed to the
  // constructors of the sensitive detectors managers in such a way that they
  // will be able to save hits in the TDS
  
  // now create the GlastDetector managers
  m_posDet = new PosDetectorManager(this, esv, gsv);
  m_intDet = new IntDetectorManager(this, esv, gsv);

  // construt local magnetic field. It seems we need to build a global magnetic
  // field first (even when it is 0) in order for the local field to work
  if(m_magFieldDes.m_magFieldVol != "") {
    m_gMagField = new G4UniformMagField(G4ThreeVector(0., 0., 0.));
    G4FieldManager* temp = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    temp->SetDetectorField(m_gMagField);
    temp->CreateChordFinder(m_gMagField);
    temp->GetChordFinder()->SetDeltaChord(0.1*mm);
  
    G4double x = 1*tesla;
    G4double y = 1*tesla;
    G4double z = 1*tesla;

    x *= m_magFieldDes.m_magFieldX;
    y *= m_magFieldDes.m_magFieldY;
    z *= m_magFieldDes.m_magFieldZ;

    m_magField = new G4UniformMagField(G4ThreeVector(x, y, z));
    m_fieldManager = new G4FieldManager(m_magField);
    m_fieldManager->CreateChordFinder(m_magField);
    m_fieldManager->GetChordFinder()->SetDeltaChord(0.1*mm);
  }
}

DetectorConstruction::~DetectorConstruction()
{ 
  delete m_posDet;
  delete m_intDet;

  if(m_magField) delete m_magField;
  if(m_fieldManager) delete m_fieldManager;
  if(m_gMagField) delete m_gMagField;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Purpose and Method: this is the method called by Geant4 to build the
  // geometry 
  // Outputs: a pointer to the world physical volume

  G4Geometry geom(m_posDet, m_intDet, &m_idMap, 
		  m_geometryMode, this);
  G4Media media;
  m_gsv->accept(media);
  m_gsv->accept(geom);

  // summary log output
  m_log << "\tDetectorConstruction created "
      << geom.getPhysicalNumber() << ", "<< geom.getLogicalNumber()<< " physical, logical volumes, using ";
  if (m_geometryMode=="" ) {
      m_log <<  "default mode from GlastDetSvc"  << std::endl;
  }else {
      m_log << " mode from G4Generator job options: " << m_geometryMode << std::endl;
  }

  return geom.getWorld();
}
















