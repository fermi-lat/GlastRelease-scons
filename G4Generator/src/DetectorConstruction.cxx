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
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SDManager.hh"

// sensitive detectors managers
#include "PosDetectorManager.h"
#include "IntDetectorManager.h"

// visitors interfaces
#include "G4Geometry.h"
#include "G4Media.h"

#include <iomanip>
#include <cassert>

DetectorConstruction::DetectorConstruction(IGlastDetSvc* gsv,
                                           IDataProviderSvc* esv,
                                           std::string geometry_mode,
                                           std::ostream& log)
                                           :m_gsv(gsv), 
                                           m_geometryMode(geometry_mode),
                                           m_log(log)
{
  // Purpose and Method: the constructor needs the GlastDetSvc and
  // DataProviderSvc; the first one is used to start the visitor mechanism using
  // G4Geometry and G4Media, while the DataProviderSvc is passed to the
  // constructors of the sensitive detectors managers in such a way that they
  // will be able to save hits in the TDS
  
  // now create the GlastDetector managers
  m_posDet = new PosDetectorManager(this, esv);
  m_intDet = new IntDetectorManager(this, esv);
}

DetectorConstruction::~DetectorConstruction()
{ 
  delete m_posDet;
  delete m_intDet;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{  
  // Purpose and Method: this is the method called by Geant4 to build the
  // geometry 
  // Outputs: a pointer to the world physical volume

  G4Geometry geom(m_posDet, m_intDet, &m_idMap, 
                                    m_geometryMode);
  G4Media media;
  m_gsv->accept(media);
  m_gsv->accept(geom);

  // summary log output
  m_log << "\tDetectorConstruction created "
      << geom.getPhysicalNumber() <<  " physical volumes, using ";
  if (m_geometryMode=="" ) {
      m_log <<  "default mode from GlastDetSvc"  << std::endl;
  }else {
      m_log << " mode from G4Generator job options: " << m_geometryMode << std::endl;
  }

  return geom.getWorld();
}
















