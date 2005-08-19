// File and Version Information:
// $Header$
//
// Description: this is a concrete implementation of the abstract interface
// IGeometry of GlastDetSvc. Its methods are called by the visitor mechanism of
// detModel encapsulated in the GlastDetSvc. The use of an abstract interface
// permits to avoid the explicit exposition of detModel interfaces to
// clients. This concrete implementation uses Geant4 classes to build the
// volumes hierarchy for G4Generator. For more detailed information, please see
// the GlastDetSvc documentation
//      
// Author(s):
//      R.Giannitrapani

#include "G4Geometry.h"
#include "PosDetectorManager.h"
#include "IntDetectorManager.h"

#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"

#include "DetectorConstruction.h"
#include "LocalMagneticFieldDes.h"

// for the display
#include "gui/GuiMgr.h"
#include "gui/DisplayControl.h"
#include "gui/DisplayRep.h"
#include "geomrep/TubeRep.h"
#include "geomrep/BoxRep.h"
#include <iomanip>
#include <cassert>

G4Geometry::G4Geometry(std::string mode) : m_mode(mode)
{
  // Inputs: the mode specify the level of details to be used in building the
  // geometry
}

G4Geometry::~G4Geometry()
{

}

IGeometry::VisitorRet 
G4Geometry::pushShape(ShapeType s, const UintVector& idvec, 
                      std::string name, std::string material, 
                      const DoubleVector& params, VolumeType type,
                      SenseType sense)
{
  // Purpose and Method: this method push a new volume in the stack of volumes
  // Inputs: for the meaning of ShapeType, the VolumeType and the params vector
  // see the GlastDetSvc documentation. T
  // Outputs: the VisitorRet can be used to stop the recursion during geometry
  // building

  // The first 6 parameters in params are the translations and rotations
  double x=params[0], y=params[1], z=params[2];
  double rx=params[3], ry=params[4], rz=params[5];  

  // the physical and logical volumes
  G4VPhysicalVolume* phys;
  G4LogicalVolume* logical;

  // Let's check if the logical is already there
  if (!(logical = m_logicals[name]))
    {
      // the shape
      G4VSolid* shape;
      // Get the material
      G4Material* ptMaterial = G4Material::GetMaterial(material);
      if (!ptMaterial) ptMaterial = G4Material::GetMaterial("Vacuum");
    
      // Build a box or a tube
      if( s==Box) {
        double dx=params[6], dy=params[7], dz=params[8];
        // if there is no actualMother, it means this is the world volume
        // we use at the moment a box of 60m*60m*60m
        if (!actualMother())
          { dx = 60000; dy = 60000; dz = 60000;}
        shape = new G4Box(name,
                          dx*mm/2,
                          dy*mm/2,
                          dz*mm/2);    
      }else if(s==Tube) {
        double dz=params[6], rmin=params[7], rmax=params[8];
        shape = new G4Tubs(name,
                           rmin*mm,
                           rmax*mm,
                           dz*mm*0.5,
                           0,2*M_PI);
      }
      // put the logical in the m_logicals vector
      logical = new G4LogicalVolume(shape,ptMaterial,name,0,0,0);
      m_logicals[name] = logical; 
    }

  // produce local magnetic field 
  if(logical && m_det && m_det->getFieldDes().m_magFieldVol == name) {
	logical->SetFieldManager(m_det->getFieldManager(), true);
  }


  /*
   * In order to avoid duplication of physical volumes in 
   * multipositioned ensemble volume (like stack or composition)
   * we check if we are in a replicated volume by looking at the
   * boolean variable m_replica and if a physical volume with the
   * same name has been already created
   */
  if ((m_replica) && (getPhysicalByName(name)))
    {
      // we push this logical as the actual mother and we exit
      pushActualMother(logical);
      return More;
    }
  
  /*
   * If this volume is already in the physicals vector, we need
   * to avoid to duplicate its contained subvolumes; for this
   * we set the boolean flag m_replica to 1 and we set the m_replicaMother
   * to the actualMother .. this can be used in the popShape method to check
   * if we have to disable the m_replica flag
   */
  if (getPhysicalByName(name)) 
    {
      m_replicaMother = actualMother();
      m_replica = 1;
    }

  // Set the rotation: note conversion to radians and change of sign
  // The sign is changed for different conventions between our XML
  // geometry database and G4 ones.
  G4RotationMatrix* rm = new G4RotationMatrix(); 
  rm->rotateX(-rx*M_PI/180);
  rm->rotateY(-ry*M_PI/180);
  rm->rotateZ(-rz*M_PI/180);

  // Create the positioned physical volume
  phys = new G4PVPlacement(rm,G4ThreeVector(x*mm,y*mm,z*mm),
                           logical,
                           name,
                           actualMother(),false,0);

  // push it in the physical volumes vector
  m_physicals.push_back(phys);

  // if there is no actualMother this is the world volume
  if (!actualMother()) setWorld(phys);

  // Set the id
  idents::VolumeIdentifier vid;
  for( UintVector::const_iterator ui=idvec.begin(); ui !=idvec.end(); ++ui) 
    {
      vid.append((unsigned int)*ui);  
    }

  // Fill the physicals-id map
  (*m_idMap)[phys] = vid;

  // Register the volume to the sensitive detector if necessary
  if (sense == posSensitive) 
    {
      m_pdm->process(logical);
    }
  else if (sense == intSensitive)
    {
      m_idm->process(logical);
    }

  // push this volume in the stack of mothers
  pushActualMother(logical);    
  return More;
}


void G4Geometry::popShape()
{
  // Purpose and Method: called to signal end of nesting 

  // pop the actual mother stack
  popActualMother();

  // if we were in a replicated volume, check if we have to set
  // the replica flag back to 0
  if (m_replica)
    if (actualMother() == m_replicaMother)
      m_replica = 0;
}


G4VPhysicalVolume* G4Geometry::getPhysicalByName(std::string name){
  // Purpose and Method: a utility method to search a physical volume by name
  
  std::vector<G4VPhysicalVolume*>::const_iterator i;
  
  for(i=m_physicals.begin();i!=m_physicals.end();i++)
    if ((std::string)(*i)->GetName() == name) return (*i);

  return 0;
}

G4LogicalVolume* G4Geometry::actualMother()const{
  // Purpose and Method: this one return the actual mother volume

  if (m_actualMother.size()) 
    return m_actualMother.back();
  else return 0;
}
