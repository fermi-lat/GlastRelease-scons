#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4VPhysicalVolume;
#ifdef WIN32
# include <float.h>
#endif
#include <string>
#include <vector>
#include <map>
#include "G4VUserDetectorConstruction.hh"
#include "idents/VolumeIdentifier.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "GaudiKernel/IDataProviderSvc.h"

class PosDetectorManager;
class IntDetectorManager;

/** 
 * @class DetectorConstruction
 *
 * @brief A class to build the detector geometry
 *
 * This class hinerits from G4VUserDetectorConstruction and it is used to build
 * the detector geometry of the G4 simulation. In particular 
 *
 *   - it builds the full detector geometry
 *   - it defines the materials table
 *   - it identifies sensitive detectors and set up callbacks
 *   - it creates a map of ids and physical volumes
 *
 * To do this it uses two other classes, G4Geometry and G4Media, that implement
 * abstract interfaces of the GlastSvc package; they setup the Visitor mechanism
 * of detModel to visit the geometry representation built from the xml
 * database. Thanks to the abstraction mechanism implemented in the GlastDetSvc,
 * this can be done without exposing internal representation of detModel.
 *  
 * @author R.Giannitrapani
 *    
 * $Header$
 */

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier > IdMap;

  /** The constructor needs from the caller (RunManager) a pointer to the
   *  IGlastDetSvc that is passed than to the DetectorManager subclasses
   *  (PosDetectorManager and IntDetectorManager); this permits to them to fill
   *  the TDS directly at the end of an event.  
   *  @param gds A pointer to the abstract interface of the GlastDetSvc 
   *  @param esv A pointer to the data provider service
   *  @param geometry_mode The mode to use for detModel
   */
  DetectorConstruction(IGlastDetSvc* gds, IDataProviderSvc* esv, std::string geometry_mode="recon");

  ~DetectorConstruction();
  
  /** Actual call to construct things. Please note that this is the concrete
   * implementation of the abstract method of G4VUserDetectorConstruction class
   * (geant4 name convention). 
   * @return The mother G4VPhysicalVolume for the full detector
   */
  G4VPhysicalVolume* Construct();
  
  //! Return the map of physical volume/id pairs.
  IdMap* idMap(){return &m_idMap;}

  //! Return the name of the volume used as mother
  const std::string & topVolumeName()const {return m_topvol;}
  
private:
  /// A pointer to the GlastDetSvc
  IGlastDetSvc* m_gsv;
  
  /// The mother volume of the full geometry
  std::string m_topvol;
  
  //! Map of physical volumes indicized by identifiers
  IdMap m_idMap;
  
  //! The sensitive detector manager of McPositionHit
  PosDetectorManager* m_posDet;

  //! The sensitive detector manager of McIntegratingHit
  IntDetectorManager* m_intDet;

  //! the mode to use when constructing the geometry
  std::string m_geometryMode;
};

#endif

