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
 * the detector geometry of the G4 simulation. In particular it defines
 *
 *   - full detector geometry
 *   - define materials
 *   - identify sensitive detectors, set up callbacks
 *   - create a map of ids and physical volumes
 *
 *  
 * @author R.Giannitrapani
 *    
 * $Header$
 */

class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier > IdMap;

  //! @param gds A pointer to the abstract interface of the GlastDetSvc
  //! @param esv A pointer to the data provider service
  DetectorConstruction(IGlastDetSvc* gds, IDataProviderSvc* esv);

  ~DetectorConstruction();
  
  /* Actual call to construct things. Please note that this is the concrete
   * implementation of the abstract method of G4VUserDetectorConstruction class
   * (geant4 name convention). 
   * @return The mother G4VPhysicalVolume for the full detector
   */
  G4VPhysicalVolume* Construct();
  
  //! Return to the map of physical volume/id pairs.
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
};

#endif

