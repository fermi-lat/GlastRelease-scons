#ifndef DetectorManager_h
#define DetectorManager_h
#ifdef WIN32 // for G4 
#include <float.h>
#endif

#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"

#include "DisplayManager.h"

#include "idents/VolumeIdentifier.h"
#include "DetectorConstruction.h"

#include "GlastEvent/MonteCarlo/McPositionHit.h"
#include "GaudiKernel/IDataProviderSvc.h"

#include <map>
class G4TouchableHistory;
class DisplayManager;
namespace mc {class McPositionHit;}

/** 
 * @class DetectorManager
 *
 * @brief An abstract class for sensitive detectors definition
 *
 * An abstract class whose concrete implementations manage sensitive detectors
 * for the simulation. Two main subclasses are provided:
 *
 *     - PosDetectorManager for dealing with sensitive detectors that must 
 *       register hits as McPositionHit (mainly silicon plane of the TKR)
 *     - IntDetectorManager for dealing with sensitive detectors that must 
 *       register hits as McIntegratingHit (ACD tiles and CAL logs)
 *  
 * @author T.Burnett and R.Giannitrapani
 *    
 * $Header$
 */
class DetectorManager : public G4VSensitiveDetector {
 public:

  /** The IDataProviderSvc is used at the end of the event to store the
   *  collected hits in the TDS
   *  @param idMap map of volume ids for all sensitive detectors
   *  @param esv the data provider service for TDS access
   *  @param name the name to be used by the sensitive detector manager
   */
  DetectorManager( DetectorConstruction::IdMap* idMap, 
                   IDataProviderSvc* esv, 
                   std::string name);
  
  ~DetectorManager();

  //! Called from DetectorConstruction to register the sensitive detector to a
  //! logical volume
  void process(G4LogicalVolume*);

 protected:
  /// This method build the volume identifier from information contained 
  /// in a G4step.
  idents::VolumeIdentifier constructId(G4Step * aStep);
  
  /// generate the box parameters to display
  void makeDisplayBox(G4TouchableHistory* touched,
                      idents::VolumeIdentifier id,
                      bool hitBox=false);
  
  /// display a G4Step, showing the enclosing volume (if a box).
  void display(G4TouchableHistory* touched, idents::VolumeIdentifier id, 
               const HepPoint3D& entry, const HepPoint3D& exit);
  
  /// @return Return the display manager ponter
  DisplayManager* displayMgr(){return m_display;} 
  
  typedef std::map<idents::VolumeIdentifier, unsigned int> DetectorList;

  /// A map to keep track of hit detectors for display
  DetectorList m_detectorList;

  /// The pointer to the IdataProviderSvc
  IDataProviderSvc* m_esv;
  
 private:
  /// The identifier indicized map of physical volume
  DetectorConstruction::IdMap* m_idMap;
  DisplayManager* m_display;
};
#endif
