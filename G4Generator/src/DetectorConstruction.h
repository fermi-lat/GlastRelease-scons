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

/** callec by the RunManager to construct:
- full detector geometry
- define materials
- identify sensitive detectors, set up callbacks
- create a map of partical ids and physical volumes
*/
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //! @param gds A pointer to the abstract interface of the GlastDetSvc

  DetectorConstruction(IGlastDetSvc* gds, IDataProviderSvc* esv);
  ~DetectorConstruction();
  
  //! actual call to construct things.
  G4VPhysicalVolume* Construct();
  
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier > IdMap;
  
  //! access to the map of physical volume /id pairs.
  IdMap* idMap(){return &m_idMap;}

  const std::string & topVolumeName()const {return m_topvol;}
  
private:
  /// A pointer to the GlastDetSvc
  IGlastDetSvc* m_gsv;
  
  std::string m_topvol;
  
  //! map with integer id's
  IdMap m_idMap;
  
  // delegate creation of detectors.
  PosDetectorManager* m_posDet;

  // delegate creation of detectors.
  IntDetectorManager* m_intDet;
};

#endif

