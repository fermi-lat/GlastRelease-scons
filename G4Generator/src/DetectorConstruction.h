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

class GlastDetectorManager;

/** callec by the RunManager to construct:
- full detector geometry
- define materials
- identify sensitive detectors, set up callbacks
- create a map of partical ids and physical volumes
*/
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  //! @param topvol Name of the topvolume, e.g. "LAT"
  //! @param mode   Visitor mode
  DetectorConstruction(std::string topvol, std::string mode);
  ~DetectorConstruction();
  
  //! actual call to construct things.
  G4VPhysicalVolume* Construct();
  
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier > IdMap;
  
  //! access to the map of physical volume /id pairs.
  IdMap* idMap(){return &m_idMap;}

  const std::string & topVolumeName()const {return m_topvol;}
  
private:
  
  std::string m_topvol;
  
  //! map with integer id's
  IdMap m_idMap;
  
  // delegate creation of detectors.
  GlastDetectorManager* m_glastdet;
};

#endif

