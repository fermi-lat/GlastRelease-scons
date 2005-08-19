#ifndef G4GEOMETRY_h
#define G4GEOMETRY_h

#include <vector>
#include <map>
#include <fstream>

#include "idents/VolumeIdentifier.h"
#include "GlastSvc/GlastDetSvc/IGeometry.h"

namespace ident{class VolumeIdentifier;}
class G4LogicalVolume;
class G4VPhysicalVolume;
class PosDetectorManager;
class IntDetectorManager;

class DetectorConstruction;

/**
 * @class G4Geometry
 * 
 * @brief Utility class to instantiate the Geant4 geometry 
 *
 * This class instanciates the Geant4 geometry from the xml representation; this
 * is done via the visitors mechanism provided by detModel. To hide the concrete
 * implementation of detModel, this is done with methods provided by the
 * abstracts GlastSvc interfaces. The G4Geometry implements one of these
 * interfaces.
 *
 * @author R.Giannitrapani
 * 
 * $Header$
 */
class G4Geometry : public IGeometry
{
 public:

  //! Id's are represented as a vector of unsigned ints
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier> IdMap;

  G4Geometry(std::string mode="propagate");
  G4Geometry(PosDetectorManager* pdm,IntDetectorManager* idm,
             IdMap* idmap, std::string mode="propagate",
	     DetectorConstruction* det=0) : 
    m_pdm(pdm),m_idm(idm),m_idMap(idmap), 
    m_worldPhys(0),m_replica(0),m_replicaMother(0), m_mode(mode),
    m_det(det) {};
  
  ~G4Geometry();

  /**     
   * This method implement the abstract pushShape of IGeometry; it must push a
   * new Geant4 shape in the stack of volumes
   *
   * @param s type of the shape
   * @param id vector of unsigned ints (maybe null)
   * @param name
   * @param material
   * @param params vector with the six transformation parameters, 
   *        followed by 3 or so dimensions
   * @param type kind of volume: simple, stack or composition
   * @param sense Sensitivity type of volume
   * @return tell caller whether to skip subvolumes or not
   */
  virtual IGeometry::VisitorRet
    pushShape(ShapeType s, 
              const UintVector& id, 
              std::string name, 
              std::string material, 
              const DoubleVector& params, 
              VolumeType type,
              SenseType  sense);
  
  /// called to signal end of nesting 
  virtual void popShape();

  /// current ID
  const UintVector& getId()const { return m_idValues; }

  /// return pointer to the mother volume
  G4VPhysicalVolume* getWorld()const{ return m_worldPhys; }

  /// set the mother volume
  void setWorld(G4VPhysicalVolume* world){ m_worldPhys = world; }

  /// Get a physical volume by name
  G4VPhysicalVolume* getPhysicalByName(std::string name);

  /// Push and pop methods for the stack of actual mother volumes
  void pushActualMother(G4LogicalVolume* vol){ m_actualMother.push_back(vol);}
  void popActualMother(){m_actualMother.pop_back();}

  /// Get the mother volume
  G4LogicalVolume* actualMother()const;
  
  /// Returns the number of physical volumes created
  unsigned int getPhysicalNumber()const{return m_physicals.size();};

  /// Returns the number of logical volumes created
  unsigned int getLogicalNumber()const{return m_logicals.size();};

  /// Need a setMode in order to implement getMode for IGeometry interface
  virtual void setMode(std::string pmode) {m_mode = pmode;}
  virtual std::string getMode() {return m_mode;}

 private:
  /// the sensitive detector managers
  PosDetectorManager* m_pdm;    
  IntDetectorManager* m_idm;    

  /// pointer to the clients map that we fill
  IdMap* m_idMap;

  /// The mother volume of the geometry
  G4VPhysicalVolume* m_worldPhys;

  /// This is a stack with the actual mother of the geometry
  std::vector <G4LogicalVolume*> m_actualMother;

  /// This is a map with all the physical volumes indicized by name
  std::vector <G4VPhysicalVolume*> m_physicals;

  /// This is a map with all the logicals indicized by name
  std::map <std::string,G4LogicalVolume*> m_logicals;

  //! vector of ids to describe current object
  UintVector m_idValues;
  
  //! vector of the number of ids for this geometry level: 0,1, or even 2
  UintVector m_idcount;

  /// a flag for optimization mechanism of the G4 geometry
  bool m_replica;

  /// a logical for optimizaion mechanism of the G4 geometry
  G4LogicalVolume* m_replicaMother;

  /// the mode of traversal of the geometry
  std::string m_mode;

  DetectorConstruction* m_det;
};
#endif
