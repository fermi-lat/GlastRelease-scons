// $Header$

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

/**
This class instantiates the G4 geometry, from detModel
*/
class G4Geometry : public IGeometry
{
public:

  //! Id's are represented at a vector of unsigned ints
  typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier> IdMap;

  G4Geometry();
  G4Geometry(PosDetectorManager* pdm,IntDetectorManager* idm,
	     IdMap* idmap):m_pdm(pdm),m_idm(idm),m_idMap(idmap),
      m_worldPhys(0),m_replica(0),m_replicaMother(0){};

  ~G4Geometry();

  /**     
   * @param s type of the shape
   * @param id vector of unsigned ints (maybe null)
   * @param name
   * @param material
   * @param params vector with the six transformation parameters, followed by 3 or so dimensions
  */
  virtual void pushShape(ShapeType s, const UintVector& id, std::string name, 
			 std::string material, const DoubleVector& params, VolumeType type);
  
  //* called to signal end of nesting */
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

  G4LogicalVolume* actualMother()const{
      if (m_actualMother.size()) return m_actualMother.back();
        else return 0;}
  
  /// Return the number of physical volumes created
  unsigned int getPhysicalNumber(){return m_physicals.size();};

 private:
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

  /// pointer to the clients map that we fill
  IdMap* m_idMap;

  /// the detector managers
  PosDetectorManager* m_pdm;    
  IntDetectorManager* m_idm;    

  /// a flag for optimization mechanism of the G4 geometry
  bool m_replica;

  /// a logical for optimizaion mechanism of the G4 geometry
  G4LogicalVolume* m_replicaMother;
};
#endif
