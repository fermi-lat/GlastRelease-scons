// $Header$
#ifndef G4SECTIONSVISITOR_H
#define G4SECTIONSVISITOR_H
#include "detModel/Management/SectionsVisitor.h"
#include <fstream>
#include <vector>
#include <map>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VisAttributes;
class Gdd;
//class G4VSensitiveDetector;
#include "idents/VolumeIdentifier.h"

namespace detModel { class Position; }
/*
 * This is a concrete implementation of a sectionsVisitor that produces
 * the Geant4 geometry. 
 */
class G4SectionsVisitor : public detModel::SectionsVisitor {

 public:


     //! Id's are represented at a vector of unsigned ints
     typedef std::map<const G4VPhysicalVolume*, idents::VolumeIdentifier> IdMap;

   //! ctor, with top volume and pointer to map of ids to fill
   G4SectionsVisitor(std::string topvol, IdMap* idmap);
  virtual ~G4SectionsVisitor();
  
  /**
   * This is the visitor for the GDDsectionsContainer 
   */
  virtual void visitGdd(detModel::Gdd*);
  
  /**
   * This is the visitor for the GDDsection 
   */
  virtual void visitSection(detModel::Section*);

  /**
   * This is the visitor for the GDDstack 
   */
  virtual void visitEnsemble(detModel::Ensemble*);

  /**
   * This is the visitor for the GDDbox 
   */
  virtual void visitBox(detModel::Box*);

  /**
   * This is the visitor for the Tube 
   */
  virtual void visitTube(detModel::Tube*);

  /**
   * This is the visitor for the GDDposXYZ 
   */
  virtual void visitPosXYZ(detModel::PosXYZ*);

  /**
   * This is the visitor for the GDDaxisMPos 
   */
  virtual void visitAxisMPos(detModel::AxisMPos*);

  /**
   * This is the visitor for the GDDaxisMPos 
   */
  virtual void visitIdField(detModel::IdField*);

  /**
   * This is the visitor for the GDDseg 
   */
  virtual void visitSeg(detModel::Seg*);

  typedef std::vector< G4LogicalVolume*> Logicals;
  const Logicals & getLogicals()const { return g4Logicals; }

  /// return iterators for the logicals list
  Logicals::const_iterator begin()const{ return g4Logicals.begin(); }
  Logicals::const_iterator end()  const{ return g4Logicals.end(); }
 

  /// make a summary  of the volumes on the stream
  void summary(std::ostream out)const;
  
  ///  return pointer to the World physical volume
  G4VPhysicalVolume* getWorld()const { return m_worldphys; }

 private:

  /// This method build the colors for the VRML file
  void makeColor();

  void setOpacity(std::string, float);

  void setActualMother(G4LogicalVolume* pactualMother){actualMother = pactualMother;};

  void setActualVolume(std::string pvol){actualVolume = pvol;};
  std::string getActualVolume(){return actualVolume;};

  G4LogicalVolume* getLogicalByName(std::string name);


  /// This map holds the opacity information of the material colors
  std::map <std::string, float> opacityMap;


  std::vector <G4LogicalVolume*> g4Logicals;  
  std::vector <G4VPhysicalVolume*> g4Physicals;  
  std::map <std::string, G4VisAttributes*> g4VisAttributes;  

  std::string actualVolume;

  G4LogicalVolume* actualMother;
  G4VPhysicalVolume* m_worldphys;

   //! private function to manage the identifiers
   //! @param pos a Position object, in practice either a PosXYZ or an AxisMpos.
   //! @param i the index for an element of a stack (only needed for AxixMpos)
   void processIds(/*const*/  detModel::Position* pos, unsigned int i=0); 

   //! a little map to count the number of physical volumes for each logical one
  std::map<std::string, int> m_physicalsPerLogical;

  // pointer to the clients map that we fill
  IdMap* m_idMap;

};
#endif //G4SECTIONSVISITOR_H








