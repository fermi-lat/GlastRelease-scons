#ifndef G4SECTIONSVISITOR_H
#define G4SECTIONSVISITOR_H
#include "detModel/Management/SectionsVisitor.h"
#include <fstream>
#include <vector>
#include <map>

class G4Box;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VisAttributes;
class Gdd;
class G4VSensitiveDetector;
/*
 * This is a concrete implementation of a sectionsVisitor that produces
 * a VRML file with the geometry. 
 */
class G4SectionsVisitor : public detModel::SectionsVisitor {

 public:

  G4SectionsVisitor();
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

  /// This method build the colors for the VRML file
  void makeColor();

  void setOpacity(std::string, float);

  void setActualMother(G4LogicalVolume* pactualMother){actualMother = pactualMother;};

  void setActualVolume(std::string pvol){actualVolume = pvol;};
  std::string getActualVolume(){return actualVolume;};

  G4LogicalVolume* getLogicalByName(std::string name);

  /// This map holds the opacity information of the material colors
  std::map <std::string, float> opacityMap;

  double compX;
  double compY;
  double compZ;
  double trans;
  double halfPrec;

  std::vector <G4Box*> g4Boxes;
  std::vector <G4LogicalVolume*> g4Logicals;  
  std::vector <G4VPhysicalVolume*> g4Physicals;  
  std::map <std::string, G4VisAttributes*> g4VisAttributes;  

  std::map <G4VPhysicalVolume*, std::string> g4Identifiers;  

  std::string actualVolume;

  G4VSensitiveDetector* sensibleDet;
  G4LogicalVolume* actualMother;
  G4VPhysicalVolume* worldphys;
  G4LogicalVolume* worldlog;
};
#endif //G4SECTIONSVISITOR_H







