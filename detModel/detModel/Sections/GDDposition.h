#ifndef GDDPOSITION_H
#define GDDPOSITION_H
#include "detModel/Sections/GDDboundable.h"
#include "detModel/Sections/GDDunits.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Utilities/GDDpurge.h"
#include <vector>

class GDDsectionsVisitor;
class GDDunits;
class GDDvolume;
/**
 * This class represents a generic positioned volume
 * @author R.Giannitrapani
 */
class GDDposition : public GDDboundable{
 public:
  
  /// A standard constructor
  GDDposition():volume(0),volumeRef(""){units = new GDDunits;}

  /// In the destructor we need to purge the idFields vector and the units
  virtual ~GDDposition(){purge(idFields); delete units;}
  virtual void Accept(GDDsectionsVisitor*) = 0;
  virtual void AcceptNotRec(GDDsectionsVisitor*) =0;

  /**
   * This methods sets and gets the volume and volumeref
   */ 
  void setVolumeRef(std::string pname){volumeRef=pname;}
  void setVolume(GDDvolume* pvolume){volume=pvolume;}
  std::string getVolumeRef()const{return volumeRef;}
  GDDvolume* getVolume()const{return volume;}

  /**
   *  This methods add an idField and retrive the vector of idField
   */ 
  void addIdField(GDDidField* pidField){idFields.push_back(pidField);}
  std::vector <GDDidField*> getIdFields(){return idFields;}

 private:
  /** @link aggregation */
  std::vector <GDDidField*> idFields;
  GDDunits* units;
  /** @link aggregation
   *   @supplierCardinality 1 */
  GDDvolume * volume;
  std::string volumeRef;
  
};
#endif //GDDPOSITION_H












