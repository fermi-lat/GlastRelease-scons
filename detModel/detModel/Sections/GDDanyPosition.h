#ifndef GDDANYPOSITION_H
#define GDDANYPOSITION_H
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDidField.h"
#include "detModel/Utilities/GDDpurge.h"
#include <vector>

class GDDsectionsVisitor;

enum positiontype{mposX,mposY,mposZ,mposR,mposPhi,posRPhi,posXYZ};

/**
 * This class represents a generic position
 * @author D.Favretto & R.Giannitrapani
 */
class GDDanyPosition {
public:
  virtual ~GDDanyPosition(){purge(idFields);}
  /**
   *
   */ 
  virtual void Accept(GDDsectionsVisitor*){};
  /**
   *
   */ 
  virtual void AcceptNotRec(GDDsectionsVisitor*){};
  /**
   *
   */ 
  void setVolumeRef(string pname){volumeRef=pname;}
  /**
   *
   */ 
  void setVolume(GDDvolume* pvolume){volume=pvolume;}
  /**
   * This method set the units for the lenghts of the shape
   */  
  void setUnitLength(char* punitLength){
    if (punitLength== "mm")
      unitLength = GDDmm;
    else if (punitLength== "cm")
      unitLength = GDDcm;
    else unitLength = GDDm;
  };

  vector <GDDidField*> getIdFields(){return idFields;}
  /**
   * This method set the units for the angles of the shape
   */  
  void setUnitAngle(char* punitAngle){
    if (punitAngle== "deg")
      unitAngle = GDDdeg;
    else unitAngle = GDDmrad;
  };
   /**
   *
   */ 
  void setXrot(double pxrot){xrot=pxrot;}
  /**
   *
   */ 
  void setYrot(double pyrot){yrot=pyrot;}
  /**
   *
   */ 
  void setZrot(double pzrot){zrot=pzrot;}
  /**
   *
   */ 
  void addIdField(GDDidField* pidField){idFields.push_back(pidField);}
    /**
   *
   */ 
  double getXrot()const{return xrot;}
  /**
   *
   */ 
  double getYrot()const{return yrot;}
  /**
   *
   */ 
  double getZrot()const{return zrot;}
  /**
   *
   */ 
  virtual double getBBX(){return 0;}
  /**
   *
   */ 
  virtual double getBBY(){return 0;}
  /**
   *
   */ 
  virtual double getBBZ(){return 0;}
  /**
   *
   */   
  void setS(double ps){S=ps;}
  /**
   *
   */ 
  double getS()const{return S;}
  /**
   *
   */ 
  string getVolumeRef()const{return volumeRef;}
  /**
   *
   */ 
  GDDvolume* getVolume()const{return volume;}
  /**
   *
   */ 
  unitLengthtype getUnitLength()const{return unitLength;}
  /**
   *
   */ 
  unitAngletype getUnitAngle()const{return unitAngle;}

  positiontype getPosType(){return typePosition;};

  /** @link aggregation */
  vector <GDDidField*> idFields;

 protected:
  GDDanyPosition(positiontype ptypePosition):
    typePosition(ptypePosition),xrot(0),yrot(0),zrot(0),S(0),volume(0),
    volumeRef(""){;}
  
 private:
  positiontype typePosition;
  double xrot; ///default 0
  double yrot; ///default 0
  double zrot; ///default 0
  double S;  ///default 0
  unitLengthtype unitLength; ///dufault mm
  unitAngletype unitAngle; ///default deg

  
  /** @link aggregation
   *   @supplierCardinality 1 */
  GDDvolume * volume;
  string volumeRef;
  
};
#endif //GDDANYPOSITION_H












