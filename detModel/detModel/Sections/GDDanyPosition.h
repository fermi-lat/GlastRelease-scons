#ifndef GDDANYPOSITION_H
#define GDDANYPOSITION_H
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDidField.h"


class GDDsectionsVisitor;

enum positiontype{mposX,mposY,mposZ,mposR,mposPhi,posRPhi,posXYZ};

/**
 * This class represents a generic position
 * @author D.Favretto & R.Giannitrapani
 */
class GDDanyPosition {
public:
  virtual ~GDDanyPosition(){if (idField) delete idField;}
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
      unitLength = mm;
    else if (punitLength== "cm")
      unitLength = cm;
    else unitLength = m;
  };

  GDDidField* getIdField(){return idField;}
  /**
   * This method set the units for the angles of the shape
   */  
  void setUnitAngle(char* punitAngle){
    if (punitAngle== "deg")
      unitAngle = deg;
    else unitAngle = mrad;
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
  void setIdField(GDDidField* pidField){idField=pidField;}
    /**
   *
   */ 
  double getXrot()const{return xrot;}
  /**
   *
   */ 
  double getYrot()const{return xrot;}
  /**
   *
   */ 
  double getZrot()const{return xrot;}
  /**
   *
   */ 
  virtual double getBBX(){;}
  /**
   *
   */ 
  virtual double getBBY(){;}
  /**
   *
   */ 
  virtual double getBBZ(){;}
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

 protected:
  GDDanyPosition(positiontype ptypePosition):
    typePosition(ptypePosition),xrot(0),yrot(0),zrot(0),S(0),volume(0),
    volumeRef(""), idField(0){;}
  
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
  
  /** @link aggregation */
  GDDidField* idField;
};
#endif //GDDANYPOSITION_H












