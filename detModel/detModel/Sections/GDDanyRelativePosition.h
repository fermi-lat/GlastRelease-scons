#ifndef GDDANYRELATIVEPOSITION_H
#define GDDANYRELATIVEPOSITION_H
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDidField.h"

class GDDsectionsVisitor;

enum axisPosType{AxisPos,AxisMpos};

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDanyRelativePosition {
public:
   ~GDDanyRelativePosition(){if (idField)delete idField;}
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
   *
   */ 
   void setIdField(GDDidField* pidField){idField=pidField;}
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
  void setDx(double pdx){dx=pdx;}
  /**
   *
   */ 
  void setDy(double pdy){dy=pdy;}
  /**
   *
   */ 
  void setDz(double pdz){dz=pdz;}
  /**
   *
   */ 
  void setRotation(double protation){rotation=protation;}
  /**
   *
   */ 
  virtual double getBBX(){};
  /**
   *
   */ 
  virtual double getBBY(){};
  /**
   *
   */ 
  virtual double getBBZ(){};
  /**
   *
   */ 
  
  double getDx(){return dx;}
  /**
   *
   */ 
  double getDy(){return dy;}
  /**
   *
   */ 
  double getDz(){return dz;}
  /**
   *
   */ 
  double getRotation(){return rotation;}
  /**
   *
   */ 
  stacktype getPosDir(){return posDir;};
  /**
   *
   */ 
  void setPosDir(stacktype pst){posDir = pst;};
  /**
   *
   */ 
  unitLengthtype getUnitLength()const{return unitLength;}
  /**
   *
   */ 
  unitAngletype getUnitAngle()const{return unitAngle;}
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
  axisPosType getPosType()const{return positionType;}
  /**
   *
   */ 
  GDDidField* getIdField()const{return idField;}

 protected:
  GDDanyRelativePosition(axisPosType ppositionType):
    positionType(ppositionType),dx(0),dy(0),dz(0),
   rotation(0), idField(0){}

 private:    
  
  /** @link aggregation 
   * @clientCardinality 1
   * @supplierCardinality 1*/
  GDDvolume* volume;
  string volumeRef;
  axisPosType positionType;
  double dx; ///dx default 0
  double dy; ///dy default 0
  double dz; ///dz default 0
  double rotation; ///default 0
  unitLengthtype unitLength; ///dufault mm
  unitAngletype unitAngle; ///default deg
  stacktype posDir;
  
  /** @link aggregation */
  GDDidField * idField;
};
#endif //GDDANYRELATIVEPOSITION_H




