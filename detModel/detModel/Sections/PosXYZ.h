#ifndef POSXYZ_H
#define POSXYZ_H
#include "detModel/Sections/SinglePos.h"
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{

  /**
   * @author D.Favretto & R.Giannitrapani
   */
  class PosXYZ :public SinglePos {
  public:
    /// A standard constructor
    PosXYZ():SinglePos(),x(0),y(0),z(0){;}
    /// A constructor with the coordinate of the positioned volume
    PosXYZ(double x, double y, double z):SinglePos(),
      x(x),y(y),z(z){;}
    
    void Accept(SectionsVisitor* v);
    virtual void AcceptNotRec(SectionsVisitor* v){v->visitPosXYZ(this);};
    
    void setX(double pX){x=pX;}
    void setY(double pY){y=pY;}
    void setZ(double pZ){z=pZ;}
    double getX()const{return x;}
    double getY()const{return y;}
    double getZ()const{return z;}
    
    virtual void buildBB();

  private:    
    double x; 
    double y; 
    double z; 
  
  };
}
#endif //POSXYZ_H



