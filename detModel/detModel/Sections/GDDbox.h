#ifndef GDDBOX_H
#define GDDBOX_H
#include <string>
#include "detModel/Management/GDDsectionsVisitor.h"
#include "detModel/Sections/GDDshape.h"
#include "detModel/Sections/GDDseg.h"

class GDDsectionsVisitor;
class GDDboundingBox;

/**
 * This class represent a simple box
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDbox : public GDDshape {
public:
  
    GDDbox(std::string pname):GDDshape(pname),x(0),y(0),z(0),seg(0){;}
    GDDbox(std::string pname, double x, double y, double z):
      GDDshape(pname),x(x), y(y), z(z), seg(0){;}

    ~GDDbox(){if(seg) delete(seg);}

    /**
     * This method set the pointer to seg 
     */
    void setSeg(GDDseg* pseg){seg=pseg;}
    /**
     * This method gives the pointer to seg
     */
    GDDseg* getSeg(){return seg;}
    /**
     * This method set the X dimension 
     */
    void setX(double px){x = px;}
    /**
     * This method set the Y dimension 
     */
    void setY(double py){y = py;}
    /**
     * This method set the Z dimension
     */
    void setZ(double pz){z = pz;}
    /**
     * This method gives the X dimension
     */
    double getX(){return x;}
    /**
     * This method gives the y dimension
     */
    double getY(){return y;}
    /**
     * This method gives the z dimension 
     */
    double getZ(){return z;}

    /**
     * This method build the bounding box information
     */
    virtual void buildBB();
    
    /**
     * This is the standard Accept method for the visitor mechanism
     */
    void Accept(GDDsectionsVisitor* v);
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    void AcceptNotRec(GDDsectionsVisitor* v){v->visitBox(this);};

 private:    
    double x; ///default 0
    double y; ///default 0
    double z; ///default 0

    /** @link aggregation 
     * @supplierCardinality 0..1*/
    GDDseg* seg;
};
#endif //GDDBOX_H







