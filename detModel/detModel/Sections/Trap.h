#ifndef TRAP_H
#define TRAP_H
#include <string>
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Shape.h"

namespace detModel{

  /**
   * This class represent a solid with two parallel and congruent trapezoidal
   * faces
   * @author J. Bogart
   */
  class Trap : public Shape {
  public:
  
    Trap(std::string pname):Shape(pname),m_x1(0),m_x2(0),m_xDiff(0), 
                            m_y(0),m_z(0){;}
    Trap(std::string pname, double x1, double x2, double xDiff,
         double y, double z):
      Shape(pname),m_x1(x1), m_x2(x2), m_xDiff(xDiff),  m_y(y), m_z(z) {;}

    ~Trap(){;}

    void setX1(double x){m_x1 = x;}
    void setX2(double x){m_x2 = x;}
    void setXDiff(double xd){m_xDiff = xd;}
    void setY(double y){m_y = y;}
    void setZ(double z){m_z = z;}

    double getX1()const{return m_x1;}
    double getX2()const{return m_x2;}
    double getXDiff()const{return m_xDiff;}
    double getY()const{return m_y;}
    double getZ()const{return m_z;}

    /**
     * This method build the bounding box information
     */
    virtual void buildBB();
    
    /**
     * This is the standard Accept method for the visitor mechanism
     */
    void Accept(SectionsVisitor* v);
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    void AcceptNotRec(SectionsVisitor* v){v->visitTrap(this);};

  private:    
    double m_x1; ///default 0; length of bottom
    double m_x2; ///default 0; length of top
    double m_xDiff; /// default 0; difference in x of centers of top and bottom
    double m_y; ///default 0
    double m_z; ///default 0

  };
}
#endif //TRAP_H







