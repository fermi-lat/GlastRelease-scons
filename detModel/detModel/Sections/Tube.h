#ifndef TUBE_H
#define TUBE_H
#include <string>
#include "detModel/Management/SectionsVisitor.h"
#include "detModel/Sections/Shape.h"


namespace detModel{

  /**
   *  This class represent a tube oriented along the Z axis
   */
  class Tube :public Shape {

  public:
    Tube(std::string pname):Shape(pname),m_rin(0),m_rout(0),m_z(0){;}
    Tube(std::string pname, double rin, double rout, double z):
      Shape(pname),m_rin(rin), m_rout(rout), m_z(z){;}

    double getRin(){return m_rin;};
    void setRin(double rin){m_rin = rin;};
    double getRout(){return m_rout;};
    void setRout(double rout){m_rout = rout;};
    double getZ(){return m_z;};
    void setZ(double z){m_z = z;};

    /**
     * This method build the bounding box information
     */
    virtual void buildBB();
    
    /**
     * This is the standard Accept method for the visitor mechanism
     */
    void Accept(SectionsVisitor* v){v->visitTube(this);};
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    void AcceptNotRec(SectionsVisitor* v){v->visitTube(this);};

  private:
    /// The internal radius
    double m_rin;
    /// The external radius
    double m_rout;
    /// The height
    double m_z;
  };
}
#endif //TUBE_H







