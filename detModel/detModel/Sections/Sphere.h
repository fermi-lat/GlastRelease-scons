#ifndef SPHERE_H
#define SPHERE_H
#include <string>
#include "detModel/Sections/Shape.h"
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{

  /**
   * This class represent a simple spheircal portion of space
   * 
   * @author R.Giannitrapani 
   */
  class Sphere :public Shape {
  
  public: 
    Sphere(std::string pname):Shape(pname),m_rin(0),
                              m_rout(0),m_phimin(0),
                              m_phimax(0),m_thetamin(0),
                              m_thetamax(0){;}
    Sphere(std::string pname, double rin, 
        double rout, double phimin, 
        double phimax, double thetamin, 
        double thetamax):Shape(pname),m_rin(rin),m_rout(rout),
                         m_phimin(phimin),m_phimax(phimax),
                         m_thetamin(thetamin),m_thetamax(thetamax){;};


    double getRin(){return m_rin;};
    void setRin(double rin){m_rin = rin;};

    double getRout(){return m_rout;};
    void setRout(double rout){m_rout = rout;};

    double getPhiMin(){return m_phimin;};
    void setPhiMin(double phimin){m_phimin = phimin;};

    double getPhiMax(){return m_phimax;};
    void setPhiMax(double phimax){m_phimax = phimax;};

    double getThetaMin(){return m_thetamin;};
    void setThetaMin(double thetamin){m_thetamin = thetamin;};

    double getThetaMax(){return m_thetamax;};
    void setThetaMax(double thetamax){m_thetamax = thetamax;};

    /**
     * This method build the bounding box information
     */
    virtual void buildBB();
    
    /**
     * This is the standard Accept method for the visitor mechanism
     */
    void Accept(SectionsVisitor* v){v->visitSphere(this);};
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    void AcceptNotRec(SectionsVisitor* v){v->visitSphere(this);};

  private:
    double m_rin;
    double m_rout;
    double m_phimin;
    double m_phimax;  
    double m_thetamin;
    double m_thetamax;
  };
}
#endif //SPHERE_H
