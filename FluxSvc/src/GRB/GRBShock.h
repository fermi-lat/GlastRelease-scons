/**
 * GRBShock: Description of a Shock between 2 Shells
 *
 * /author Nicola Omodei nicola.omodei@pi.infn.it 
 * /author Johann Cohen Tanugi johann.cohen@pi.infn.it
 *
 * \b revision:
 *   02/15/2002
 *
 */

#include "GRBShell.h"
#include "FluxSvc/mainpage.h"

#ifndef GRBSHOCK_H
#define GRBSHOCK_H 1


class GRBShock 
{
  
 public:
  //Constructors and destructors:
  GRBShock(GRBShell*, GRBShell*);
  ~GRBShock() { }
  
  //Accessors:
  inline double time() {return m_time;}
  inline double tobs() {return m_tobs;}
  inline double Radius() {return m_radius;}
  inline double Eint() {return m_Eint;}
  inline double gf() {return m_gf;}
  inline double VolCom() {return m_VolCom;}
  inline double npart() {return m_npart;}
  inline double Beq() {return m_Beq;}
  inline double Sum() {return m_sum;}
  
  //set functions:
  inline void setTime(double value) {m_time = value;m_tobs=m_time-m_radius/cst::c;}
  inline void setTobs(double value) {m_tobs = value;}
  inline void setSum(double value)  {m_sum=value;}
  //high level functions:
  void Write();
  double Esyn(double gammae);
  double Eic(double gamme);
  double OpticalDepht();
  double tsyn(double gammae);
  double fred(double ee, double tt);
  double Fsyn(double ee, double tt);
  double Fic(double ee, double tt);
  
 private:
  GRBShell* Sh1;
  GRBShell* Sh2;
  double m_time;
  double m_tobs;
  double m_radius;
  double m_mass;
  double m_Eint;
  double m_gsh;
  double m_VolCom;
  double m_npart;
  
  double m_gemin;
  double m_gemax;
  double m_gecoo;
  
  double m_tsyn;
  double m_riset;
  double m_sum;
  
  double m_gf;
  double m_Beq;
  
};
#endif
