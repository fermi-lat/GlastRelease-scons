#ifndef GDDIDFIELD_H
#define GDDIDFIELD_H
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */

class GDDidField {
public:

  GDDidField():value(0),step(0){}

  virtual ~GDDidField(){};

  /**
   *
   */ 
  virtual void Accept(GDDsectionsVisitor* v){v->visitIdField(this);}
  /**
   *
   */ 
  void setValue(double pvalue){value=pvalue;}
  /**
   *
   */ 
  void setStep(double pstep){step=pstep;}
  /**
   *
   */ 
  void setName(std::string pname);
  /**
   *
   */ 
  double getValue()const{return value;}
  /**
   *
   */ 
  double getStep()const{return step;}
  /**
   *
   */ 
  std::string getName();
  
 private:
  std::string name;
  double value; ///default 0 NMTOKEN
  double step; ///default 0 NMTOKEN
};
#endif //GDDIDFIELD_H





