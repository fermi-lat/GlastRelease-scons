#ifndef IDFIELD_H
#define IDFIELD_H
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{

  /**
   * @author D.Favretto & R.Giannitrapani
   */
  class IdField {
  public:

    IdField():value(0),step(0){}

    virtual ~IdField(){};

    /**
     *
     */ 
    virtual void Accept(SectionsVisitor* v){v->visitIdField(this);}
    /**
     *
     */ 
    void setValue(int pvalue){value=pvalue;}
    /**
     *
     */ 
    void setStep(int pstep){step=pstep;}
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
    std::string getName()const;
  
  private:
    std::string name;
    int value; ///default 0 NMTOKEN
    int step; ///default 0 NMTOKEN
  };
}
#endif //IDFIELD_H





