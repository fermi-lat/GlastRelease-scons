#ifndef GDDCONSTANTS_H
#define GDDCONSTANTS_H

/**
 * \todo All the constants package has to be implemented
 * @author D.Favretto & R.Giannitrapani 
 */
#include <vector.h>
#include <string.h>
#include "detModel/GDDpurge.h"
#include "detModel/Constants/GDDconstCategory.h"

class GDDconstantsVisitor;

class GDDconstants{
 public:
  
  GDDconstants():majVersion(""),minVersion(""){;}  
  virtual ~GDDconstants(){purge(constCategories);}
  
  void setVersion(string pmajVersion,string pminVersion)
    {majVersion=pmajVersion ; minVersion=pminVersion;}  
  string getVersion()const{return majVersion+"."+minVersion;}
  
  void addConstantCategory(GDDconstCategory* c){constCategories.push_back(c);};
  int getConstantCategoryNumber()const{return constCategories.size();};
  vector< GDDconstCategory* > getCategories()const{return constCategories;};

  /**
   * This is the standard Accept method for the visitor mechanism
   */
  void Accept(GDDconstantsVisitor*);
  /**
   *  This method is the non recursive Accept for the visitor pattern
   */
  void AcceptNotRec(GDDconstantsVisitor*);
  
 private:    
  /** @link aggregation */
  string majVersion;
  string minVersion;
  
  /** @link aggregation */
  vector < GDDconstCategory * > constCategories;

};
#endif //GDDCONSTANTS_H







