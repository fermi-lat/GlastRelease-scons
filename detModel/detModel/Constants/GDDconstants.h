#ifndef GDDCONSTANTS_H
#define GDDCONSTANTS_H

/**
 * \todo All the constants package has to be implemented
 * @author D.Favretto & R.Giannitrapani 
 */
#include <vector>
#include <string>
#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Constants/GDDconstCategory.h"

class GDDconstantsVisitor;

class GDDconstants{
 public:
  
  GDDconstants():majVersion(""),minVersion(""){;}  
  virtual ~GDDconstants(){purge(constCategories);}
  
  void setVersion(std::string pmajVersion,std::string pminVersion)
    {majVersion=pmajVersion ; minVersion=pminVersion;}  
  std::string getVersion()const{return majVersion+"."+minVersion;}
  
  void addConstantCategory(GDDconstCategory* c){constCategories.push_back(c);};
  int getConstantCategoryNumber()const{return constCategories.size();};
  std::vector< GDDconstCategory* > getCategories()const{return constCategories;};

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
  std::string majVersion;
  std::string minVersion;
  
  /** @link aggregation */
  std::vector < GDDconstCategory * > constCategories;

};
#endif //GDDCONSTANTS_H







