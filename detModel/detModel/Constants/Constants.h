#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <vector>
#include <string>
#include "detModel/Utilities/Purge.h"
#include "detModel/Constants/ConstCategory.h"

namespace detModel{

  class ConstantsVisitor;

  /**
   * 
   * @author D.Favretto & R.Giannitrapani 
   */
  class Constants{
  public:
  
    Constants():majVersion(""),minVersion(""){;}  
    virtual ~Constants(){purge(constCategories);}
  
    void setVersion(std::string pmajVersion,std::string pminVersion)
      {majVersion=pmajVersion ; minVersion=pminVersion;}  
    std::string getVersion()const{return majVersion+"."+minVersion;}
  
    void addConstantCategory(ConstCategory* c){constCategories.push_back(c);};
    int getConstantCategoryNumber()const{return constCategories.size();};
    std::vector< ConstCategory* > getCategories()const{return constCategories;};

    /**
     * This is the standard Accept method for the visitor mechanism
     */
    void Accept(ConstantsVisitor*);
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    void AcceptNotRec(ConstantsVisitor*);
  
  private:    
    /** @link aggregation */
    std::string majVersion;
    std::string minVersion;
  
    /** @link aggregation */
    std::vector < ConstCategory * > constCategories;

  };
}
#endif //GDDCONSTANTS_H







