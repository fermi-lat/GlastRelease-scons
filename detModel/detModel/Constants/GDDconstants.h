#ifndef GDDCONSTANTS_H
#define GDDCONSTANTS_H

/**
 * \todo All the constants package has to be implemented
 * @author D.Favretto & R.Giannitrapani 
 */
#include <vector.h>
#include <string.h>

class GDDconstCategory;

class GDDconstants{
 public:

  GDDconstants():majVersion(""),minVersion(""){;}
  ~GDDconstants();//use purge.h

  void setVersion(string pmajVersion,string pminVersion)
    {majVersion=pmajVersion;minVersion=pminVersion;}  
  string getVersion()const{return concat(majversion,minVersion);}

  void addConstantCategory(GDDconstCategory* c){constCategories.push_back(c);}
  int getConstantCategoryNumber()const{return constCategories.size();}  
  vector<GDDconst*> getConstantCategories()const{return constsCategories;}

 private:    
  /** @link aggregation */
  string majVersion;
  string minVersion;

  /** @link aggregation */
  vector < GDDconstCategory * > constCategories;
};
#endif //GDDCONSTANTS_H



