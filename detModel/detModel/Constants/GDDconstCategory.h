#ifndef GDDDERCATEGORY_H
#define GDDDERCATEGORY_H

#include<vector.h>
#include<string.h>

#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Constants/GDDconst.h"

class GDDconstantsVisitor;

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDconstCategory {
 public:

  GDDconstCategory():categoryName(""),overview(""),primary(true),save(false){;} 

  virtual ~GDDconstCategory(){purge(consts);}

  void setName(string name){categoryName=name;}
  string getName()const{return categoryName;}

  void setOverview(string poverview){overview=poverview;}
  string getOverview()const{return overview;}

  void addConstant(GDDconst* c){consts.push_back(c);}
  int getConstantNumber()const{return consts.size();}  
  vector<GDDconst*>getConsts()const{return consts;}

  bool getPrimary()const{return primary;}
  void setPrimary(bool p){primary=p;}

  bool getSave()const{return save;}
  void setSave(bool p){save=p;}

  /**
   * This is the standard Accept method for the visitor mechanism
   */
  void Accept(GDDconstantsVisitor*);
  /**
   *  This method is the non recursive Accept for the visitor pattern
   */
  void AcceptNotRec(GDDconstantsVisitor*);

 private:
  string categoryName;
  string overview;
  /** @link aggregation */
  vector < GDDconst * > consts;
  bool primary;/// this attribute  expres if the category is primary=1 or derived=0
  bool save;

};
#endif //GDDDERCATEGORY_H








