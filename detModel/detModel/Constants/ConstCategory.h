#ifndef CONSTCATEGORY_H
#define CONSTCATEGORY_H
#include<vector>
#include<string>
#include "detModel/Utilities/Purge.h"
#include "detModel/Constants/Const.h"

namespace detModel{

  class ConstantsVisitor;

  /**
   * @author D.Favretto & R.Giannitrapani 
   */
  class ConstCategory {
  public:

    ConstCategory():categoryName(""),overview(""),primary(true),save(false){;} 

    virtual ~ConstCategory(){purge(consts);}

    void setName(std::string name){categoryName=name;}
    std::string getName()const{return categoryName;}

    void setOverview(std::string poverview){overview=poverview;}
    std::string getOverview()const{return overview;}

    void addConstant(Const* c){consts.push_back(c);}
    int getConstantNumber()const{return consts.size();}  
    std::vector<Const*>getConsts()const{return consts;}

    bool getPrimary()const{return primary;}
    void setPrimary(bool p){primary=p;}

    bool getSave()const{return save;}
    void setSave(bool p){save=p;}

    /**
     * This is the standard Accept method for the visitor mechanism
     */
    void Accept(ConstantsVisitor*);
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    void AcceptNotRec(ConstantsVisitor*);

  private:
    std::string categoryName;
    std::string overview;
    /** @link aggregation */
    std::vector < Const * > consts;
    bool primary;/// this attribute  expres if the category is primary=1 or derived=0
    bool save;

  };
}
#endif //CONSTCATEGORY_H








