#ifndef CONST_H
#define CONST_H
#include <string>

namespace detModel{

  /**
   * @author D.Favretto & R.Giannitrapani 
   */

  class ConstantsVisitor;


  class Const {
  public:
    enum constMeaning{length,angle,coun,mat,energy,none};

    Const():name(""),note(""),unitLength("mm"),unitAngle("deg"),
            unitEnergy("MeV"){;}

    virtual ~Const(){}

    /**
     * This is the standard Accept method for the visitor mechanism
     */
    virtual void Accept(ConstantsVisitor*) =0;
    /**
     *  This method is the non recursive Accept for the visitor pattern
     */
    virtual void AcceptNotRec(ConstantsVisitor*) = 0;
  
    void setName(std::string pname){name=pname;}
    void setConstMeaning(std::string puType);
    void setNote(std::string pnote){note=pnote;}

    std::string getName()const{return name;}
    constMeaning getConstMeaning()const{return uType;}
    std::string getNote()const{return note;}
    std::string getUnitLength()const{return unitLength;}
    std::string getUnitAngle()const{return unitAngle;}
    std::string getUnitEnergy()const {return unitEnergy;}

  private:
    std::string name;
    constMeaning uType;
    std::string note;
    std::string unitLength;
    std::string unitAngle;
    std::string unitEnergy;
  };
}
#endif //CONST_H
