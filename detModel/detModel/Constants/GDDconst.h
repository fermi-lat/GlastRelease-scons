#ifndef GDDCONST_H
#define GDDCONST_H

#include <string>

/**
 * @author D.Favretto & R.Giannitrapani 
 */

class GDDconstantsVisitor;


class GDDconst {
 public:
  enum constMeaning{length,angle,coun,mat,none};

  GDDconst():name(""),note(""),unitLength("mm"),unitAngle("deg"){;}

  virtual ~GDDconst(){}

  /**
   * This is the standard Accept method for the visitor mechanism
   */
  virtual void Accept(GDDconstantsVisitor*) =0;
  /**
   *  This method is the non recursive Accept for the visitor pattern
   */
  virtual void AcceptNotRec(GDDconstantsVisitor*) = 0;
  
  void setName(std::string pname){name=pname;}
  void setConstMeaning(std::string puType);
  void setNote(std::string pnote){note=pnote;}

  std::string getName()const{return name;}
  constMeaning getConstMeaning()const{return uType;}
  std::string getNote()const{return note;}
  std::string getUnitLength()const{return unitLength;}
  std::string getUnitAngle()const{return unitAngle;}

 private:
  std::string name;
  constMeaning uType;
  std::string note;
  std::string unitLength;
  std::string unitAngle;

};
#endif //GDDCONST_H

























