#ifndef GDDCONST_H
#define GDDCONST_H

#include <string>

/**
 * @author D.Favretto & R.Giannitrapani 
 */

class GDDconstantsVisitor;

enum constType{i,f,d,s};///int float double string
enum constMeaning{length,angle,coun,mat,none};

class GDDconst {
 public:

  virtual ~GDDconst(){}

  void setName(string pname){name=pname;}
  void setConstType(constType pct){ct=pct;}
  void setConstMeaning(string puType);
  void setNote(string pnote){note=pnote;}


  string getName()const{return name;}
  constType getConstType()const{return ct;}
  constMeaning getConstMeaning()const{return uType;}
  string getNote()const{return note;}
  string getUnitLength()const{return unitLength;}
  string getUnitAngle()const{return unitAngle;}

  /**
   * This is the standard Accept method for the visitor mechanism
   */
  void Accept(GDDconstantsVisitor*);
  /**
   *  This method is the non recursive Accept for the visitor pattern
   */
  void AcceptNotRec(GDDconstantsVisitor*);

 private:
  string name;
  constType ct;
  constMeaning uType;
  string note;
  string unitLength;
  string unitAngle;
 protected:
  GDDconst():name(""),note(""),unitLength("mm"),unitAngle("deg"){;}
};
#endif //GDDCONST_H

























