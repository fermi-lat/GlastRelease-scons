#ifndef GDDCONST_H
#define GDDCONST_H

#include <string.h>

/**
 * @author D.Favretto & R.Giannitrapani 
 */
///\todo rimane da risolvere il problema del tipo di costante

enum constType{i,f,d,s};///int float double string
enum constMeaning{length,angle,count,mat,none};

class GDDconst {
 public:
  GDDconst();
  ~GDDconst();

  void setName(string pname){name=pname;}
  void setValue();

  void setConstType(string pconstType);
  void setconstMeaning(string puType);

  void setNote(string pnote){note=pnote;}
  void setPrimary(bool p){primary=p;}

  string getName()const{return name;}
  void* getValue();
  constType getConstType()const{return constType;}
  constMeaning getConstMeaning()const{return uType;}
  string getNote()const{return note;}
  bool getPrimary()const{return primary;}
 private:
  string name;
  double value;
  constType constType;
  constMeaning uType;
  string note;
  bool primary;/// this attribute  expres if the contant is primary=1 or derived=0
};
#endif //GDDCONST_H








