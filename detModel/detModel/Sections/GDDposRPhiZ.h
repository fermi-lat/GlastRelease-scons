#ifndef GDDRPHIZ_H
#define GDDRPHIZ_H
#include "detModel/Sections/GDDsinglePos.h"

/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDposRPhiZ : public GDDsinglePos {
public:

  GDDRPhiZ(positiontype ptypePosition);
  /**
   *
   */ 
  void setRot(bool pimpliedRot){impliedRot=pimpliedRot;}
  /**
   *
   */ 
  void setR(double pR){R=PR;}
  /**
   *
   */ 
  void setPhi(double pPhi){pPhi=Phi;}
  /**
   *
   */ 
  void setZ(double pZ){Z=pZ;}
  /**
   *
   */ 
  void setRRef(string pR){RREF=pR;}
  /**
   *
   */ 
  void setPhiRef(string pPhi){PhiREF=pPhi;}
  /**
   *
   */ 
  void setZRef(string pZ){ZREF=pZ;}
  /**
   *
   */ 
  double getR()const{return R;}
  /**
   *
   */ 
  double getPhi()const{return Phi;}
  /**
   *
   */ 
  double getZ()const{return Z;}
  /**
   *
   */ 
  double getRot()const{return impliedRot;}
  /**
   *
   */ 
  string getRRef()const{return RREF;}
  /**
   *
   */ 
  string getPhiRef()const{return PhiREF;}
  /**
   *
   */ 
  string getZRef()const{return ZREF;}
 private:
  bool impliedRot;///default true
  double R; ///default 0
  double Phi; ///default 0
  double Z; ///default 0
  string RREF;
  string PhiREF;
  string ZREF;
};
#endif //GDDRPHIZ_H





