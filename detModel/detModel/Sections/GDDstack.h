#ifndef GDDSTACK_H
#define GDDSTACK_H
#include <vector>

#include "detModel/Sections/GDDvolume.h"
#include "detModel/Sections/GDDanyRelativePosition.h"

enum stackorigin{atStart,atCentre};

class GDDsectionsVisitor;

/**
 * @author D.Favretto & R.Giannitrapani
 */

class GDDstack :public GDDvolume {
 public:
  GDDstack(stacktype pStype):
    GDDvolume(stack),Stype(pStype),origin(atStart){}
  virtual ~GDDstack();
  /**
   *
   */ 
  virtual void Accept(GDDsectionsVisitor*);
  /**
   *
   */ 
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitStack(this);};
  /**
   *
   */ 
  void addPosition(GDDanyRelativePosition *pPosition){anyRelativePosition.push_back(pPosition);}
  /**
   *
   */ 
  void setOrigin(stackorigin porigin){origin=porigin;}
  /**
   *
   */ 
  stackorigin getOrigin()const{return origin;}
  /**
   *
   */ 
  stacktype getStackType()const{return Stype;}
  /**
   *
   */ 
  std::vector <GDDanyRelativePosition*> getPositions()
  {return anyRelativePosition;};
  /**
   *
   */
  virtual void constructBB();
 private:
  stacktype Stype;
  stackorigin origin; ///default atStart
  /** @link aggregation
   *  @supplierCardinality 1..n */
  std::vector < GDDanyRelativePosition * > anyRelativePosition;
};
#endif //GDDSTACK_H








