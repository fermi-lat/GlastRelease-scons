#ifndef GDDENSEMBLE_H
#define GDDENSEMBLE_H
#include <vector>
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Sections/GDDposition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author R.Giannitrapani
 */
class GDDensemble : public GDDvolume {

 public:
  GDDensemble(std::string pName):GDDvolume(pName){;};
  GDDensemble():GDDvolume(){;};

  virtual ~GDDensemble(){purge(positions);};

  /// This method is the recursive Accept for the visitor pattern
  virtual void Accept(GDDsectionsVisitor*);
  /// This method is the non recursive Accept for the visitor pattern
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitEnsemble(this);};

  /** 
   *  This is a Template Method used in the addPosition method; it can be
   *  reimplemented by ensemble subclasses in order to provide check mechanisms 
   *  before adding a position.
   */
  bool checkPosition(GDDposition*){return 1;};
  
  /// This method add a new position to the positions vector 
  void addPosition(GDDposition* ppos){if(checkPosition(ppos)) positions.push_back(ppos);}
  /// This methods return the total number of positions
  int getNumPos(){return positions.size();}
  /// This method returns the vector of positions of the ensamble
  std::vector<GDDposition*> getPositions(){return positions;};

 private:
  /** @link aggregation
   *  @supplierCardinality 1..n */
  std::vector < GDDposition * > positions;
};
#endif //GDDENSEMBLE_H


