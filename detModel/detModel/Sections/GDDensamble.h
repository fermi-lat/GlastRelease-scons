#ifndef GDDENSAMBLE_H
#define GDDENSAMBLE_H
#include <vector>
#include "detModel/Sections/GDDvolume.h"
#include "detModel/Utilities/GDDpurge.h"
#include "detModel/Sections/GDDposition.h"
#include "detModel/Management/GDDsectionsVisitor.h"

/**
 * @author R.Giannitrapani
 */
class GDDensamble : public GDDvolume {

 public:
  GDDensamble(std::string pName):GDDvolume(pName){;};
  GDDensamble():GDDvolume(){;};

  virtual ~GDDensamble(){purge(positions);};

  /// This method is the recursive Accept for the visitor pattern
  virtual void Accept(GDDsectionsVisitor*);
  /// This method is the non recursive Accept for the visitor pattern
  virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitEnsamble(this);};
  
  /// This method add a new position to the positions vector 
  void addPosition(GDDposition* ppos){positions.push_back(ppos);}
  /// This methods return the total number of positions
  int getNumPos(){return positions.size();}
  /// This method returns the vector of positions of the ensamble
  std::vector<GDDposition*> getPositions(){return positions;};

 private:
  /** @link aggregation
   *  @supplierCardinality 1..n */
  std::vector < GDDposition * > positions;
};
#endif //GDDENSAMBLE_H


