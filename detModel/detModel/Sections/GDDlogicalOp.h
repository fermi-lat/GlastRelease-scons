#ifndef GDDLOGICALOP_H
#define GDDLOGICALOP_H
#include "GDDensemble.h"

/**
 * @author D.Favretto & R.Giannitrapani 
 */
class GDDlogicalOp : GDDenvelope {
 public:

  GDDlogicalOp(std::string pname);
  
 private:    
  /** @link aggregation
   *  @supplierCardinality 1..n */
  vector < GDDanySinglePosition * > anySinglePosition;
};

#endif //GDDLOGICALOP_H



