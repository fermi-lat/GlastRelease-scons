#ifndef LOGICALOP_H
#define LOGICALOP_H
#include "detModel/Sections/Ensemble.h"

namespace detModel{

  /**
   * @author D.Favretto & R.Giannitrapani 
   */
  class LogicalOp : Envelope {
  public:
    
    LogicalOp(std::string pname);
    
  private:    
    /** @link aggregation
     *  @supplierCardinality 1..n */
    vector < SinglePosition * > SinglePosition;
  };
}
#endif //LOGICALOP_H



