#ifndef ENSEMBLE_H
#define ENSEMBLE_H
#include <vector>
#include "detModel/Sections/Volume.h"
#include "detModel/Utilities/Purge.h"
#include "detModel/Sections/Position.h"
#include "detModel/Management/SectionsVisitor.h"

namespace detModel{

  /**
   * @author R.Giannitrapani
   */
  class Ensemble : public Volume {

  public:
    Ensemble(std::string pName):Volume(pName){;};
    Ensemble():Volume(){;};

    virtual ~Ensemble(){purge(positions);};

    /// This method is the recursive Accept for the visitor pattern
    virtual void Accept(SectionsVisitor*);
    /// This method is the non recursive Accept for the visitor pattern
    virtual void AcceptNotRec(SectionsVisitor* v){v->visitEnsemble(this);};

    /** 
     *  This is a Template Method used in the addPosition method; it can be
     *  reimplemented by ensemble subclasses in order to provide check mechanisms 
     *  before adding a position.
     */
    bool checkPosition(Position*){return 1;};
  
    /// This method add a new position to the positions vector 
    void addPosition(Position* ppos){if(checkPosition(ppos)) positions.push_back(ppos);}
    /// This methods return the total number of positions
    int getNumPos(){return positions.size();}
    /// This method returns the vector of positions of the ensamble
    std::vector<Position*> getPositions(){return positions;};

  private:
    /** @link aggregation
     *  @supplierCardinality 1..n */
    std::vector < Position * > positions;
  };
}
#endif //ENSEMBLE_H


