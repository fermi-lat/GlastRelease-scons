#ifndef COMPOSITION_H
#define COMPOSITION_H
#include <vector>
#include "detModel/Sections/Ensemble.h"

namespace detModel{

  class Shape;

  /**
   * @author R.Giannitrapani & D.Favretto
   */
  class Composition : public Ensemble {

  public:
    Composition(std::string pName):Ensemble(pName),envelope(0){;}
    Composition():Ensemble(),envelope(0){;}
    ~Composition();

    /**
     * This method set the envelopeRef of the compisition
     * The envolopeRef is used by resolvereference method to find the 
     */
    void setEnvelopeRef(std::string pname){envelopeRef=pname;}
    /**
     * This method set the envelope of the composition
     */
    void setEnvelope(Shape* pvolume){envelope = pvolume;}
    /**
     *
     */
    std::string getEnvelopeRef(){return envelopeRef;}
    /**
     *
     */
    Shape* getEnvelope(){return envelope;}
    /**
     *
     */
    virtual void buildBB();

    /// Composition use this Template Method to put a check on position to be added
    bool checkPosition(Position* p);
  
  private:

    /** without loose of generality we can presume that the envelove is always a simple shape.
     * requirement:
     *  xlm file has to satisfy the proposal document for the geometry description
     * @link aggregation
     * @clientCardinality 1
     * @supplierCardinality 1
     * @label envelope*/
    Shape* envelope;
    std::string envelopeRef;
  };
}
#endif //COMPOSITION_H
