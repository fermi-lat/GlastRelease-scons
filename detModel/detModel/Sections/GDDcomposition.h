#ifndef GDDCOMPOSITION_H
#define GDDCOMPOSITION_H
#include <vector>
#include "GDDensemble.h"
#include "GDDensemble.h"

class GDDshape;

/**
 * @author R.Giannitrapani & D.Favretto
 */
class GDDcomposition : public GDDensemble {
public:
  GDDcomposition(std::string pName):GDDensemble(pName),envelope(0){;}
  GDDcomposition():GDDensemble(),envelope(0){;}
  ~GDDcomposition();

  /**
   * This method set the envelopeRef of the compisition
   * The envolopeRef is used by resolvereference method to find the 
   */
  void setEnvelopeRef(std::string pname){envelopeRef=pname;}
  /**
   * This method set the envelope of the composition
   */
  void setEnvelope(GDDshape* pvolume){envelope = pvolume;}
  /**
   *
   */
  std::string getEnvelopeRef(){return envelopeRef;}
  /**
   *
   */
  GDDshape* getEnvelope(){return envelope;}
  /**
   *
   */
  virtual void buildBB();
  
 private:

  /** without loose of generality we can presume that the envelove is always a simple shape.
   * requirement:
   *  xlm file has to satisfy the proposal document for the geometry description
   * @link aggregation
   * @clientCardinality 1
   * @supplierCardinality 1
   * @label envelope*/
  GDDshape* envelope;
  std::string envelopeRef;
};
#endif //GDDCOMPOSITION_H
