#ifndef GDDCOMPOSITION_H
#define GDDCOMPOSITION_H
#include <vector>
#include "detModel/Sections/GDDshape.h"
#include "detModel/Sections/GDDanyPosition.h"
#include "detModel/Sections/GDDvolume.h"


/**
 * @author D.Favretto & R.Giannitrapani
 */
class GDDcomposition : public GDDvolume {
public:
    GDDcomposition():GDDvolume(composition),envelope(0){;}
    ~GDDcomposition();

    /// This method is the recursive Accept for the visitor pattern
    virtual void Accept(GDDsectionsVisitor*);
    /// This method is the non recursive Accept for the visitor pattern
    virtual void AcceptNotRec(GDDsectionsVisitor* v){v->visitComposition(this);};
    /**
     * This method add a pointer to new position object to the anyPosition vector 
     */
    void addPosition(GDDanyPosition* ppos){anyPosition.push_back(ppos);}
    /**
     * This method set the envelopeRef of the compisition
     * The envolopeRef is used by resolvereference method to find the 
     */
    void setEnvelopeRef(string pname){envelopeRef=pname;}
    /**
     * This method set the envelope of the composition
     */
    void setEnvelope(GDDvolume* pvolume){envelope=static_cast<GDDshape*>(pvolume);}
    /**
     *
     */
    string getEnvelopeRef(){return envelopeRef;}
    /**
     *
     */
    GDDvolume* getEnvelope(){return envelope;}
    /**
     *
     */
    int getNumPos(){return anyPosition.size();}
    /**
     *
     */
    virtual double getBBX();
    /**
     *
     */
    virtual double getBBY();
    /**
     *
     */
    virtual double getBBZ();
    /**
     *
     */
    vector<GDDanyPosition*> getPositions(){return anyPosition;};

 private:

    /** @link aggregation
     *  @supplierCardinality 1..n */
    vector < GDDanyPosition * > anyPosition;

    /** without loose of generality we can presume that the envelove is always a simple shape.
     * requirement:
     *  xlm file has to satisfy the proposal document for the geometry description
     * @link aggregation
     * @clientCardinality 1
     * @supplierCardinality 1
     * @label envelope*/
    GDDshape* envelope;
//    GDDvolume* prova;
    string envelopeRef;
};
#endif //GDDCOMPOSITION_H
