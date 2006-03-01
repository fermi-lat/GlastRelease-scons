#ifndef SECTIONSVISITOR_H
#define SECTIONSVISITOR_H
#include "detModel/Management/Visitor.h"

namespace detModel{

  class Gdd;
  class Section;
  class Box;
  class Tube;
  class Sphere;
  class Ensemble;
  class PosXYZ;
  class AxisMPos;
  class IdField;
  class Seg;

  /**
   * This is the abstract class for sections visitors, i.e. visitors
   * that have to traverse the sections of the XML file */
  class SectionsVisitor : public Visitor {

  public:

    /**
     * This is the visitor for the SectionsContainer
     */
    virtual void visitGdd(Gdd*) = 0;
    /**
     * This is the visitor for the Section 
     */
    virtual void visitSection(Section*) = 0;
    /**
     * This is the visitor for a generic Ensemble
     */
    virtual void visitEnsemble(Ensemble*) = 0;
    /**
     * This is the visitor for the Box 
     */
    virtual void visitBox(Box*) = 0;
    /**
     * This is the visitor for the Tube 
     */
    virtual void visitTube(Tube*) = 0;
    /**
     * This is the visitor for the Tube 
     */
    virtual void visitSphere(Sphere*) = 0;
    /**
     * This is the visitor for the PosXYZ 
     */
    virtual void visitPosXYZ(PosXYZ*) = 0;
    /**
     * This is the visitor for the AxisMPos 
     */
    virtual void visitAxisMPos(AxisMPos*) = 0;
    /**
     * This is the visitor for the IdField 
     */
    virtual void visitIdField(IdField *) = 0;
    /**
     * This is the visitor for the Seg 
     */
    virtual void visitSeg(Seg *) = 0;

    // Do we really need to keep this class pure virtual?  Or
    // would it be better to go ahead and implement setMode 
    // and getMode here to avoid replication in all section visitors?
    // Since Visitor is already not pure virtual, might as
    // well go ahead and implemen...

    /// Allow client to set the mode for Choice elements
    virtual void setMode(std::string pmode) {m_mode = pmode;}
    virtual std::string getMode() {return m_mode;}

    // Need virtual destructor since we have some virtual functions
    // which are not pure virutal

    virtual ~SectionsVisitor() {};
    
  protected:       
    SectionsVisitor() : m_mode(""){};
    std::string m_mode;
  };

}
#endif //SECTIONSVISITOR_H


