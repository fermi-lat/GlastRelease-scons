#ifndef GDDSECTIONSVISITOR_H
#define GDDSECTIONSVISITOR_H
#include "detModel/Management/GDDvisitor.h"

class GDD;
class GDDsection;
class GDDbox;
class GDDensemble;
class GDDposXYZ;
class GDDaxisPos;
class GDDaxisMPos;
class GDDidField;
class GDDseg;

/**
 * This is the abstract class for sections visitors, i.e. visitors
 * that have to traverse the sections of the XML file */
class GDDsectionsVisitor : public GDDvisitor {

 public:

  /**
   * This is the visitor for the GDDsectionsContainer
   */
  virtual void visitGDD(GDD*) = 0;
  /**
   * This is the visitor for the GDDsection 
   */
  virtual void visitSection(GDDsection*) = 0;
  /**
   * This is the visitor for a generic GDDensamble
   */
  virtual void visitEnsemble(GDDensemble*) = 0;
  /**
   * This is the visitor for the GDDbox 
   */
  virtual void visitBox(GDDbox*) = 0;
  /**
   * This is the visitor for the GDDposXYZ 
   */
  virtual void visitPosXYZ(GDDposXYZ*) = 0;
  /**
   * This is the visitor for the GDDaxisPos 
   */
  virtual void visitAxisPos(GDDaxisPos*) = 0;
  /**
   * This is the visitor for the GDDaxisMPos 
   */
  virtual void visitAxisMPos(GDDaxisMPos*) = 0;
  /**
   * This is the visitor for the GDDidField 
   */
  virtual void visitIdField(GDDidField *) = 0;
  /**
   * This is the visitor for the GDDseg 
   */
  virtual void visitSeg(GDDseg *) = 0;
};
#endif //GDDSECTIONSVISITOR_H


