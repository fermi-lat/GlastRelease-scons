#ifndef GDDXERCESBUILDER_H
#define GDDXERCESBUILDER_H
#include "detModel/Management/GDDbuilder.h"
#include "xml/XmlParser.h"



class GDDsection;
class GDDbox;
class GDDcomposition;
class GDDsinglePos;
class GDDstackedPos;
class GDDstack;
class GDDchoice;
class GDDcontants;
class GDDconst;
/**
 * This class is a concrete GDDbuilder that use the Xerces parser and the
 * DOM functionalities to parse the xml file and build the generic model
 * hierarchy
 * @author D.Favretto & R.Giannitrapani
 */
class GDDXercesBuilder : public GDDbuilder {
 public:

  /**
   * This is the constructor. It defines the dom tree (domfile), uses
   * the class Substitute to evaluate all the derived constants and substitute
   * them in the sections. It accepts a string for the name of the file  to open
   */
  GDDXercesBuilder();
  
  /** This method initialize the parser 
   */
  virtual void parseFile(char* nameFile);
  /**
   * This method implement the virutal method of GDDbuilder and builds the
   * sections part of the hierarchy, putting it in the private pointer currentSections
   */
  virtual void buildSections();
  /**
   * This method implement the virtual method of GDDbuilder and buil the constant part of
   * the hierarchy.
   */
  virtual void  buildConstants();

 private:
  /**
   * This method build the GDDsection object and return a pointer to it
   */
  GDDsection* buildSection(DOM_Node* e);
  /**
   * This method build the GDDconst object and return a pointer to it
   */
  GDDconst* buildConst(DOM_Node* e);
  /**
   * This method build a GDDchoice object and return a pointer to it
   */
  GDDchoice* buildChoice(DOM_Node* e);
  /**
   * This method build a GDDcomposition object and return a pointer to it
   */
  GDDcomposition* buildComposition(DOM_Node* e);
  /**
   * This method build a GDDStack object and return a pointer to it
   */
  GDDstack* buildStack(DOM_Node* e);
  /**
   * This method build a GDDBox object and return a pointer to it
   */
  GDDbox* buildBox(DOM_Node* e);
  /**
   * This method build a GDDposition object and return a pointer to it
   */
  GDDsinglePos* buildPosition(DOM_Node* e);
  GDDstackedPos* buildRelativePosition(DOM_Node* e);

  void setAttributePosition(GDDsinglePos* pos, DOM_Node* e);
  void setAttributeRelativePosition(GDDstackedPos* pos, DOM_Node* e);

 private:
  /**
   * This private variable is the root element of the DOM hiearchy
   */
  DOM_Document domfile;
  /**  used to verify if a DOM_Node is a comment */
  enum{Comment=8}; 
};
#endif //GDDXERCESBUILDER_H







