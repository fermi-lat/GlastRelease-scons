#ifndef GDDXERCESBUILDER_H
#define GDDXERCESBUILDER_H
#include "detModel/Management/GDDbuilder.h"
#include "xml/XmlParser.h"



class GDDsection;
class GDDbox;
class GDDcomposition;
class GDDanyPosition;
class GDDstack;
class GDDanyRelativePosition;
class GDDchoice;

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
  GDDXercesBuilder(char* nameFile);

  /**
   * This method implement the virutal method of GDDbuilder and builds the
   * sections part of the hierarchy, putting it in the private pointer currentSections
   */
  virtual void buildSections();

 private:
  /**
   * This method build the GDDsection object and return a pointer to it
   */
  GDDsection* buildSection(DOM_Node* e);
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
   * This method build a GDDanyPosition object and return a pointer to it
   */
  GDDanyPosition* buildPosition(DOM_Node* e);
  /**
   * This method build a GDDanyrelativePosition object and return a pointer to it
   */
  GDDanyRelativePosition* buildRelativePosition(DOM_Node* e);
  /**
   * This method permits to specify the generic attributes of a position
   */  
  void setAttributePosition(GDDanyPosition* pos, DOM_Node* e);
  /**
   * This method permits to specify the generic attributes of a relative position
   */  
  void setAttributeRelativePosition(GDDanyRelativePosition* pos, DOM_Node* e);

 private:
  /**
   * This private variable is the root element of the DOM hiearchy
   */
  DOM_Document domfile;
  /**  used to verify if a DOM_Node is a comment */
  enum{Comment=8}; 
};
#endif //GDDXERCESBUILDER_H




