#ifndef GDDFAKEBUILDER_H
#define GDDFAKEBUILDER_H
#include "detModel/Management/GDDbuilder.h"

class GDDsection;
class GDDbox;
class GDDcomposition;
class GDDanyPosition;
class GDDstack;
class GDDanyRelativePosition;
class GDDchoice;
class GDDcontants;
class GDDconst;

/**
 * This is a concrete implementation of GDDbuilder It does not read
 * any file, it simply build a very simple fake geometry hard coded in
 * the implementation. This is useful to test detModel without XML
 * stuff; using this builder there is no need of any external package
 * or library apart from detModel itself. Can be useful for working on
 * a computer on which it is impossible, for one motivation or
 * another, to install Xerces and other GLAST packages. 
 * @author R.Giannitrapani 
 */
class GDDFakeBuilder : public GDDbuilder {
 public:

  /**
   * This is the constructor. 
   */
  GDDFakeBuilder();
  
  /** This method initialize the parser; in this case it does nothing. 
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

};
#endif //GDDFAKEBUILDER_H




