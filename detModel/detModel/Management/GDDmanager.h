#ifndef GDDMANAGER_H
#define GDDMANAGER_H
#include <map>

#include "detModel/Management/GDDmanager.h"
#include "detModel/Management/GDDvisitor.h"
#include "detModel/GDD.h"

class GDDchoice;
class GDDvolume;
class GDDbuilder;

enum buildType{all,constants,sections};

/**
 * This class implement the manager for the sections part of the XML file. 
 * It is a singleton.
 * @author D.Favretto & R.Giannitrapani
 */
class GDDmanager {

public:
  
  /** The destructor must be invoched by the client before leaving the
     application. It will destroy explicitely the GDD object and this
     will start the deallocation of all the objects hierarchy */
  ~GDDmanager();

  /** This method is used to access the pointer to the singleton */
  static GDDmanager* getPointer();
  
  /** This method is used to specify a builder to be used by the manager. */
  void setBuilder(GDDbuilder * b);

  /** This method start the parse of the xml file by the builder.
   */
  void build(buildType);

  /// This method gives back the actual mode for choice elements
  std::string getMode(){return manMode;};

  /// This method sets the actual mode for choice elements
  void setMode(std::string pmode){manMode = pmode;};

  /// This method starts the visit of a hierarchy with a specified visitor
  void startVisitor(GDDvisitor*);

  void cleanGDD();

  /** This method gives back the pointer to the GDD object. It coule be
     used by expert clients to access information without the need of
     a visitor */
  GDD* getGDD(){return manGDD;};

  void setNameFile(char* pname){nameFile = pname;};
  char* getNameFile(){return nameFile;};

 protected:
  /** The constructor is protected; in such a way it is forbidden to *
      build directely a manager. Here the new GDD object is explicitely
   built */
  GDDmanager(){manGDD = new GDD;};

 private:

  /** 
      This is a private pointer to the GDDsectionsContainer that
      provide the entry point to all the generic model hierarchy of
      sections */
  GDD * manGDD;  
  /**
     This is the pointer to the only manager that can be istanciated */
  static GDDmanager* pointer;
  /** This is a private pointer to the builder used by the manager to
      build the hierarchy from the XML file */
  GDDbuilder * manBuilder;
  /// This string represent the mode for choices elements
  std::string manMode;

  char* nameFile;
};
#endif //GDDMANAGER_H
