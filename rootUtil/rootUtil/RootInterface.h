#ifndef ROOTINTERFACE_H
#define ROOTINTERFACE_H 1

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"

#include <string>

/** @class RootInterface
    @brief An interface to low-level ROOT file access.

    @author Heather Kelly
    $Header$

  */
class RootInterface {
public:
    RootInterface();
    RootInterface(std::string fileName, std::string openFlag="READ");

    ~RootInterface();

    /// Locate a ROOT object stored in a ROOT file using its name
    TObject* retrieveObjByName(std::string objName);

    /// Locate the first ROOT object stored in a ROOT file according to type
    TObject* retrieveObjByType(std::string objType);

    /// Calls TSystem::Load to dynamically load a library
    int loadDynamicLibrary(std::string libName);

    /// Returns a pointer the opened TFile
    const TFile* getTFile() { return m_rootFile; };

    std::string getFileName() { return m_fileName; };
    int getOpenFlag() { return m_openFlag; };

private:

    /// Private recursive routine used for searching
    TObject* searchKeys(TIter &listofkeys, std::string name, bool &done);

    TFile *m_rootFile;
    std::string m_fileName;

    int m_status;
    int m_openFlag; // 0 = READ, 1 = WRITE, 2 = UPDATE
};

#endif