#ifndef ROOTINTERFACE_H
#define ROOTINTERFACE_H 1

#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"

#include <string>

class RootInterface {
public:
    RootInterface();
    RootInterface(std::string fileName, std::string openFlag="READ");

    ~RootInterface();

    TObject* retrieveObjByName(std::string objName);
    TObject* retrieveObjByType(std::string objType);

    int loadDynamicLibrary(std::string libName);

    const TFile* getTFile() { return m_rootFile; };
    std::string getFileName() { return m_fileName; };
    int getOpenFlag() { return m_openFlag; };

private:
    TObject* searchKeys(TIter &listofkeys, std::string name, bool &done);

    TFile *m_rootFile;
    std::string m_fileName;

    int m_status;
    int m_openFlag; // 0 = READ, 1 = WRITE, 2 = UPDATE
};

#endif