/** @file RootInterface.cxx
    @brief Implementation of a low-level interface to ROOT files.

    @author Heather Kelly

    $Header$

  */

#ifndef ROOTINTERFACE_CXX 
#define ROOTINTERFACE_CXX 1

#include "rootUtil/RootInterface.h"
#include "TKey.h"

#include "TSystem.h"

RootInterface::RootInterface() : m_rootFile(0), m_status(0) {
    
}

RootInterface::RootInterface(std::string fileName, std::string openFlag) 
: m_rootFile(0), m_status(0) {
    
    m_rootFile = new TFile(fileName.c_str(), openFlag.c_str());
    if (!m_rootFile) m_status = -1;
}

RootInterface::~RootInterface() {
    if (m_rootFile) {
        m_rootFile->Close();
        delete m_rootFile;
        m_rootFile = 0;
    }
}

int RootInterface::loadDynamicLibrary(std::string libName) {
    // Purpose and Method:  Calls ROOT's TSystem::Load routine
    return (gSystem->Load(libName.c_str()));
}


TObject* RootInterface::retrieveObjByName(std::string objName) {
    // Purpose and Method:  Searches a TFile for an object having the name, objName
    // Input:  objName - the name of the object being searched for
    // Output: A pointer to the TObject if found, NULL otherwise

    TIter nextkey(m_rootFile->GetListOfKeys());
    bool done = false;

    return searchKeys(nextkey, objName, done);
}


TObject* RootInterface::retrieveObjByType(std::string objType) {
    // Purpose and Method:  Searches a TFile for the first occurance of a particular type.
    // Input:  objType - the name of the type being searched for
    // Output:  A pointer to a TObject or NULL if not such object is found.

    TIter nextkey(m_rootFile->GetListOfKeys());
    TKey *key;

    while ((key=(TKey*)nextkey())) {
        std::string keyType = key->ReadObj()->ClassName();
        if (keyType == objType)
            return key->ReadObj();
    }

    return 0;
}

TObject* RootInterface::searchKeys(TIter &listofkeys, std::string name, bool &done) {
    // Purpose and Method:  A recursive routine used to search for a particular name
    // Inputs:  listofkeys - A list of TKeys; name - name of the object being searched for
    //          done - a boolean flag denoting we're finished
    // Output:  A TObject pointer or NULL if no such object was found.

    TKey *key;
    
    while( (key = (TKey*)listofkeys()) && !done ) {
        std::string keyName = key->GetName();
        if (keyName == name) {
            done = kTRUE;
            return key->ReadObj();
        }
        std::string keyType = key->GetClassName();
        if (keyType == "TDirectory") {
            TDirectory* dir = (TDirectory*)key->ReadObj();
            TIter keyList(dir->GetListOfKeys());
            return searchKeys(keyList, name, done);
        }
        
    }
    
    return 0;
    
}


#endif