/**
 *
 * @file runTreeMaker
 *
 * @brief Simple program to use the TreeMaker macro in compiled code
 *
 * @authors Michael Kuss, Nicola Omodei, Gloria Spandre
 *
 * $Header$
 *
 */

#include "TreeMaker.h"
#include "TSystem.h"
#include <iostream>
#include <string>

int main(int argn, char** argc) {
  
#ifdef WIN32
    gSystem->Load("libTree.dll");
    gSystem->Load("reconRootData.dll");
#endif

    std::string digiFileName;
    std::string reconFileName;
    std::string mcFileName;
    std::string treeFileName;
    unsigned int numEvents = 10000000;

    if ( argn > 1 ) digiFileName  = argc[1];
    if ( argn > 2 ) reconFileName = argc[2];
    if ( argn > 3 ) mcFileName    = argc[3];
    if ( argn > 4 ) treeFileName  = argc[4];
    if ( argn > 5 ) numEvents = atoi(argc[5]);
    
    if ( gSystem->AccessPathName(digiFileName.c_str(), kFileExists) ) {
        std::cout << " ===>>> Wrong Digifile name!!!" << std::endl;
        return 1;
    }
    
    if ( reconFileName != "" && gSystem->AccessPathName(reconFileName.c_str(), kFileExists) ) {
        std::cout << " ===>>> Reconfile not existing!!!" << std::endl;
        return 2;
    }

    if ( mcFileName != "" && gSystem->AccessPathName(mcFileName.c_str(), kFileExists) ) {
        std::cout << " ===>>> MCfile not existing!!!" << std::endl;
        return 3;
    }
    if ( treeFileName == "" )
        treeFileName = "MyRootFile.root";

    TreeMaker r(digiFileName.c_str(), reconFileName.c_str(), mcFileName.c_str(), treeFileName.c_str());
    
    std::cout << "Setup for Processing" << std::endl;
    std::cout << "\"" << digiFileName.c_str()
              << "\" \"" << reconFileName.c_str()
              << "\" \"" << mcFileName.c_str()
              << "\" \"" << treeFileName.c_str() << "\""
              << std::endl;
    r.CreateTree(numEvents);
    std::cout << "Done Processing " << numEvents << " Events" << std::endl;
    return 0;
}
