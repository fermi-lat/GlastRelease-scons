/**
 *
 * @file runTreeMaker
 *
 * @brief Simple program to use the TreeMaker macro in compiled code
 *
 * @authors Michael Kuss, Nicola Omodei
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
    unsigned int numEvents = 10000000;
    const char* path = ::getenv("ROOTANALYSISROOT");

    std::string digiFileName(path);

    digiFileName += "/src/test/digi.root";

    std::string reconFileName(path);
    reconFileName += "/src/test/recon.root";

    std::string mcFileName(path);
    mcFileName += "/src/test/mc.root";
    
    std::string TreeFileName(path);
    TreeFileName += "/src/test/MyRootFile.root";
    
    if (argn > 1) digiFileName  = argc[1];
    if (argn > 2) reconFileName = argc[2];
    if (argn > 3) mcFileName    = argc[3];
    if (argn > 4) TreeFileName  = argc[4];
    if ( argn > 5) numEvents = atoi(argc[5]);
    
    
    TreeMaker r(digiFileName.c_str(), reconFileName.c_str(),
		mcFileName.c_str(), (char*)TreeFileName.c_str());
    
    std::cout << "Setup for Processing" << std::endl;
    std::cout << "\"" << digiFileName.c_str()
              << "\" \"" << reconFileName.c_str()
              << "\" \"" << mcFileName.c_str()
              << "\" \"" << TreeFileName.c_str() << "\""
              << std::endl;
    r.CreateTree(numEvents);
    std::cout << "Done Processing " << numEvents << " Events" << std::endl;
    return 0;
}
