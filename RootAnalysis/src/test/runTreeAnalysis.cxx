/**
 *
 * @file runTreeAnalysis
 *
 * @brief Simple test program to exercise RootTreeAnalysis
 *
 * @author no clue
 *
 * $Header$
 *
 */

#include "RootTreeAnalysis.h"
#include "TSystem.h"
#include <iostream>
#include <string>

int main(int argn, char** argc) {

#ifdef WIN32
    gSystem->Load("libTree.dll");
    gSystem->Load("reconRootData.dll");
#endif
    unsigned int numEvents = 100;
    if ( argn > 5) numEvents = atoi(argc[5]);
    const char* path = ::getenv("ROOTANALYSISROOT");
    std::string digiFileName(path);
    digiFileName += "/src/test/digi.root";
    std::string mcFileName(path);
    mcFileName += "/src/test/mc.root";
    std::string reconFileName(path);
    reconFileName += "/src/test/recon.root";
    if (argn > 1) mcFileName = argc[1];
    if (argn > 2) digiFileName = argc[2];
    if (argn > 3) reconFileName = argc[3];
    std::cout << "Setup for Processing" << std::endl;
    RootTreeAnalysis r(digiFileName.c_str(), reconFileName.c_str(), mcFileName.c_str());
    if (argn > 4 && strlen(argc[4]) > 0)
        r.SetHistFileName(argc[4]);
    r.Go(numEvents);
    std::cout << "Done Processing " << numEvents << " Events" << std::endl;
    r.WriteHist();

    return 0;
}
