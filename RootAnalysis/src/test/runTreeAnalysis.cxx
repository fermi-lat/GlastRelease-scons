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

#include "../RootTreeAnalysis/RootTreeAnalysis.h"
#include "TSystem.h"
#include <iostream>
#include <string>
#include "facilities/commonUtilities.h"
int main(int argn, char** argc) {

#ifdef WIN32
    gSystem->Load("libTree.dll");
    gSystem->Load("libTreePlayer.dll");
    gSystem->Load("reconRootData.dll");
#endif
    Long64_t numEvents = 25;
    std::string dataPath = 
      facilities::commonUtilities::getDataPath("rootTestData");
    if ( argn > 5) numEvents = atoi(argc[5]);
    std::string digiFileName(dataPath);
    digiFileName += "/default/digi.root";
    std::string mcFileName(dataPath);
    mcFileName += "/default/mc.root";
    std::string reconFileName(dataPath);
    reconFileName += "/default/recon.root";
    if (argn > 1) mcFileName = argc[1];
    if (argn > 2) digiFileName = argc[2];
    if (argn > 3) reconFileName = argc[3];
    std::string histFileName = "Histograms.root";
    if (argn > 4 && strlen(argc[4]) > 0)
        histFileName = argc[4];
    std::cout << "Setup for Processing" << std::endl;
    RootTreeAnalysis r(digiFileName.c_str(), reconFileName.c_str(),
                       mcFileName.c_str(),
                       const_cast<char*>(histFileName.c_str()));
    std::cout << "        histogram: " << r.m_histFileName << std::endl;
    r.Go(numEvents);
    std::cout << "Done Processing " << numEvents << " Events" << std::endl;
    r.WriteHist();

    return 0;
}
