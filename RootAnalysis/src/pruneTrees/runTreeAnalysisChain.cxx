#include "./src/pruneTrees/PruneTree.h"
#include "TSystem.h"
#include <iostream>
#include <string>
#include <stdio.h>

/** @file runTreeAnalysisChain
@brief Run modified RootTreeAnalysis that prunes sets of mc, digi, 
recon files,

Input params are:
text file containing list of mc files to chain for pruning
text file containing list of recon files to chain for pruning
text file containing list of digi files to chain for pruning
*/
int main(int argn, char** argc) {
    
  std::cout << "Chain Pruner with " << argn << " args " << std::endl; 

#ifdef WIN32
    gSystem->Load("libTree.dll");
    gSystem->Load("reconRootData.dll");
#endif
    unsigned int numEvents = 0;
    const char* path = ::getenv("ROOTANALYSISROOT");
    std::string digiFileName(path);
    digiFileName += "/src/test/digi*";
    std::string mcFileName(path);
    mcFileName += "/src/test/mc*";
    std::string reconFileName(path);
    reconFileName += "/src/test/recon*";
    std::string anaTupFileName(path);
    anaTupFileName += "/src/test/ntuple*";
    TChain* mc=0;
    TChain* digi=0;
    TChain* recon=0;
    TChain* anaTup=0;

    if (argn > 1) {
        //mcFileName = std::string(argc[1]) + "*";
	//	std::cout << "MC fn " << mcFileName  << std::endl; 
        FILE *ifstr;
        ifstr = fopen(argc[1], "r");
        char name[2048]="";
        std::string namestr(name);
        
	mc = new TChain("Mc");
        while(fscanf(ifstr,"%s", name)!=EOF) {
          namestr+=name;
	  mc->Add(namestr.c_str());
        }
    }
    if (argn > 2) {
        //reconFileName = std::string(argc[2]) + "*";
	//	std::cout << "Recon fn " << reconFileName  << std::endl; 
        FILE *ifstr;
        char name[2048]="";
        std::string namestr(name);
        ifstr = fopen(argc[2], "r");
	recon = new TChain("Recon");
        while(fscanf(ifstr,"%s", name)!=EOF) {
          namestr+=name;
	  recon->Add(namestr.c_str());
        }
    }
    if (argn > 3) {
        //digiFileName = std::string(argc[3]) + "*";
	//	std::cout << "digi fn " << digiFileName  << std::endl; 
        FILE *ifstr;
        ifstr = fopen(argc[3], "r");
        char name[2048]="";
        std::string namestr(name);
	digi = new TChain("Digi");
        while(fscanf(ifstr,"%s", name)!=EOF) {
          namestr+=name;
	  digi->Add(namestr.c_str());
        }
    }

    if (argn > 4) {
        anaTupFileName = std::string(argc[4]) + "*";
	//	std::cout << "AnaTup fn " << anaTupFileName  << std::endl; 
	anaTup = new TChain("1");
	anaTup->Add(anaTupFileName.c_str());
    }

    std::cout << "Setup for Processing" << std::endl;
    RootTreeAnalysis r(digi, recon, mc, anaTup);
    r.CopyTrees();
    r.Go();
    std::cout << "Done Processing " << numEvents << " Events" << std::endl;
    r.WriteHist();

    return 0;
}







