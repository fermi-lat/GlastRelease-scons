#if !defined(__CINT__)
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

TChain* chainTrees(int numFiles, char* list[], char* treePath, bool verbose=false) {

    if (numFiles <= 0) {
        printf("number of Files is <= 0!\n");
        return 0;
    }
    else if (verbose){
      printf("\nRequesting %d files to be merged\n",numFiles);
    }

    TChain *chainedTree = new TChain(treePath);

    for(int i=0; i < numFiles; i++){

      if (verbose) printf("Opening file %s\n",list[i]);
      
      chainedTree->Add(list[i]);

    }

    return chainedTree;

}

void copyTree(TFile *f, TChain *chainedTree) {

    f->cd();
    TTree *tNew = chainedTree->CopyTree("");
    f->Write(0, TObject::kWriteDelete);
}

void help()
{
  std::cout<<"usage : executable <OutputFile> <InputFile> "<<std::endl;
  std::cout<<"where InputFile is a text file with the path to the input files, one per line. "<<std::endl;
}

std::vector<std::string> parseInputFiles(char* inputFileName)
{
  std::string s(inputFileName);
  std::vector<std::string> listFiles;

  std::ifstream f1(inputFileName);
  if (!f1.is_open()) 
    {
      std::cout<<"Did not find Input file "<<inputFileName<<std::endl;
      exit(1);
    }
  char buf[200];
  //  f1.getline(buf,100);
  while(f1.getline(buf,100))
    {
      std::cout<<buf<<std::endl;
      listFiles.push_back(std::string(buf));
    } 
  f1.close();
  return listFiles;
}

/// Merge any number of ROOT files with the same structure.  The files can 
/// contain any number of TTrees (Other ROOT objects are not currently
/// handled.  To use:
/// run mergeRootFiles.exe and provide two input parameters:
/// name of new output ROOT files and the name of the text file containing
/// the list of ROOT files to merge.
int main(int argn, char** argc) {

  printf("reading entries...\n");
  char *outFileName("mergedFile.root");
  if (argn < 2) {
    help();
    exit(1);
  }
  outFileName = argc[1];
  printf("creating output file with name: %s\n",outFileName);
  char *inputFileName = argc[2];

  std::vector<std::string> listFiles = parseInputFiles(inputFileName);
  
  
  int numFiles=listFiles.size();
  //char **fileList;//[numFiles];
  char **fileList = (char**)calloc(numFiles, sizeof(char *));
  int i=0;
  for(i=0;i<numFiles;i++)
    {
      fileList[i] =  const_cast<char*> (listFiles[i].c_str());
    }

  ///  can handle up to 50 TTrees per file
  char *listTrees[50];
  int numTrees = 0;
  
  TFile *fOut = new TFile(outFileName, "RECREATE");
  
  // Use one input file to figure out the names of all the TTrees that we 
  // want to chain
  //    TFile f(oneInputFileName, "READ");
  printf("Reading %s to retrieve content information...\n",listFiles[0].c_str());
  TFile f(listFiles.front().c_str(), "READ");
  TList *keyList = f.GetListOfKeys( );
  TIter keyIter(keyList);
  while (TKey *key = (TKey*)keyIter.Next()) {
    const char* keyClassName=key->GetClassName( );
    char* atree="TTree";
    //       if (strcmp(keyClassName,atree)==0)  {
    printf(" Found %s %d : %s;",keyClassName,numTrees,key->GetName());
    listTrees[numTrees] = const_cast<char*>(key->GetName());
    numTrees++;
    
  }
  
  // now loop oer all the trees, make a chain..and then copy it into the 
  // output file.
  printf("Performing merge - this may take a few minutes\n");
  for (int i = 0; i < numTrees; i++) {
    TChain *chain = chainTrees(numFiles, fileList, listTrees[i]);
    copyTree(fOut, chain);
  }

  printf("Merge Complete\n");
  
  return 0;
}

