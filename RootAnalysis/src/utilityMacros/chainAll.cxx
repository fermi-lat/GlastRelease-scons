/** @file chainAll.cxx
 @fn TChain* chainAll(char* dir, char* treePath)
 @brief Chains all TTrees in all ROOT files located in a directory
 @param dir string containing the directory path of the ROOT files
 @param treePath string containing the name and path of the TTree in the ROOT files
          NOTE:  it is expected that the treePath is the same for ALL files
  OUTPUT:  a new TChain object - containing all TTrees in files in the specified directory
*/
TChain* chainAll(char* dir, char* treePath) {

    TChain *chainedTree = new TChain(treePath);

    TString searchPath(dir);
    if (searchPath.EndsWith("/")) {
        searchPath.Append("*.root");
    } else {
        searchPath.Append("/*.root");
    }
    const char *str = searchPath;
    // chains all TTrees - if a file does not contain the appropriate treePath 
    // the TTree is not added is not added to the chain.
    chainedTree->Add(str);

  return chainedTree;

}

