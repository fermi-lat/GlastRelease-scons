/** @file chainTrees.cxx
@fn TChain* chainTrees(int numFiles, char* list[], char* treePath) 
@brief  This short macro demonstrates how to create a chain of TTrees contained
in multiple files.  

The resulting TChain pointer is returned.  Users should
check to be sure that pointer is not NULL.  The chain can then be passed to 
another analysis routine, such as RootTreeAnalysis.
  @param numFiles number of files to chain
  @param list list of file names
  @param treePath name of tree or if in a subdirectory, the path to the tree i.e. PDR/t1

  OUTPUT:  Pointer to a TChain

  To Run:
  .L chainTrees.cxx
  char *list[] = { "myFile1.root", "myFile2.root"};
  int numFiles = sizeof(list)/sizeof(char*);
  TChain *myNewChain = chainTrees(numFiles, list, "t1");
  */

TChain* chainTrees(int numFiles, char* list[], char* treePath) {
    
    if (numFiles <= 0) {
        printf("number of Files is <= 0!\n");
        return 0;
    }
    
    TChain *chainedTree = new TChain(treePath);
    
    for(int i=0; i < numFiles; i++){
        
        printf("Opening file %s\n",list[i]);
        
        chainedTree->Add(list[i]);
        
    }
    
    return chainedTree;
    
}