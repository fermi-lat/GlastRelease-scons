
/** @file pruneTree.cxx
@fn void pruneTree (char *orgFileName, char *treePath, TEventList *elist, char *newFileName) 
@brief Creates a new file, containing just those events contained in the input TEventList
  @param orgFileName original file name
  @param treePath tree name
  @param elist TEventList pointer
  @param newFileName name of the new file to be created.
  OUTPUT: A new file, named newFileName that contains all the events that were contained in the TEventList
  
    TO Run:
    .L pruneTree.cxx
    pruneTree("myOrgFile.root", "PDR/t1", elist, "myNewFile.root");
    
*/

void pruneTree (char *orgFileName, char *treePath, TEventList *elist, char *newFileName) {
    
    TFile *fOrg = new TFile(orgFileName);
    TTree *tOrg = (TTree*)fOrg->Get(treePath);
    tOrg->SetEventList(elist);
    TFile *fNew = new TFile(newFileName, "RECREATE");
    TTree *tNew = tOrg->CopyTree("");
    tNew->Write();
    tNew->Print();
    fNew->Close();
    fOrg->Close();
}

/** @fn void pruneTree (int numFiles, char* list[], char* treePath, TCut cut, char *newFileName) 
@brief This version handles a list of files
NOTE:  CopyTree for TChains was not available until ROOT 3.01.
  @param numFiles Number of files to process
  @param list list of files names
  @param treePath directory and tree name,
  @param cut cut to apply when pruning
  @param newFileName name of the new file to contain the results.
  OUTPUT:  numFiles+1 will be created:
            all files will be pruned separately, creating new files 
            then all of the new files will be chained and merged into a new ROOT file.
  */

void pruneTree (int numFiles, char* list[], char* treePath, TCut cut, char *newFileName) {
    
	TChain *chain = chainTrees(numFiles, list, treePath);
	chain->Draw(">>elist", cut);
    chain->SetEventList(elist);
    TFile *fNew = new TFile(newFileName, "RECREATE");
    TTree *tNew = chain->CopyTree("");
    tNew->Write();
    tNew->Print();
    fNew->Close();
}
