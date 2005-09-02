/** @file copyTree.cxx
@fn void copyTree(char *orgFileName, char *newFileName, char *newDir, char *treeName)
@brief Copies a TTree from one file into another file and storing into a specified directory
  @param orgFileName file containing the TTree we want to copy
  @param newFileName new file to be created to contain copy of the TTree
  @param newDir string containing the directory to store the TTree
  @param treeName string containing the path and name of the TTree to copy
  OUTPUT:  a new file, called newFileName containing a copy of the TTree stored in orgFileName

*/
void copyTree(char *orgFileName, char *newFileName, char *newDir, char *treeName)
{
  TFile fOrg(orgFileName);
  TTree *tOrg = (TTree*)fOrg.Get(treeName);
  TFile fNew(newFileName, "RECREATE");
  fNew.mkdir(newDir);
  fNew.cd(newDir);
  TTree *tNew = tOrg->CopyTree("");
  fNew.Write(0, TObject::kWriteDelete);
  fOrg.Close();
  fNew.Close();
}
