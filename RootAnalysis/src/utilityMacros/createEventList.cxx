/**  @file createEventList.cxx
@fn TEventList* createEventList(TChain *chain, char *cutStr) 
@brief Create a ROOT TEventList.
  @param chain TChain pointer
  @param cutStr string containing the cuts to apply
  
  OUTPUT: TEventList pointer

  To Run:
  .L createEventList.cxx
  TEventList *elist = createEventList(myChain, "Trig_Bits > 3.0");
    
*/


TEventList* createEventList(TChain *chain, char *cutStr) {
    
    
    chain->Draw(">>elist", cutStr);
    TEventList *elist = (TEventList*)gDirectory->Get("elist");
    return elist; 
}


/*! \fn TEventList* createEventList(char *fileName, char *treePath, char *cutStr) 
\brief Create a ROOT TEventList
  \param fileName name of the ROOT file
  \param treePath path and name of the TTree contained in the ROOT file
  \param cutStr string describing the cut to apply when creating the eventlist
  OUTPUT:  pointer to a new TEventList
  
    To Run:
    .L createEventList.cxx
    TEventList *elist = createEventList("myFile.root", "PDR/t1", "Trig_Bits > 3.0");
*/

TEventList* createEventList(char *fileName, char *treePath, char *cutStr) {
    
    TFile *f = new TFile(fileName);
    TTree *t = (TTree*)f->Get(treePath);
    t->Draw(">>elist", cutStr);
    TEventList *elist = (TEventList*)gDirectory->Get("elist");
    return elist;
}