{

/** @file chainNtuples.c
Short Macro that demonstrates the use of the chainTree routine
on ntuples
*/


     gROOT->Reset();    // roll back CINT context to last Save
     gROOT->LoadMacro("chainTrees.cxx");  // load your code  

     char * list[] = {
       "backgndmaxpdr100000TUP0000.root",
       "backgndmaxpdr100000TUP0001.root",
       "backgndmaxpdr100000TUP0004.root"};
 
     int numbFiles = sizeof(list)/sizeof(char*);
     TChain* myChain = chainTrees(numbFiles,list, "t1");

}
