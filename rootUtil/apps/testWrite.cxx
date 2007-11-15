
#include <rootUtil/TestData.h>
#include "TFile.h"
#include "TTree.h"
#include "Riostream.h"

int main(int argc, char **argv)
 {
  char * baseName = "test" ;
  int n =1 ;
  if (argc > 1)
   { baseName = argv[n++] ; } 

  int sc = 0 ;
  try 
   {
    TRandom random ;
    unsigned int numRuns = 2 ;
    unsigned int numEvents = 3 ;
    Long64_t runID ;
    for ( runID = 0 ; runID<numRuns ; ++runID )
     {
	  sc = testWrite<TestDigiLabel>(baseName,runID,numEvents,random) ;
	  sc = testWrite<TestReconLabel>(baseName,runID,numEvents,random) ;
     }
   }
  catch (...)
   {
    std::cout<<"AN UNKNOWN EXCEPTION HAS BEEN RAISED"<<std::endl ;
    sc = 1 ;
   }
         
  if (sc == 0)
   { std::cout<<"TEST ROOT file writing succeeded!"<<std::endl ; }
  else
   { std::cout<<"FAILED writing"<<std::endl ; }

  return(sc) ;
 }


