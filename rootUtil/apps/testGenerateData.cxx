
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
    unsigned int numFilesByRun = 3 ;
    unsigned int numEventsByFile = 4 ;
    Long64_t runID, fileIndex, firstEvent ;
    for ( runID = 0 ; runID<numRuns ; ++runID  )
     {
      for ( fileIndex = 0, firstEvent = 0 ;
            fileIndex<numFilesByRun ;
            ++fileIndex, firstEvent += numEventsByFile )
       { sc = testWrite<TestDigiLabel>(baseName,runID,firstEvent,firstEvent+numEventsByFile-1,random) ; }
     }
    for ( runID = 0 ; runID<numRuns ; ++runID )
     {
      for ( fileIndex = 0, firstEvent = 0 ;
            fileIndex<numFilesByRun ;
            ++fileIndex, firstEvent += numEventsByFile )
       { sc = testWrite<TestReconLabel>(baseName,runID,firstEvent,firstEvent+numEventsByFile-1,random) ; }
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


