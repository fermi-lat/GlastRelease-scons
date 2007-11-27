
#include <rootUtil/TestData.h>
#include <rootUtil/CompositeEventList.h>
#include "Riostream.h"

int main( int argc, char ** argv )
 {
  char * baseName = "test" ;
  int argi = 1 ;
  if (argc > 1)
   { baseName = argv[argi++] ; } 

  int sc = 0 ;
  try 
   {
	  
    TestReader reader ;
    reader.add(baseName,new TestDigiComponent) ;
    reader.add(baseName,new TestReconComponent) ;
    reader.showByComponent() ;
    reader.showByEvent() ;
    
    CompositeEventList cel ;
    TString fileName = baseName ;
    fileName += ".cel.root" ;
    cel.openCelFile(fileName) ;
    cel.printout() ;
    cel.closeCelFile() ;
    
   }
  catch (...)
   {
    std::cout<<"AN UNKNOWN EXCEPTION HAS BEEN RAISED"<<std::endl ;
    sc = 1 ;
   }
         
  if (sc == 0)
   { std::cout << "TEST ROOT files reading succeeded!" << std::endl ; }
  else
   { std::cout << "FAILED recon reading" << std::endl ; }

  return(sc) ;
 }
