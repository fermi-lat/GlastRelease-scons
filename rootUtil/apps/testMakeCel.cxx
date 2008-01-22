
#include <rootUtil/TestData.h>
#include <rootUtil/CompositeEventList.h>
#include <TObjString.h>
#include <Riostream.h>

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
       
    TString fileName = baseName ;
    fileName += ".cel.root" ;
    TObjArray componentNames ;
    componentNames.Add(new TObjString(TestDigiComponent::name())) ;
    componentNames.Add(new TObjString(TestReconComponent::name())) ;
    
    CompositeEventList cel(fileName,"RECREATE",&componentNames) ;
    if (cel.isOk()==kFALSE)
     {
      std::cerr
        <<"[TestReader] failed to create file "
        <<fileName
        <<std::endl ;
      sc = 1 ;
     }
    else
     {    
      reader.resetEvent() ;
      std::cout<<"[TestReader] Looping on events"<<std::endl ;
      TRandom random ;
      while (reader.nextEvent())
       {
        reader.showEvent() ;
        cel.fillEvent(reader.getForest()) ;
        if (random.Uniform()<0.2)
         { cel.fillFileAndTreeSet() ; }
       }
      cel.fillFileAndTreeSet() ;
      cel.writeAndClose() ;
     }
   }
  catch (...)
   {
    std::cout<<"AN UNKNOWN EXCEPTION HAS BEEN RAISED"<<std::endl ;
    sc = 1 ;
   }
         
  if (sc == 0)
   { std::cout << "TEST ROOT files listing succeeded!" << std::endl ; }
  else
   { std::cout << "FAILED recon reading" << std::endl ; }

  return(sc) ;
 }
