
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
    std::cout
      << "=================================\n"
      << "MAKING CEL\n"
      << std::flush ;
    
    TestReader reader ;
    reader.add(baseName,new TestDigiComponent) ;
    reader.add(baseName,new TestReconComponent) ;
       
    TString fileName = baseName ;
    fileName += ".cel.root" ;
    // components can be declared when constructing the
    // cel, and/or one by one afterwards, and before any
    // call to fill methods.
    // below, so to quickly test both methods, I declare
    // one constructor to the constructor, and one afterwards.
    TObjArray componentNames ;
    componentNames.Add(new TObjString(TestDigiComponent::name())) ;
    //componentNames.Add(new TObjString(TestReconComponent::name())) ;
    CompositeEventList cel(fileName,"RECREATE",&componentNames) ;
    cel.declareComponent(TestReconComponent::name()) ;
    if (cel.isOk()==kFALSE)
     {
      std::cerr
        <<"[TestReader] failed to create file "
        <<fileName
        <<std::endl ;
      throw ;
     }
        
    // loop on events
    reader.resetEvent() ;
    std::cout<<"[TestReader] Looping on events"<<std::endl ;
    TRandom random ;
    Long64_t eventIndex, setIndex ;
    while (reader.nextEvent())
     {
      reader.showEvent() ;
      eventIndex = cel.fillEvent(reader.getForest()) ;
      if (eventIndex==-1)
       {
        std::cerr<<"[TestReader] event not filled"<<std::endl ;
        throw ;
       }
      if (random.Uniform()<0.2)
       {
        setIndex = cel.fillFileAndTreeSet() ;
        if (setIndex==-1)
         {
          std::cerr<<"[TestReader] set not filled"<<std::endl ;
          throw ;
         }
       }
     }
    setIndex = cel.fillFileAndTreeSet() ;
    if (setIndex==-1)
     {
      std::cerr<<"[TestReader] set not filled"<<std::endl ;
      throw ;
     }
    cel.writeAndClose() ;
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
