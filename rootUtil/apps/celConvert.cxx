
#include <rootUtil/CompositeEventList.h>
#include <rootUtil/ComponentsInfoGlast.h>
#include <rootUtil/RuInvertedIndex.h>
#include <rootUtil/RuChain.h>

#if !defined(__CINT__) || defined(__MAKECINT__)

#  include <TFile.h>
#  include <TTree.h>
#  include <TSystem.h>
#  include <TMap.h>
#  include <TObjArray.h>
#  include <TObjString.h>
#  include <TChainIndex.h>

#  include <stdlib.h>
#  include <fstream>
#  include <iostream>
#  include <vector>

#endif

#ifndef WIN32
#include <unistd.h>
#else
#include "apps/winutil/XGetopt.h"
#endif

std::string theApp("celConvert.exe") ;

  
//====================================================================
// 
//====================================================================



//====================================================================
// 
//====================================================================

void print_cel_header( const char* fileName )
 {
  std::cout
    <<"\n"
    <<"#\n"
    <<"#! CEL TXT 0.1\n"
    <<"#! Converted from ROOT CEL file : "<<fileName<<'\n'
    <<"#\n"
    <<std::flush ;
 }

void print_files
 ( const char* fileName,
   unsigned long int beginEvent,
   unsigned long int numEvents )
 {
  std::cout<<"\n#! SECTION Files\n"<<std::flush ;
  
  CompositeEventList cel(fileName) ;
  if ( cel.isOk()==kFALSE )
   {
    std::cerr<<"Can't open file "<<fileName<<std::endl ;
    return ;
   }
  if (numEvents==0)
   { numEvents = cel.numEvents()-beginEvent ; }
  Long64_t endEvent = beginEvent+numEvents ;
  
  // discover components and prepare arrays
  cel.shallowRead(beginEvent) ; // enforce discovery of components
  UInt_t iComponent, numComponents = cel.numComponents() ;
  TObjArray singleFileNames(numComponents) ;
  TObjArray orderedFileNames(numComponents) ;
  TObjArray treeNames(numComponents) ;
  for ( iComponent = 0 ; iComponent < numComponents ; ++iComponent )
   {
    singleFileNames[iComponent] = new TMap ;        // used as a set
    orderedFileNames[iComponent] = new TObjArray ; // preserve order
    treeNames[iComponent] = new TObjString(*(cel.treeName(iComponent))) ;
   }
  
  // loop on events
  const char* cFileName ;
  TObjString * objFileName ;
  Long64_t iEvent ;
  for ( iEvent = beginEvent ; iEvent < endEvent ; ++iEvent )
   {
    cel.shallowRead(iEvent) ;
    for ( iComponent = 0 ; iComponent < numComponents ; ++iComponent )
     {
      cFileName = cel.fileName(iComponent)->GetName() ;
      if ((((TMap*)singleFileNames[iComponent])->FindObject(cFileName))==0)
       {
        objFileName = new TObjString(cFileName) ;
        ((TMap*)singleFileNames[iComponent])->Add(objFileName,new TObject) ;
        ((TObjArray*)orderedFileNames[iComponent])->Add(objFileName) ;
        if (cel.treeName(iComponent)->IsEqual(treeNames[iComponent])!=kTRUE)
         {
          std::cerr<<"Mutable tree names not yet supported"<<std::endl ;
          return ;
         }
       }
     }
   }  
  
  // print file paths
  TIter * fileNameIter ;
  for ( iComponent = 0 ; iComponent < numComponents ; ++iComponent )
   {
    fileNameIter = new TIter((TObjArray*)orderedFileNames[iComponent]) ;
    while ((objFileName=((TObjString *)fileNameIter->Next())))
     {
      std::cout
        <<'('<<cel.componentName(iComponent)<<')'
        <<objFileName->GetName()
        <<std::endl ;
     }
    delete fileNameIter ;
   }
 }

void print_events
 ( const char* fileName,
   unsigned long int beginEvent,
   unsigned long int numEvents )
 {
  std::cout<<"\n#! SECTION Events\n"<<std::flush ;
  CompositeEventList cel(fileName) ;
  if ( cel.isOk()==kFALSE )
   {
    std::cerr<<"Can't open file "<<fileName<<std::endl ;
    return ;
   }
  if (numEvents==0)
   { numEvents = cel.numEvents()-beginEvent ; }
  Long64_t endEvent = beginEvent+numEvents ;
  
  // discover components and prepare info
  ComponentsInfoGlast componentsInfo ;
  cel.shallowRead(beginEvent) ; // enforce discovery of components
  UInt_t iComponent, numComponents = cel.numComponents() ;
  std::vector<const ComponentInfo *> componentsInfoArray(numComponents) ;
  std::vector<RuChain *> componentsChains(numComponents) ;
  std::vector<RuInvertedIndex *> componentsIndexes(numComponents) ;
  const ComponentInfo * info ;
  RuChain * chain ;
  RuInvertedIndex * index ;
  for ( iComponent = 0 ; iComponent < numComponents ; ++iComponent )
   {
    info = componentsInfo.getInfo(cel.componentName(iComponent)) ;
    componentsInfoArray[iComponent] = info ;
    if ((info->runIdBranchName=="")||(info->eventIdBranchName==""))
     { continue ; }
    chain = cel.newChain(iComponent) ;
    componentsChains[iComponent] = chain ;
    index = new RuInvertedIndex(chain,info->runIdBranchName,info->eventIdBranchName) ;
    componentsIndexes[iComponent] = index ;
   }

  // loop on events
  Int_t runID1, runID2 ;
  Int_t eventID1, eventID2 ;
  Bool_t found1, found2 ;
  Long64_t iEvent, iEntry ;
  for ( iEvent = beginEvent ; iEvent < endEvent ; ++iEvent )
   {
    found1 = kFALSE ;
    cel.shallowRead(iEvent) ;
    for ( iComponent = 0 ; iComponent < numComponents ; ++iComponent )
     {
      info = (const ComponentInfo *)componentsInfoArray[iComponent] ;
      if ((info->runIdBranchName=="")||(info->eventIdBranchName==""))
       { continue ; }
      iEntry = cel.entryIndex(iComponent) ;
      chain = componentsChains[iComponent] ;
      index = componentsIndexes[iComponent] ;
      if (found1==kFALSE)
       {
        found1 = index->GetMajorMinorWithEntry(iEntry,runID1,eventID1) ;
        if (found1==kFALSE)
         {
          std::cerr
            <<"Can't find event "<<iEvent
            <<" which is entry "<<iEntry
            <<" in reverse index of component "<<cel.componentName(iComponent)
            <<std::endl ;
          return ;
         }
       }
      else
       {
        found2 = index->GetMajorMinorWithEntry(iEntry,runID2,eventID2) ;
        if (found2==kFALSE)
         {
          std::cerr
            <<"Can't find event "<<iEvent
            <<"which is entry "<<iEntry
            <<" in reverse index of component "<<cel.componentName(iComponent)
            <<std::endl ;
          return ;
         }
        if ((runID2!=runID1)||(eventID2!=eventID1))
         {
          std::cerr
            <<"Inconsistent run/event IDs for event "<<iEvent
            <<std::endl ;
          return ;
         }
       }
     }
    std::cout<<runID1<<" "<<eventID1<<std::endl ;
   }  
 }

void usage()
 {

  std::cout << std::endl
       << theApp << std::endl
       << "Display on the ouput channel the textual flavor of a ROOT Composite Event List" << std::endl
       << std::endl ;
  
  std::cout << "Usage:" << std::endl
       << "\t" << theApp << " [options] RootCelFile" << std::endl 
       << std::endl
       << std::endl
       << "\tOPTIONS for all jobs" << std::endl
       << "\t   -h                : print this message" << std::endl
       << "\t   -s <flags>        : what sections to convert ([fe]+)"      << std::endl
       << "\t   -n <nEvents>      : run over <nEvents>" << std::endl
       << "\t   -b <beginEvent>   : begin with event <beginEvent>" << std::endl
       << std::endl ;
 }
  

int main(int argn, char** argc)
 {
# ifdef WIN32
    gSystem->Load("libTreePlayer.dll");
# endif

  std::string sections("fe") ;
  unsigned long int beginEvent(0) ;
  unsigned long int nEvent(0) ;

  // parse options
  char * endPtr ;  
  int opt ;
  while ( (opt = getopt(argn, argc, "hp:n:s:")) != EOF )
   {
    switch (opt)
     {
      case 'h':   // help      
        usage();
        return 1;
      case 's':   //  output
        sections = std::string(optarg) ;
        break ;
      case 'n':   // number of events
        nEvent = strtoul(optarg,&endPtr,0) ;
        break ;
      case 'b':   // begin event
        beginEvent = strtoul(optarg,&endPtr,0) ;
        break;
      case '?' :
        usage() ;
        return 2 ;
      default:
        std::cerr<<opt<<"not parsable..."<<std::endl ;
        return 2 ;      
     }
   }

  if ( argn - optind == 0 )
   {
    std::cerr << theApp << " require an input file" << std::endl
         << "Try " << theApp << " -h" << std::endl ;
    return 3 ;
   }

  if ( argn - optind > 1 )
   {
    std::cerr << theApp << " accept only one input file" << std::endl
         << "Try " << theApp << " -h" << std::endl ;
    return 4 ;
   }

  const char * fileName = argc[optind] ;
  print_cel_header(fileName) ;
  std::size_t iChar, nChar = sections.size() ;
  for ( iChar = 0 ; iChar<nChar ; ++iChar )
   {
    switch (sections[iChar])
     {
      case 'f' :
        print_files(fileName,beginEvent,nEvent) ;
        break ;
      case 'e' :
        print_events(fileName,beginEvent,nEvent) ;
        break ;
      default:
        std::cerr<<opt<<"not parsable..."<<std::endl ;
        return 5 ;      
     }
   }
  std::cout<<std::endl ;
    
  return 0 ;
 }
