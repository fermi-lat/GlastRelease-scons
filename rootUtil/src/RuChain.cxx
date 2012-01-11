//
//========================================================
// File: RuChain.cxx
// Main developer: D.Chamont.
//========================================================
//

#include <rootUtil/RuChain.h>
#include <TObjString.h>
#include <TList.h>
#include <TLeaf.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <Riostream.h>
#include <cstdlib>


RuChain::Observer::Observer( RuChain * chain )
 : chain_(chain), zombie_(kFALSE)
 { chain_->AddObserver(this) ; }

RuChain::Observer::~Observer()
 {}

void RuChain::Observer::makeZombie()
 { zombie_ = kTRUE ; }

Bool_t RuChain::Observer::isZombie()
 { return zombie_ ; }

RuChain::RuChain( const char* name, const char* title )
 : TChain(name,title),
   unsignedChainRunIdType_(false), unsignedChainEventIdType_(false),
   lastTreeNumber_(-1)
 {}

RuChain::~RuChain()
 {
  std::vector<Observer *>::iterator observer ;
  for ( observer = observers_.begin() ; observer != observers_.end() ; ++observer )
   { (*observer)->makeZombie() ; }
 }

void RuChain::AddObserver( Observer * observer )
 { observers_.push_back(observer) ; }

Long64_t RuChain::LoadTree( Long64_t entry )
 {
  Long64_t res = TChain::LoadTree(entry) ;
  Int_t treeNumber = GetTreeNumber() ;
  if (lastTreeNumber_!=treeNumber)
   {
    lastTreeNumber_=treeNumber ;
    std::vector<Observer *>::iterator observer ;
    for ( observer = observers_.begin() ; observer != observers_.end() ; ++observer )
     { (*observer)->newTreeLoaded() ; }
   }
  return res ;
 }


//====================================================================
// 
//====================================================================

void RuChain::SetInfo( const ComponentInfo & info )
 { info_ = info ; }

const ComponentInfo & RuChain::GetInfo() const
 { return info_ ; }

void RuChain::CheckRunEventInfo()
 {
  if ( info_.componentName == "" )
   {
    std::cerr<<"[RuChain::CheckRunEventInfo] unknown component"<<std::endl ;
    exit(-1) ;
   }
  if ( info_.runIdBranchName == "" )
   {
    std::cerr
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"unknown branch name for run id"<<std::endl ;
    exit(-1) ;
   }
  if ( info_.eventIdBranchName == "" )
   {
    std::cerr
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"unknown branch name for event id"<<std::endl ;
    exit(-1) ;
   }
  
  TLeaf * runIdLeaf = GetLeaf(info_.runIdBranchName) ;
  TLeaf * eventIdLeaf = GetLeaf(info_.eventIdBranchName) ;
  if ( !runIdLeaf )
   {
    std::cerr
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"cannot find branch "<<info_.runIdBranchName<<std::endl ;
    exit(-1) ;
   }
  if ( !eventIdLeaf )
   {
    std::cerr
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"cannot find branch "<<info_.eventIdBranchName<<std::endl ;
    exit(-1) ;
   }
  
  TString runIdBranchType = runIdLeaf->GetTypeName() ;
  TString eventIdBranchType = eventIdLeaf->GetTypeName() ;
  if ( runIdBranchType == "UInt_t" )
   {
    std::cout
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"unsigned int for run id"<<std::endl ;
    unsignedChainRunIdType_ = true ;
   }
  else if ( runIdBranchType != "Int_t" )
   {
    std::cerr
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"unsupported type for run id: "<<runIdBranchType<<std::endl ;
    exit(-1) ;
   }
  else
   { 
    std::cerr<<"[RuChain::CheckRunEventInfo] OK for "<<info_.componentName<<"."<<info_.runIdBranchName<<std::endl ;
   }
  if ( eventIdBranchType == "UInt_t" )
   {
    std::cout
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"unsigned int for event id"<<std::endl ;
    unsignedChainEventIdType_ = true ;
   }
  else if ( eventIdBranchType != "Int_t" && eventIdBranchType != "UInt_t" )
   {
    std::cerr
      <<"[RuChain::CheckRunEventInfo] ("<<info_.componentName<<") "
      <<"unsupported type for event id: "<<eventIdBranchType<<std::endl ;
    exit(-1) ;
   }
  else
   { 
    std::cerr<<"[RuChain::CheckRunEventInfo] OK for "<<info_.componentName<<"."<<info_.eventIdBranchName<<std::endl ;
   }
 }

void RuChain::CheckRunEventIds( Int_t runId, Int_t eventId )
 {
  if ( unsignedChainRunIdType_ && (runId<0) )
   {
    std::cerr<<"[RuChain::CheckRunEventInfo] too big run id: "<<((UInt_t)runId)<<std::endl ;
    exit(-1) ;
   }
  if ( unsignedChainEventIdType_ && (eventId<0) )
   {
    std::cerr<<"[RuChain::CheckRunEventInfo] too big event id: "<<((UInt_t)eventId)<<std::endl ;
    exit(-1) ;
   }
 }
    
//====================================================================
// 
//====================================================================

struct MyObjString : public TObjString
 {
  MyObjString( const char * str ) : TObjString(str) {}
  virtual Int_t Compare( const TObject *  ) const ;
 } ;

Int_t MyObjString::Compare( const TObject * obj ) const
 {
   if (this == obj) return 0 ;
   if (MyObjString::Class() != obj->IsA()) return -1 ;
   
   TString string1 = GetString() ;
   TString string2 = ((MyObjString*)obj)->GetString() ;
   string1.Remove(string1.Index(".root")) ;
   string2.Remove(string2.Index(".root")) ;
   string1.Remove(0,1+string1.Last('_')) ;
   string2.Remove(0,1+string2.Last('_')) ;
   
   Int_t int1 = string1.Atoi() ;
   Int_t int2 = string2.Atoi() ;
   return (int1<int2?-1:+1) ;
 }


// return the number of files
Int_t RuChain::AddFileSet( const char * firstFileName, Long64_t nentries /* = kBigNumber */ )
 {
  Int_t nf = AddFile(firstFileName,nentries) ;
  TString basename(firstFileName) ;
  
  // extract tree specification
  Int_t dotslashpos = basename.Index(".root/");
  TString behind_dot_root ;
  if (dotslashpos>=0)
   {
    // Copy the tree name specification
    behind_dot_root = basename(dotslashpos+6,basename.Length()-dotslashpos+6) ;
    // and remove it from basename
    basename.Remove(dotslashpos+5) ;
   }

  // extract directory
  Int_t slashpos = basename.Last('/') ;
  TString directory;
  if (slashpos>=0)
   {
    directory = basename(0,slashpos); // Copy the directory name
    basename.Remove(0,slashpos+1);      // and remove it from basename
   }
  else
   {
    directory = gSystem->UnixPathName(gSystem->WorkingDirectory()) ;
   }
  
  // prepare regular expression
  dotslashpos = basename.Index(".root") ;
  TString sReg(basename) ;
  sReg.Remove(dotslashpos) ;
  sReg += "_[0-9]+.root" ;
  TRegexp re(sReg,kTRUE) ;
  
  const char * file ;
  void * dir = gSystem->OpenDirectory(gSystem->ExpandPathName(directory.Data())) ;
  if (dir)
   {
    // establish the list
    TList l ;
    while ((file = gSystem->GetDirEntry(dir)))
     {
      if (!strcmp(file,".") || !strcmp(file,"..")) continue ;
      TString s = file ;
      if (s.Index(re) == kNPOS) continue ;
      l.Add(new MyObjString(file)) ;
     }
    gSystem->FreeDirectory(dir) ;
    
    //sort the files in alphanumeric order
    l.Sort() ;
    TIter next(&l) ;
    TObjString *obj ;
    while ((obj = (TObjString*)next()))
     {
       file = obj->GetName() ;
       if (behind_dot_root.Length() != 0)
          nf += AddFile(Form("%s/%s/%s",directory.Data(),file,behind_dot_root.Data()),kBigNumber) ;
       else
          nf += AddFile(Form("%s/%s",directory.Data(),file),kBigNumber);
     }
    l.Delete() ;
   }
   
   return nf ;
 }

 
 
