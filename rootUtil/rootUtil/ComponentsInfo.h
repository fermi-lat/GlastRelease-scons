//
//========================================================
// File: ComponentsInfo.h
// Main developer: D.Chamont.
//========================================================
//

#ifndef rootUtil_ComponentsInfo_h
#define rootUtil_ComponentsInfo_h 1

#include "../rootUtil/RuStringSet.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TString.h>
#  include <TObjArray.h>
#endif

class ComponentInfo : public TObject
 {
  public :
    TString componentName ;  
    TString shortName ;  
    TString libName ;  
    TString treeName ;  
    TString topBranchName ;  
    TString topBranchType ;  
    TString runIdBranchName ;  
    TString eventIdBranchName ;
 } ;

class ComponentsInfo
 {
  public :
  
    ComponentsInfo() ;
    ~ComponentsInfo() ;
    Bool_t needLibrary( const RuStringSet & componentNames ) const ;
    const ComponentInfo * getInfo( const TString & componentName ) const ;
    TIter getInfoIter() const ;
      
  protected :
  
    ComponentInfo * newInfo() ;
    
  private :

    //ComponentsInfo( const ComponentsInfo & ) ;
    //ComponentsInfo & operator=( const ComponentsInfo & ) ;
   
    TObjArray info_ ;
  
 } ;

#endif
