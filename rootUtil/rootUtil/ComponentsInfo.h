//
//========================================================
// File: ComponentsInfo.h
// Main developer: D.Chamont.
//========================================================
//


#ifndef rootUtil_ComponentsInfo_h
#define rootUtil_ComponentsInfo_h 1

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TObject.h>
#  include <TString.h>
#  include <TObjArray.h>
#  include <TMap.h>
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
    Bool_t needLibrary( const TString & componentNames ) const ;
    const ComponentInfo * getInfo( const TString & componentName ) const ;
    TIter getInfoIter() const ;
      
  protected :
  
    ComponentInfo * newInfo() ;
    
  private :

    ComponentsInfo( const ComponentsInfo & ) ;
    ComponentsInfo & operator=( const ComponentsInfo & ) ;
   
    TObjArray info_ ;
  
 } ;

#endif
