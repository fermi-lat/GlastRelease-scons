//
//========================================================
// File: ComponentsInfo.cxx
// Main developer: D.Chamont.
//========================================================
//

#include <rootUtil/ComponentsInfo.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <Riostream.h>
#  include <TNamed.h>
#  include <TObjArray.h>
#  include <TObjString.h>
#  include <TString.h>
#endif

ComponentsInfo::ComponentsInfo()
 {}
 
ComponentsInfo::~ComponentsInfo()
 { info_.Delete() ; }
 
ComponentInfo * ComponentsInfo::newInfo()
 {
  ComponentInfo * newInfo = new ComponentInfo ;
  info_.Add(newInfo) ;
  return newInfo ;
 }

TIter ComponentsInfo::getInfoIter() const
 { return TIter(&info_) ; }

const ComponentInfo * ComponentsInfo::getInfo( const TString & componentName ) const
 {
  ComponentInfo * currentInfo ;
  TIter infoIter = getInfoIter() ;
  while ((currentInfo=(ComponentInfo *)infoIter()))
   {
    if (componentName==currentInfo->componentName)
     { return currentInfo ; }
   }
  return 0 ;
 }

Bool_t ComponentsInfo::needLibrary( const TString & componentNames ) const
 {
  TObjArray * componentNamesArray = componentNames.Tokenize(':') ;
  TIter componentNamesIter(componentNamesArray) ;
  TObjString * componentName ;
  TIter * infoIter ;
  ComponentInfo * currentInfo ;
  while ((componentName=(TObjString *)componentNamesIter()))
   {
    infoIter = new TIter(getInfoIter()) ;
    while ((currentInfo=(ComponentInfo *)infoIter->Next()))
     {
      if ((componentName->String()==currentInfo->componentName)&&
          (currentInfo->libName!=""))
       {
        delete infoIter ;
        delete componentNamesArray ;
        return kTRUE ;
       }
     }
    delete infoIter ;
    infoIter = 0 ;
   }
  delete componentNamesArray ; 
  return kFALSE ;
 }
 
 

