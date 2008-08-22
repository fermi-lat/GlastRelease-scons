//
//========================================================
// File: ComponentsInfo.cxx
// Main developer: D.Chamont.
//========================================================
//

// Unhappily enough, I need to use a given include syntax
// for the skimmer interpretation, and another one for
// GlastRelease compilation
#if !defined(__CINT__) || defined(__MAKECINT__)

#  include <rootUtil/ComponentsInfo.h>
#  include <Riostream.h>
#  include <TNamed.h>
#  include <TObjArray.h>
#  include <TObjString.h>
#  include <TString.h>

#else

#  include "ComponentsInfo.h"

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

Bool_t ComponentsInfo::needLibrary( const RuStringSet & componentNames ) const
 {
  TIter componentNamesIter = componentNames.MakeIterator() ;
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
        return kTRUE ;
       }
     }
    delete infoIter ;
    infoIter = 0 ;
   }
  return kFALSE ;
 }
 
 

