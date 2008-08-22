//
//========================================================
// File: RuInvertedIndex.cxx
// Main developer: D.Chamont.
//========================================================
//


#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <rootUtil/RuStringSet.h>
#  include <TSystem.h>
#  include <Riostream.h>
#else
#  include "RuStringSet.h"
#endif


//====================================================================
// 
//====================================================================

RuStringSet::RuStringSet()
 : elements_(0)
 {}

RuStringSet::RuStringSet( const TString & globalValue )
 : globalValue_(globalValue), elements_(0)
 {}

RuStringSet::RuStringSet( const RuStringSet & stringSet )
 : globalValue_(stringSet.globalValue_), elements_(0)
 {}

const TString & RuStringSet::GetString() const
 { return globalValue_ ; }

RuStringSet & RuStringSet::operator=( const RuStringSet & stringSet )
 { return ((*this)=stringSet.globalValue_) ; }

RuStringSet & RuStringSet::operator=( const TString & stringSet )
 {
  globalValue_ = stringSet ;
  delete elements_ ;
  elements_ = 0 ;
  return (*this) ;
 }

RuStringSet::~RuStringSet()
 { delete elements_ ; }
    

//====================================================================
// 
//====================================================================

TIter RuStringSet::MakeIterator() const
 {
//  if (globalValue_.Length()==0)
//   {
//    std::cerr<<"[RuStringSet::MakeIterator] does not make sense for empty set"<<std::endl ;
//    exit(1) ;
//   }
  if (elements_==0)
   { elements_ = globalValue_.Tokenize(':') ; }
  return TIter(elements_) ;
 }

Bool_t RuStringSet::Empty()
 { return (globalValue_.Length()==0?kTRUE:kFALSE) ; }

Int_t RuStringSet::NbElements() const
 {
  if (globalValue_.Length()==0)
   { return 0 ; }
  if (elements_==0)
   { elements_ = globalValue_.Tokenize(':') ; }
  return elements_->GetEntries() ;
 }

Bool_t RuStringSet::Contains( const TString & e )
 {
  if (globalValue_.Length()==0)
   { return kFALSE ; }
  if (elements_==0)
   { elements_ = globalValue_.Tokenize(':') ; }
  TIter elementsItr(elements_) ;
  TObjString * element ;
  while ((element=(TObjString *)elementsItr()))
   {
    if (element->GetString().CompareTo(e)==0)
     { return kTRUE ; }
   }
  return kFALSE ;
 }

void RuStringSet::Add( const TString & e )
 {
  if (Contains(e)==kTRUE)
   { return ; }
  globalValue_ = globalValue_ + ":" + e ;
  if (elements_==0)
   { elements_ = globalValue_.Tokenize(':') ; }
  else
   { elements_->Add(new TObjString(e)) ; }
 }

void RuStringSet::ExpandPathNames()
 {
  TString tmp ;
  if (globalValue_.Length()==0)
   { return ; }
  if (elements_==0)
   { elements_ = globalValue_.Tokenize(':') ; }
  TIter elementsItr(elements_) ;
  TObjString * element ;
  while ((element=(TObjString *)elementsItr()))
   {
    tmp = element->GetString() ;
    gSystem->ExpandPathName(tmp) ;
    element->SetString(tmp) ;
   }
 }


//====================================================================
// 
//====================================================================

std::ostream & operator<<( std::ostream & os, const RuStringSet & set )
 { return (os<<set.GetString()) ; }

 
