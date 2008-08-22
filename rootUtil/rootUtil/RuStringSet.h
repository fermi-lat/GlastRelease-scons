//
//========================================================
// File: RuStringSet.h
// Main developer: D.Chamont.
//========================================================
//


#ifndef rootUtil_RuStringSet_h
#define rootUtil_RuStringSet_h

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TObjString.h>
#  include <TObjArray.h>
#  include <Riosfwd.h>
#else
#endif

// A set of strings, which are concatenated with
// separator ":" when written or read.

class RuStringSet
 {
  public :
   
    RuStringSet() ;
    explicit RuStringSet( const TString & globalValue ) ;
    RuStringSet( const RuStringSet & ) ;
    RuStringSet & operator=( const RuStringSet & ) ;
    RuStringSet & operator=( const TString & ) ;
    const TString & GetString() const ;
    ~RuStringSet() ;
    
    TIter MakeIterator() const ;
    Bool_t Empty() ;
    Int_t NbElements() const ;
    Bool_t Contains( const TString & element ) ;
    void Add( const TString & element ) ;
    void ExpandPathNames() ;
    
    friend std::ostream & operator<<( std::ostream &, const RuStringSet & ) ;
    
  private :
    
    TString globalValue_ ;
    mutable TObjArray * elements_ ;
 } ;

#endif
