
#ifndef rootUtil_RuUtil_h
#define rootUtil_RuUtil_h

/*
* Project: GLAST
* Package: rootUtil
* File: RuUtil.h
* Authors:
*   DC, David Chamont, LLR, chamont@llr.in2p3.fr
*   EC, Eric Charles , SLAC, echarles@slac.stanford.edu
*/

// ROOT types
#include <Rtypes.h>

// forward declares
class TTree ;
class TFile ;
class TDirectory ;


namespace rootUtil {
  
  class AutoCd
   {
    public :
      explicit AutoCd( TDirectory * newDir =0 ) ;
      ~AutoCd() ;
    private :
      AutoCd( const AutoCd & ) ;
      AutoCd & operator=( const AutoCd & ) ;
      TDirectory * _oldDir ;  
   } ;

  // use OxFFFF as a guard value to persist NULL pointers
  const UShort_t NOKEY(0xFFFF) ; 
  // Open a file.  Return 0 if failed
  TFile * openFile( const char * fileName ) ;
  // Get the TFile associated to a TTree.  Return 0 if failed
  TFile * getFile( TTree & tree ) ;
  
  // Checks to see if the character 'opt' appears in the string 'options'
  Bool_t has_option( const char * options, char opt) ;
  
}

#endif
