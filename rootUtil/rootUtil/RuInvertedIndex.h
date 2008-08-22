//
//========================================================
// File: RuInvertedIndex.h
// Main developer: D.Chamont.
//========================================================
//


#ifndef rootUtil_RuInvertedIndex_h
#define rootUtil_RuInvertedIndex_h

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TTreeIndex.h>
#  include <TTreeFormula.h>
#  include <vector>
#endif

class RuChain ;

// A new kind of chain index, meant to do the reverse transformation from
// an entry number to a major/minor pair.
// The implementation favors quickness over memory size : when relevant
// indices already exists in the underlying trees, they are inverted and
// reused ; otherwise the values are read throw a TTreeFormula.

class RuInvertedIndex
 {
  public :
   
    RuInvertedIndex( RuChain * ) ;
    ~RuInvertedIndex() ;
    Bool_t GetMajorMinorWithEntry                      // Affect the current
     ( Long64_t entryNumber,                              // entry of the chain.
       Int_t & majorValue, Int_t & minorValue ) ;
    
  private :
    
    // utility class
    class InvertedTreeIndex
     {
      public :
        InvertedTreeIndex() : fNbEntries(0), fMajorMinorLongs(0) {}
        Long64_t fNbEntries ;
        Long64_t * fMajorMinorLongs ;
     } ;

    // forbid copy
    RuInvertedIndex( const RuInvertedIndex & ) ;
    RuInvertedIndex & operator=( const RuInvertedIndex & ) ;
    
    TTreeIndex * GetTreeIndex( Long64_t index ) ;
    
    RuChain * fChain ;
    std::vector<InvertedTreeIndex> fInvertedIndices ;
    TString fMajorName ; // Index major name
    TString fMinorName ; // Index minor name
    TTreeFormula * fMajorFormula ;
    TTreeFormula * fMinorFormula ;
    Bool_t fZombie ;
 } ;

#endif
