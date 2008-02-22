//
//========================================================
// File: RuInvertedIndex.cxx
// Main developer: D.Chamont.
//========================================================
//

#include <rootUtil/RuInvertedIndex.h>
#include <rootUtil/RuChain.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TTreeIndex.h>
#  include <TDirectory.h>
#  include <TMath.h>
#  include <Riostream.h>
#endif


//====================================================================
// 
//====================================================================

class MajorMinorPair
 {
  public :
    MajorMinorPair( Long64_t value )
     : longValue(value), majorValue(value>>31), minorValue(value<<33>>33)
     {}
    friend std::ostream & operator<<( std::ostream & os, const MajorMinorPair & mm )
     { return (os<<mm.majorValue<<'|'<<mm.minorValue<<" ("<<mm.longValue<<')') ; }
    Long64_t longValue ;
    Int_t majorValue ;
    Int_t minorValue ;
  } ;
  
class MyFormula : public TTreeFormula, public RuChain::Observer
 {
  public :
    MyFormula( const char *name, const char *formula, RuChain *chain)
     : TTreeFormula(name,formula,chain), RuChain::Observer(chain)
     { SetQuickLoad(kTRUE) ; }
    virtual void newTreeLoaded()
     { UpdateFormulaLeaves() ; }
 } ;

 
//====================================================================
// 
//====================================================================

RuInvertedIndex::RuInvertedIndex
 ( TChain * chain, const char* majorname, const char* minorname )
 : fChain(chain), fMajorName(majorname), fMinorName(minorname),
   fMajorFormula(0), fMinorFormula(0),
   fZombie(kFALSE)
 {
  TString methodName("[RuInvertedIndex::RuInvertedIndex] ") ;

  // Preconditions
  if (fChain->GetNtrees()==0)
   {
    std::cerr<<methodName<<"empty chain"<<std::endl ;
    fZombie = kTRUE ;
    return ;
   }
  
  // Invert existing indices
  TTreeIndex * index ;
  TTree * tree ;
  InvertedTreeIndex invertedIndex ;
  Int_t i ;
  for ( i = 0 ; i < chain->GetNtrees() ; ++i )
   {
    if (fZombie==kTRUE)
     { break ; }
    
    invertedIndex.fNbEntries = 0 ;
    invertedIndex.fMajorMinorLongs = 0 ;
    
    fChain->LoadTree((chain->GetTreeOffset())[i]) ;
    tree = fChain->GetTree() ;
    index = (TTreeIndex *)(tree->GetTreeIndex()) ;
    
    if (index)
     {
      if ((index->IsZombie())||(index->GetN()==0))
       {
        std::cerr<<methodName<<"index problem"<<std::endl ;
        fZombie = kTRUE ;
        break ;
       }
      invertedIndex.fNbEntries = index->GetN() ;
      invertedIndex.fMajorMinorLongs = new Long64_t [invertedIndex.fNbEntries] ;
      Long64_t * indexEntryNumbers = index->GetIndex() ;
      Long64_t * indexMajorMinorLongs = index->GetIndexValues() ;
      Long64_t iEntry, iOriginalEntry ;
      for ( iEntry = 0 ; iEntry < invertedIndex.fNbEntries ; ++iEntry )
       {
        iOriginalEntry = indexEntryNumbers[iEntry] ;
        if ((iOriginalEntry<0)||(iOriginalEntry>=invertedIndex.fNbEntries))
         {
          std::cerr
            <<methodName<<"unexpected entry number "
            <<iOriginalEntry<<std::endl ;
          fZombie = kTRUE ;
          break ;
         }
        invertedIndex.fMajorMinorLongs[iOriginalEntry]
         = indexMajorMinorLongs[iEntry] ;
       }
     }
    
    fInvertedIndices.push_back(invertedIndex) ;
   }
 }

RuInvertedIndex::~RuInvertedIndex()
 {
  unsigned int  i ;
  for ( i = 0 ; i < fInvertedIndices.size() ; ++i )
   { delete [] fInvertedIndices[i].fMajorMinorLongs ; }
  delete fMajorFormula ;
  delete fMinorFormula ;
 }

Bool_t RuInvertedIndex::GetMajorMinorWithEntry
 ( Long64_t entryNumber, 
   Int_t & majorValue, Int_t & minorValue )
 {
  // Preconditions
  TString methodName("[RuInvertedIndex::GetMajorMinorWithEntry] ") ;
  if (fZombie==kTRUE)
   {
    std::cerr<<methodName<<"this reverse index is a zombie"<<std::endl ;
    return kFALSE ;
   }
  if (fInvertedIndices.size() == 0)
   {
    std::cerr<<methodName<<"no reverse index"<<std::endl ;
    return kFALSE ;
   }
  if ((entryNumber<0)||(entryNumber>=fChain->GetEntries()))
   {
    std::cerr
      <<methodName<<"unexpected entry number "
      <<entryNumber<<std::endl ;
    return kFALSE ;
   }

  // Get the majorname/minorname
  Long64_t treeEntryNumber = fChain->LoadTree(entryNumber) ;
  Int_t treeNo = fChain->GetTreeNumber() ;
  if (fInvertedIndices[treeNo].fNbEntries!=0)
   {
    MajorMinorPair mm(fInvertedIndices[treeNo].fMajorMinorLongs[treeEntryNumber]) ;
    majorValue = mm.majorValue ;
    minorValue = mm.minorValue ;
   }
  else
   {
    if (fMajorFormula==0)
     {
      fMajorFormula = new MyFormula("Major",fMajorName,fChain) ;
      fMinorFormula = new MyFormula("Minor",fMinorName,fChain) ;
     }
    majorValue = (Long64_t)fMajorFormula->EvalInstance() ;
    minorValue = (Long64_t)fMinorFormula->EvalInstance() ;
   }
  
  return kTRUE ;
 }
