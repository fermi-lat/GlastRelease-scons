//
//========================================================
// File: RuChain.h
// Main developer: D.Chamont.
//========================================================
//


#ifndef rootUtil_RuChain_h
#define rootUtil_RuChain_h

#if !defined(__CINT__) || defined(__MAKECINT__)
#  include <TChain.h>
#  include <vector>
#endif



// Extension of TChain, trying to make it observable
// for few characteristics.

class RuChain : public TChain
 {
  public :
   
    // cannot inherit from TObject, or we will encounter problems
    // of repeated inheritance when trying to derive exiting ROOT
    // classes so to make them observers
    class Observer
     {
      public :
        Observer( RuChain * chain ) : chain_(chain), zombie_(kFALSE)
         { chain_->AddObserver(this) ; }
        virtual ~Observer() {}
        virtual void newTreeLoaded() =0 ;
        void makeZombie() { zombie_ = kTRUE ; }
        Bool_t isZombie() { return zombie_ ; }
      protected :
        RuChain * chain_ ;
      private :
        Bool_t zombie_ ;
     } ;
  
    // Construction/destruction
    RuChain( const char* name, const char* title = "")
     : TChain(name,title), lastTreeNumber_(-1)
     {}
    virtual ~RuChain()
     {
      std::vector<Observer *>::iterator observer ;
      for ( observer = observers_.begin() ; observer != observers_.end() ; ++observer )
       { (*observer)->makeZombie() ; }
     }
    
    // Observers
    void AddObserver( Observer * observer )
     { observers_.push_back(observer) ; }
    
    // Observed TChain fonctions
    virtual Long64_t LoadTree( Long64_t entry )
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
    
  private :
   
    Int_t lastTreeNumber_ ;
    std::vector<Observer *> observers_ ;
    
 } ;


 
#endif

 
 
