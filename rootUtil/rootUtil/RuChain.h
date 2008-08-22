//
//========================================================
// File: RuChain.h
// Main developer: D.Chamont.
//========================================================
//


#ifndef rootUtil_RuChain_h
#define rootUtil_RuChain_h

#include <rootUtil/ComponentsInfo.h>
#include <TChain.h>
#include <vector>


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
        Observer( RuChain * chain ) ;
        virtual ~Observer() ;
        virtual void newTreeLoaded() =0 ;
        void makeZombie() ;
        Bool_t isZombie() ;
      protected :
        RuChain * chain_ ;
      private :
        Bool_t zombie_ ;
     } ;
  
    // Construction/destruction
    RuChain( const char* name, const char* title = "" ) ;
    virtual ~RuChain() ;
    
    // Add files generated due to maximum size
    Int_t AddFileSet( const char * firstFileName, Long64_t nentries = kBigNumber ) ;
    
    // Observers
    void AddObserver( Observer * observer ) ;
    
    // Observed TChain fonctions
    virtual Long64_t LoadTree( Long64_t entry ) ;
    
    // Associated component info
    void SetInfo( const ComponentInfo & info ) ;
    const ComponentInfo & GetInfo() const ;
    void CheckRunEventInfo() ;
    void CheckRunEventIds( Int_t runId, Int_t eventId ) ;
    
  private :
   
    ComponentInfo info_ ;
    bool unsignedChainRunIdType_ ;
    bool unsignedChainEventIdType_ ;
    
    Int_t lastTreeNumber_ ;     
    std::vector<Observer *> observers_ ;
     
 } ;


 
#endif

 
 
