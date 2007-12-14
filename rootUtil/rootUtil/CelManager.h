
#ifndef CelManager_h
#define CelManager_h

/** 
* @file CelManager.h
* @brief The class CelManager is used to set up and handle the composite event list (cel) root files
*
* Originally made to connect the CEL world to the needs of RootIoSvc.
* Meant to disappear, once we have understood what is lacking in CompositeEventList,
* and what shoudl return to RootIoSvc.
* 
* File: $Header$
* Authors:
*   HK, Heather Kelly, heather@lheapop.gsfc.nasa.gov
*   DC, David Chamont, LLR, chamont@poly.in2p3.fr
*/

#include "CompositeEventList.h"
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <vector>
#include <map>


class CelManager
 {
  public :
  
    CelManager( Bool_t verbose =false ) ;
    ~CelManager() ; 

    /// Writing
    Bool_t initWrite
     ( const TString & celFileName ="cel.root",
       const TString & options ="RECREATE" ) ;
    Bool_t addComponent( const TString & compName, TTree * t ) ;
    Bool_t fillEvent() ; 
    Bool_t fillFileAndTreeSet() ;

    /// Reading
    Bool_t initRead( const TString & celFileName ="cel.root" ) ;
    Long64_t getNumEvents() ;
    Long64_t getEventIndexInTree( const TString & treeName, Long64_t index ) ;
    TChain * getChainByType( const TString & treeName ) ;
    int setIndex() ;


  private :

	Bool_t delayedInitWrite() ;

	/// utility data
	Bool_t m_verbose ;  /// set the chattiness for debug statements
	  
	/// writing data
	TString m_fileNameWrite ;
	TString m_outputOptions ;
	Bool_t m_initWriteDone ;
    //TFile * m_fileWrite ;
    CompositeEventList * m_celWrite ;
    Long64_t m_eventCounter ;  // Count number of events filled to the TTree so far
    TObjArray m_componentNamesCol ;
    std::vector<TTree*> m_treeCol ;

    /// reading data
	TString m_fileNameRead ;
	//TFile * m_fileRead ;
    CompositeEventList * m_celRead ;
    TObjArray * m_compChainCol ; // List of component TChains for reading
    TChain * m_masterChain ;  // Master TChain for reading
    std::map<TString, TVirtualIndex*> m_chainIndexCol ;
    
 } ;



#endif



