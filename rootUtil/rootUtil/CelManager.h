/** 
* @file CelManager.h
* @brief The class CelManager is used to set up and handle the composite event list (cel) root files
*
*  $Header$
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#ifndef CelManager_h
#define CelManager_h

#include "rootUtil/CompositeEventList.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <string>
#include <vector>
#include <map>


class CelManager
 {
  public :
  
      CelManager() : m_fileNameWrite(""),m_fileNameRead(""),
          m_fileWrite(0),m_fileRead(0),m_verbose(false), m_eventCounter(0),
          m_compChainCol(0), m_masterChain(0)  
      { };

      ~CelManager() ; 

      /// Writing 
      bool initWrite(const std::string &fileName="cel.root",
          const std::string &options="RECREATE", bool verbose=false);

      bool initOutputFile();

    UInt_t addComponent(const std::string &compName, TTree *t);

	bool fillEvent(); 
    
    bool fillFileAndTreeSet();

    /// Reading
    bool initRead(const std::string &fileName="cel.root", bool verbose=false);
    Long64_t getNumEntries() { return (m_celRead.entries()); };
    Long64_t getEventIndex(const std::string &treeName, Long64_t index);
    TChain* getChainByType(const std::string &treeName);
    int setIndex();


  private :

	TString m_fileNameWrite ; /// The name of the original output file opened for writing
	TString m_fileNameRead;
    std::string m_outputOptions;
    TFile *m_fileWrite, *m_fileRead;
    CompositeEventList m_celWrite, m_celRead;
    std::vector<TTree*> m_treeCol;
	bool m_verbose;  /// set the chattiness for debug statements
    Long64_t m_eventCounter;  // Count number of events filled to the TTree so far
    TObjArray *m_compChainCol; // List of component TChains for reading
    TChain *m_masterChain;  // Master TChain for reading
    std::map<std::string, TVirtualIndex*> m_chainIndexCol;
    
 } ;



#endif



