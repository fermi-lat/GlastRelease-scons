
/** 
* @file CelManager.cxx
* @brief definition of the class CelManager
*
* File: $Header$
* Authors:
*   HK, Heather Kelly, heather@lheapop.gsfc.nasa.gov
*   DC, David Chamont, LLR, chamont@poly.in2p3.fr
*/

#include <rootUtil/CelManager.h>
#include <rootUtil/CelEventComponent.h>
#include <TSystem.h>
#include <Riostream.h>


CelManager::CelManager()
 : m_verbose(false), 
   m_fileNameWrite(""), m_fileWrite(0), m_eventCounter(0),
   m_fileNameRead(""),  m_fileRead(0), m_compChainCol(0), m_masterChain(0)  
 {}

CelManager::~CelManager() 
 {
  m_treeCol.clear() ;
  if (m_compChainCol) m_compChainCol->Delete() ;
  if (m_masterChain) delete m_masterChain ;
 }

// Writing Methods 

Bool_t CelManager::initWrite( const TString & fileName, const TString & options, bool verbose)
 {
  Bool_t stat = kTRUE ;
  m_fileNameWrite = fileName ;
  if (gSystem->ExpandPathName(m_fileNameWrite)==kTRUE)
   {
    std::cout << "Failed to expand env variable in filename" << std::endl;
    return kFALSE ;
   }
  m_outputOptions = options ;
  m_verbose = verbose ;
  return stat;
 }

Bool_t CelManager::initOutputFile()
 {
  Bool_t stat = kTRUE ;
  m_fileWrite = m_celWrite.makeFile(m_fileNameWrite,m_outputOptions.Data()) ;
  if (!m_fileWrite->IsOpen())
   {
    stat = kFALSE ;
    std::cout << "CelManager Error opening ROOT file " << m_fileNameWrite << std::endl;
    return stat ;
   }
  return stat ;
 }


UInt_t CelManager::addComponent( const TString & compName, TTree * t )
 {
  unsigned int ret = m_celWrite.addComponent(compName) ;
  m_treeCol.push_back(t) ;
  return ret ;
 }

Bool_t CelManager::fillEvent()
 {
  Bool_t stat = kTRUE ;
  try
   {
    TDirectory * saveDir = gDirectory ;
    // Need to call CompositeEventList::makeFile after the AddComponent calls
    if ((!m_fileWrite)&&(!initOutputFile()))
     { throw ; }
    if (m_fileWrite->TestBits(TFile::kWriteError))
     { throw ; }
    m_fileWrite->cd() ;
    std::vector<TTree*>::iterator treeItr ;
    Long64_t numBytes ;
    for ( treeItr=m_treeCol.begin() ; treeItr != m_treeCol.end() ; treeItr++ )
     { numBytes = (*treeItr)->LoadTree(m_eventCounter) ; }
    m_celWrite.fillEvent(m_treeCol) ;
    saveDir->cd() ;
    ++m_eventCounter ;
   }
  catch(...)
   {
    std::cerr << "Error filling Meta ROOT file" << std::endl ;
    std::cerr.flush() ;
    throw ;
   }
  return stat ;
 }

Bool_t CelManager::fillFileAndTreeSet()
 {
  Bool_t stat = kTRUE ;
  try
   {
    TDirectory * saveDir = gDirectory ;
    m_fileWrite->cd() ;
    m_celWrite.fillFileAndTreeSet() ;
    m_fileWrite->Write(0,TObject::kWriteDelete) ;
    m_fileWrite->Close() ;
    saveDir->cd() ;
   }
  catch(...)
   {
    std::cerr << "Failed final write to meta ROOT file" << std::endl ; 
    std::cerr.flush(); 
    throw ;
   }
  return stat ;
 }



/// Reading Methods

Bool_t CelManager::initRead( const TString & fileName, Bool_t verbose )
 {
  Bool_t stat = kTRUE ;
  m_fileNameRead = fileName ;
  if (gSystem->ExpandPathName(m_fileNameRead)==kTRUE)
   {
    std::cout << "Failed to expand env variable in filename" << std::endl ;
    return kFALSE ;
   }
  m_verbose = verbose ;
  m_fileRead = m_celRead.openFile(m_fileNameRead) ;
  if (!m_fileRead->IsOpen())
   {
    stat = kFALSE ;
    std::cout << "CelManager Error opening ROOT file " << m_fileNameRead << std::endl ;
    return stat ;
   }

  // Set up our master TChain and the component TChains
  if (!m_compChainCol) m_compChainCol = new TObjArray();
  m_masterChain = m_celRead.buildLinks(m_compChainCol);
  TIter itr(m_compChainCol);
  TChain *curChain;
  while ( (curChain = (TChain*)itr.Next()) )
   {
    curChain->SetBranchStatus("*", 1) ;
    m_chainIndexCol[curChain->GetName()] = curChain->GetTreeIndex() ;
   }
  return stat ;
 }

Long64_t CelManager::getEventIndex( const TString & treeName, Long64_t index)
 {
  // make sure the cel index is the active TVirtualIndex
  setIndex() ;
  Long64_t retVal = m_masterChain->LoadTree(index) ;  // Causes all component TChains to be loaded correctly via CelIndex
  if (retVal<0) return retVal ;
  TChain * compChain = getChainByType(treeName) ;
  if (!compChain) return -1 ;
  return (compChain->GetReadEntry()) ;
 }

int CelManager::setIndex()
 {
  int numSet=0 ;
  std::map<TString, TVirtualIndex*>::iterator mapIt ;
  for ( mapIt=m_chainIndexCol.begin() ; mapIt != m_chainIndexCol.end() ; mapIt++ )
   {
    TChain * curChain = getChainByType(mapIt->first) ;
    if (curChain)
     {
      curChain->SetTreeIndex(mapIt->second) ;
      numSet++ ;
     }
   }
  return numSet ;
 }

TChain * CelManager::getChainByType( const TString & treeName )
 {
  TIter itr(m_compChainCol) ;
  TChain * curChain ;
  while ((curChain = (TChain*)itr.Next()))
   { if (curChain->GetName() == treeName) return curChain ; }
  return 0 ;
 }



