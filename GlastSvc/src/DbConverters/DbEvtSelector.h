//====================================================================
//	DbEvtSelector.h
//--------------------------------------------------------------------
//
//	Package    : DbEvtSelector  (The LHCb Event Selector Package)
//
//  Description: The DbEvtSelector component is able
//               to produce a list of event references given a set of "selection
//               criteria".
//
//	Author     : P. Mato
//  Created    : 4/01/99
//	Changes    : 
//
//
//====================================================================
// $Header$

#ifndef DBCNV_DBEVTSELECTOR_H
#define DBCNV_DBEVTSELECTOR_H 1

// Include files
#include "Gaudi/Interfaces/IEvtSelector.h"
#include "Gaudi/Interfaces/IProperty.h"
#include "Gaudi/Interfaces/ISvcLocator.h"
#include "Gaudi/Kernel/Service.h"

#include <list>
#include <fstream>

// Forward declarations
class IAddressCreator;
class DbEventIterator;

#ifdef WIN32
// Disable warning C4786: identifier was truncated to '255' characters in the debug information
#pragma warning ( disable : 4786 )
#endif


//--------------------------------------------------------------------
// Event Selector 
//--------------------------------------------------------------------

class DbEvtSelector : public Service, virtual public IEvtSelector {
public:
  // IEvtSelector implementation
  virtual StatusCode setCriteria( const std::string& criteria );
  virtual StatusCode setCriteria( const SelectionCriteria& criteria );  
  virtual IEvtSelector::Iterator* begin() const;
  virtual IEvtSelector::Iterator* begin(int iPile) const;
  virtual IEvtSelector::Iterator* end() const;
  virtual IEvtSelector::Iterator& next(IEvtSelector::Iterator& it) const;
  virtual IEvtSelector::Iterator& next(IEvtSelector::Iterator& it, int iPile) const;
  virtual IEvtSelector::Iterator& previous(IEvtSelector::Iterator& it) const; 
  IOpaqueAddress* reference(const IEvtSelector::Iterator& it) const;

  /// IService implementation: Db event selector override
  virtual StatusCode initialize();
  // IInterface implementation: query interfaces
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);

  // Constructors & Destructors
  DbEvtSelector( const std::string& name, ISvcLocator* svcloc );
  ~DbEvtSelector();

private:
  // New types
  typedef std::list<std::string>   ListName;

  StatusCode parseFileCriteria( const std::string& namelist, ListName* flist );

  // Data Members
  std::string               m_criteria;
  std::string               m_cnvService;
  ListName*                 m_filenameList;
  IDataBaseMgr*             m_dbManager;
  mutable DbEventIterator*  m_emptyIterator;
};
#endif  // DBCNV_DBEVTSELECTOR_H
