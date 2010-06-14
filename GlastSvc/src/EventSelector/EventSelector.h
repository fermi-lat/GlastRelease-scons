// $Header$
//====================================================================
//  GLASTEventSelector.h
//--------------------------------------------------------------------
//
//  Package    : GlastEventSelector  (The Glast Event Selector Package)
//
//  Description: The GlastEventSelector component is able
//               to produce a list of event references given a set of "selection
//               criteria".
//
//
//====================================================================

#ifndef GLASTCNV_GLASTEVENTSELECTOR_H
#define GLASTCNV_GLASTEVENTSELECTOR_H 1

// Include files
#include "Gaudi/Interfaces/IEvtSelector.h"
#include "Gaudi/Interfaces/IProperty.h"
#include "Gaudi/Interfaces/ISvcLocator.h"
#include "Gaudi/Kernel/Service.h"
#include "GlastEvtIterator.h"

#include <list>
#include <vector>
#include <string>
#include <fstream>

// forward declaration
class IGlastDetSvc;
class IDataProviderSvc;


//--------------------------------------------------------------------
// Event Selector 
//--------------------------------------------------------------------

/*! Select events for Glast I/O.
 * Copied from SicBEvent.

*/
class EventSelector : virtual public Service, 
                      virtual public IEvtSelector,
                      virtual public IProperty {
public:
  static IService* createClassObject( const std::string& svcname, ISvcLocator* svcloc );
  // IEvtSelector implementation
  virtual StatusCode initialize();
  virtual StatusCode setCriteria( const std::string& criteria );
  virtual StatusCode setCriteria( const SelectionCriteria& criteria ); 

  virtual IEvtSelector::Iterator* begin() const;
  virtual IEvtSelector::Iterator* end() const;
  virtual IEvtSelector::Iterator& next(IEvtSelector::Iterator& it) const;
  virtual IEvtSelector::Iterator& previous(IEvtSelector::Iterator& it) const; 
  IOpaqueAddress* reference(const IEvtSelector::Iterator& it) const;

  // IProperty implementation
   virtual StatusCode setProperty(const Property& property);
   virtual StatusCode getProperty( Property* property ) const;

  // IInterface implementation
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);

  // Constructors & Destructors
  EventSelector(  const std::string& name, ISvcLocator* svcloc );
  ~EventSelector();

private:
  // New types
  enum    CriteriaType { IRFFILE };
  typedef std::list<std::string>   ListName;
  

  StatusCode EventSelector::parseStringInList(const std::string& namelist, ListName* flist );
  StatusCode EventSelector::getFileName(ListName::const_iterator* inputIt, std::string& fName) const;

  std::string                 m_criteria;
  std::string                 m_jobInput;
  CriteriaType                m_criteriaType;
  GlastEvtIterator            m_evtEnd;
  ListName*                   m_inputDataList;
  IGlastDetSvc*       m_detSvc;
  IDataProviderSvc* m_eventDataSvc;         // Reference to the Event Data Service

};

#endif  // GLASTCNV_GLASTEVENTSELECTOR_H



















