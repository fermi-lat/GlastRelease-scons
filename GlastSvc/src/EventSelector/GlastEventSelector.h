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

/*! GlastEventSelector performes the function of controlling the 
    ApplicationMgr loop. Upgraded on (02/28/2001) to be more general
    then just handeling the IRF input. The input property can now be set 
    to NONE and the loop will run until it hits EvtMax. A new class
    variable was added m_evtMax which is a local copy of the AppMgrs 
    property. You can change the number of events run my changing the
    the EvtMax property of this Svc. Examples of how to do this can be
    found in the GuiSvc.
*/
class GlastEventSelector : virtual public Service, 
                      virtual public IEvtSelector,
                      virtual public IProperty {
public:
  static IService* createClassObject( const std::string& svcname, ISvcLocator* svcloc );
  // IEvtSelector implementation
  virtual StatusCode initialize();
  
  //! The string from Input parameter gets passed in here
  virtual StatusCode setCriteria( const std::string& criteria );
  virtual StatusCode setCriteria( const SelectionCriteria& criteria ); 

  virtual IEvtSelector::Iterator* begin() const;
  virtual IEvtSelector::Iterator* end() const;
  virtual IEvtSelector::Iterator& next(IEvtSelector::Iterator& it) const;
  virtual IEvtSelector::Iterator& previous(IEvtSelector::Iterator& it) const; 
  IOpaqueAddress* reference(const IEvtSelector::Iterator& it) const;

  // IInterface implementation
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);

  // Constructors & Destructors
  GlastEventSelector(  const std::string& name, ISvcLocator* svcloc );
  ~GlastEventSelector();

private:
  // New types
  enum    CriteriaType { 
      IRFFILE, //When the input is a IRF file
      NONE     //When there is no input
  };

  typedef std::list<std::string>   ListName;
  
  //! Used internally to get the MaxEvent from the ApplicationMgr
  StatusCode GlastEventSelector::getMaxEvent();
  //! Parses out single environment variables and well as multiple files
  StatusCode GlastEventSelector::parseStringInList(const std::string& namelist, ListName* flist );
  //! Takes a null string and sticks into it the current file name being worked on.
  StatusCode GlastEventSelector::getFileName(ListName::const_iterator* inputIt, std::string& fName) const;

  std::string           m_criteria;
  std::string           m_jobInput;
  CriteriaType          m_criteriaType;
  GlastEvtIterator      m_evtEnd;
  ListName*             m_inputDataList;
  IGlastDetSvc*         m_detSvc;
  IDataProviderSvc*     m_eventDataSvc;
  GlastEvtIterator*     m_it; 
  int                   m_evtMax;
  

};

#endif  // GLASTCNV_GLASTEVENTSELECTOR_H



















