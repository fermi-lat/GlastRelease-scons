#ifndef LdfEventSelector_H
#define LdfEventSelector_H 1

#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Service.h"
#include "src/LdfEvtIterator.h"

#include "ldfReader/LdfParser.h"

#include <list>
#include <vector>
#include <string>
#include <fstream>

class IAddressCreator;
class IDataProviderSvc;

/** @class LdfEventSelector
* @brief LdfEventSelector performes the function of controlling the 
* ApplicationMgr loop. 
*
* Upgraded on (02/28/2001) to be more general then just handling the 
* IRF input. The input property can now be set to NONE and the loop 
* will run until it hits EvtMax. A new class variable was added m_evtMax 
* which is a local copy of the AppMgrs property. You can change the 
* number of events run my changing the EvtMax property of this Svc. 
* Examples of how to do this can be found in the GuiSvc.
*
* $Header$
*/
class LdfEventSelector : virtual public Service, 
                         virtual public IEvtSelector,
                         virtual public IProperty {
public:
  static IService* createClassObject( const std::string& svcname, ISvcLocator* svcloc );
  // IEvtSelector implementation
  virtual StatusCode initialize();

  // All the information is in instance variables so don't really need
  // to pass it in, but have to have a setCriteria with string arg.
  //! The string from Input parameter gets passed in here
  virtual StatusCode setCriteria(const std::string& storageType);

  // Also needed for IEvtSelector interface
  virtual StatusCode setCriteria( const SelectionCriteria& criteria ); 
    
  virtual IEvtSelector::Iterator* begin() const;
  virtual IEvtSelector::Iterator* end() const;
  virtual IEvtSelector::Iterator& next(IEvtSelector::Iterator& it) const;
  virtual IEvtSelector::Iterator& previous(IEvtSelector::Iterator& it) const; 
  IOpaqueAddress* reference(const IEvtSelector::Iterator& it) const;
    
  // IInterface implementation
  virtual StatusCode queryInterface(const IID& riid, void** ppvInterface);
    
  LdfEventSelector(const std::string& name, ISvcLocator* svcloc );
  ~LdfEventSelector();
    
private:
  enum CriteriaType { 
    EBFFILE, //When the input is a raw file
    EBFFITS, //When the input is fitsified raw ebf
    NONE     //When there is no input
  };
    
  typedef std::list<std::string>   ListName;
    
  //! Used internally to get the MaxEvent from the ApplicationMgr
  StatusCode LdfEventSelector::getMaxEvent();

  /// Job options parameter. Options are EBFFILE, EBFFITS
  StringProperty        m_storageType;

  /// Job options parameter. List of input data files
  StringArrayProperty   m_inputList;

  /// same info as m_storageType, but as enum
  CriteriaType          m_criteriaType;

  LdfEvtIterator        m_evtEnd;
  ListName*             m_inputDataList;
  IDataProviderSvc*     m_eventDataSvc;
  LdfEvtIterator*       m_it; 
  IntegerProperty       m_evtMax;
  std::string           m_instrument;  // for now should be EM or LAT
  std::string           m_fileName;
  ldfReader::LdfParser*            m_ebfParser;
  IntegerProperty       m_ebfDebugLevel;  // If 0 (the default) no debug.
  IAddressCreator*      m_addrCreator;
  CLID                  m_rootCLID;
};

#endif  // LdfEventSelector_H
