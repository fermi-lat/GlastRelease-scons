#ifndef GLASTCNV_GLASTEVENTSELECTOR_H
#define GLASTCNV_GLASTEVENTSELECTOR_H 1

#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Service.h"
#include "GlastEvtIterator.h"

#include <list>
#include <vector>
#include <string>
#include <fstream>

class IAddressCreator;
class IGlastDetSvc;
class IDataProviderSvc;

/** @class GlastEventSelector
* @brief GlastEventSelector performes the function of controlling the 
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
    
    GlastEventSelector(  const std::string& name, ISvcLocator* svcloc );
    ~GlastEventSelector();
    
private:
    enum CriteriaType { 
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
    StringProperty         m_jobInput;
    CriteriaType          m_criteriaType;
    GlastEvtIterator      m_evtEnd;
    ListName*             m_inputDataList;
    IGlastDetSvc*         m_detSvc;
    IDataProviderSvc*     m_eventDataSvc;
    GlastEvtIterator*     m_it; 
    IntegerProperty           m_evtMax;
    
    IAddressCreator*      m_addrCreator;
    CLID                  m_rootCLID;
};

#endif  // GLASTCNV_GLASTEVENTSELECTOR_H









