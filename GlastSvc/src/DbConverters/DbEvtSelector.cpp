//====================================================================
//	DbEvtSelector.cpp
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

// Include files
#include "Gaudi/Kernel/StatusCode.h"
#include "Gaudi/Kernel/SvcFactory.h"
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Interfaces/IAddressCreator.h"

#include "DbCnv/DbOOMs.h"
#include "DbCnv/DbObject.h"
#include "DbCnv/IDataBaseMgr.h"
#include "DbCnv/DbIter.h"
#include "src/DbCnv/DbAddress.h"
#include "src/Services/DbEvtSelector.h"

extern const CLID& CLID_Event;
extern unsigned char DBOOMS_StorageType;

static const SvcFactory<DbEvtSelector>    s_DbEvtSelectorFactory;
const ISvcFactory& DbEvtSelectorFactory = s_DbEvtSelectorFactory;

class DbEventIterator: public IEvtSelector::Iterator {
private:
  dbIter<DbObject>*                        m_iterator;
  const IEvtSelector*                      m_selector;
  DbEvtSelector::ListName::const_iterator  m_files;
public:
  friend class DbEvtSelector;
  DbEventIterator( const IEvtSelector* selector, dbIter<DbObject>* it, DbEvtSelector::ListName::const_iterator fit)
  : m_iterator(it), m_selector(selector), m_files(fit)
  {
  }
  virtual ~DbEventIterator() { 
    if ( 0 != m_iterator )    {
      delete m_iterator;
    }
  }
  virtual void release()                          { 
    delete this;  
  }
  dbIter<DbObject>* iter()                        {
    return m_iterator;
  }
  void replaceObjIterator(dbIter<DbObject>* new_iter)   {
    if ( m_iterator ) delete m_iterator;
    m_iterator = new_iter;
  }
  DbEvtSelector::ListName::const_iterator fileIterator()    {
    return m_files;
  }
  void replaceFileIterator(DbEvtSelector::ListName::const_iterator new_iter)    {
    m_files = new_iter;
  }
  virtual IOpaqueAddress* operator*() const       { 
    return m_selector->reference(*this);
  }
  virtual IOpaqueAddress* operator->() const      { 
    return m_selector->reference(*this);
  }
  virtual IEvtSelector::Iterator& operator++()    { 
    return m_selector->next(*this);
  }
  virtual IEvtSelector::Iterator& operator--()    { 
    return m_selector->previous(*this);
  }
  virtual IEvtSelector::Iterator& operator++(int) { 
    return m_selector->next(*this);
  }
  virtual IEvtSelector::Iterator& operator--(int) { 
    return m_selector->previous(*this);
  }
  virtual bool operator==(const IEvtSelector::Iterator& it) const {
    const DbEventIterator* pIt = dynamic_cast<const DbEventIterator*>(&it);
    if ( m_iterator && pIt->m_iterator )    {
      return (*m_iterator == *(pIt->m_iterator)) && (m_files == pIt->m_files);
    }
    return (m_iterator == pIt->m_iterator) && (m_files == pIt->m_files);
  }
  virtual bool operator!=(const IEvtSelector::Iterator& it) const {
    return !(operator ==(it));
  }
};


DbEvtSelector::DbEvtSelector(const std::string& name, ISvcLocator* svcloc ) 
  : Service( name, svcloc),
    m_dbManager(0),
    m_emptyIterator(0)
{
  declareProperty("CnvService", m_cnvService);
  m_filenameList  = new ListName;
  m_emptyIterator = new DbEventIterator(this, 0, m_filenameList->end());
}

DbEvtSelector::~DbEvtSelector() {
  delete m_filenameList;
  delete m_emptyIterator;
}

/// IService implementation: Db event selector override
StatusCode DbEvtSelector::initialize()    {
  // Initialize base class
  StatusCode status = Service::initialize();
  MsgStream logger(messageService(), name());
  if ( !status.isSuccess() )    {
    logger << "Error initializing base class Service!" << endreq;
    return status;
  }
  // Get properties from job options service
  status = setProperties();
  if ( !status.isSuccess() )    {
    logger << "Error retrieving properties!" << endreq;
    return status;
  }
  // Retrive conversion service handling event iteration
  status = serviceLocator()->getService(m_cnvService, IID_IDataBaseMgr, (IInterface*&)m_dbManager);
  if( !status.isSuccess() ) {
    logger << "Unable to localize interface IID_IDataBaseMgr from service:" << m_cnvService << endreq;
    return status;
  }
  return status;
}

// IEvtSelector::setCriteria
StatusCode DbEvtSelector::setCriteria( const std::string& criteria ) {

  StatusCode sc = StatusCode::FAILURE;
  MsgStream logger(messageService(), name());
  m_criteria = criteria;
  // It's now time to parse the criteria for the event selection
  // The format for the criteria is:
  //        FILE  filename1, filename2 ...
  //        JOBID number1-number2, number3, ...
  std::string rest = criteria;
  std::string type;
  int ipos  = rest.find_first_not_of(' ');
  rest      = rest.substr(ipos, -1);            // remove blanks before
  int lpos  = rest.find_first_of(' ');          // locate next blank
  type      = rest.substr(0,lpos );
  rest      = rest.substr(lpos, -1);
  if( type == "FILE" )  {
    m_filenameList->erase(m_filenameList->begin(), m_filenameList->end());
    sc = parseFileCriteria( rest, m_filenameList );    
    if( sc.isSuccess() ) {
      return sc;
    }
  }
  logger << MSG::ERROR << "Invalid Event Selection Criteria: " << criteria << endreq;
  return sc;
}

StatusCode DbEvtSelector::parseFileCriteria( const std::string& namelist, ListName* filelist ) {
  std::string rest = namelist;
  while(true) {
    int ipos = rest.find_first_not_of(" ,");
    if (ipos == -1 ) break;
    rest = rest.substr(ipos, -1);               // remove blanks before
    int lpos  = rest.find_first_of(" ,");       // locate next blank
    filelist->push_back( rest.substr(0,lpos )); // insert in list
    if (lpos == -1 ) break;
    rest = rest.substr(lpos, -1);               // get the rest
  }
  return StatusCode::SUCCESS;
}

// IEvtSelector::setCriteria
StatusCode DbEvtSelector::setCriteria( const SelectionCriteria& /* criteria */ ) {
  return StatusCode::SUCCESS;
}

// IEvtSelector::begin()
IEvtSelector::Iterator* DbEvtSelector::begin() const {
  return begin(1);
}

// IEvtSelector::begin()
IEvtSelector::Iterator* DbEvtSelector::begin(int iPile) const   {
  DbEventIterator  *pIt  = m_emptyIterator;
  if ( iPile > 1 )    {
    MsgStream log(msgSvc(), name());
    log << MSG::FATAL << " EventSelector::begin pileup is not yet not supported " << endreq;
  }
  else  {
    ListName::const_iterator fileit = m_filenameList->begin();
    if ( fileit != m_filenameList->end() )    {
      dbIter<DbObject>* iter = m_dbManager->createIterator( *fileit, "/Event");
      pIt = new DbEventIterator(this, iter, fileit);
    }
  }
  return pIt;
}

IEvtSelector::Iterator& DbEvtSelector::next(IEvtSelector::Iterator& it) const   {
  return next(it,1);
}

IEvtSelector::Iterator& DbEvtSelector::next(IEvtSelector::Iterator& it, int iPile) const  {
  DbEventIterator  *pIt  = dynamic_cast<DbEventIterator*>(&it);
  dbIter<DbObject> *iter = pIt->iter(), *new_iter = 0;
  if ( iter && iPile >= 1 )   {
    if( !iter->next() )    {
      ListName::const_iterator fileit = pIt->fileIterator();
      fileit++;
      if( fileit != m_filenameList->end() )   {
        new_iter = m_dbManager->createIterator ( *(fileit), "/Event" );
      }
      pIt->replaceObjIterator(new_iter);
      pIt->replaceFileIterator(fileit);
    }
  }
  else    {
    pIt->replaceObjIterator(new_iter);
    pIt->replaceFileIterator(m_filenameList->end());
  }
  return it;
}

// IEvtSelector::end()
IEvtSelector::Iterator* DbEvtSelector::end() const {
  return m_emptyIterator;
}

// IEvtSelector::previous
IEvtSelector::Iterator& DbEvtSelector::previous(IEvtSelector::Iterator& it) const {
  MsgStream log(msgSvc(), name());
  log << MSG::FATAL << " EventSelector Iterator, operator -- not supported " << endreq;
  return it;
}

//  IEvtSelector::reference
IOpaqueAddress* DbEvtSelector::reference(const IEvtSelector::Iterator& it) const {
  const DbEventIterator* cpIt = dynamic_cast<const DbEventIterator*>(&it);
  DbEventIterator* pIt = const_cast<DbEventIterator*>(cpIt);
  dbIter<DbObject>* iter = pIt->iter();
  if ( 0 != iter )    {
    dbHandle<DbObject>& objH = iter->operator->();
    if ( objH.isValid() )   {
      DbAddress* pAddress = new DbAddress(objH, CLID_Event);
      return pAddress;
    }
  }
  return 0;
}

// IInterface::queryInterface
StatusCode DbEvtSelector::queryInterface(const IID& riid, void** ppvInterface)  {
  if ( riid == IID_IEvtSelector )  {
    *ppvInterface = (IEvtSelector*)this;
  }
  else if ( riid == IID_IProperty )  {
    *ppvInterface = (IProperty*)this;
  }
  else   {
    return Service::queryInterface( riid, ppvInterface );
  }
  addRef();
  return SUCCESS;
}
