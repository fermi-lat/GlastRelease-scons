// $Id$
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
#include "GaudiKernel/ClassID.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Tokenizer.h"
#include "GaudiKernel/SvcFactory.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/ConversionSvc.h"
// C++ include files
#include <list>

#include <list>
#include <fstream>

//--------------------------------------------------------------------
// Event Selector 
//--------------------------------------------------------------------

class DbEvtSelector : public Service, virtual public IEvtSelector 
{
public:
  /// IService implementation: Db event selector override
  virtual StatusCode initialize();

  /// IService implementation: Service finalization
  virtual StatusCode finalize();

  // IInterface implementation: query interfaces
  virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);

  /// Create a new event loop context
  /** @param refpCtxt   [IN/OUT]  Reference to pointer to store the context
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode createContext(Context*& refpCtxt) const;

  /** Access last item in the iteration
    * @param refContext [IN/OUT] Reference to the Context object.
    */
  virtual StatusCode last(Context& refContext) const;

  /// Get next iteration item from the event loop context
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode next(Context& refCtxt) const;

  /// Get next iteration item from the event loop context, but skip jump elements
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode next(Context& refCtxt,int jump) const;

  /// Get previous iteration item from the event loop context
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * @param jump      [IN]      Number of events to be skipped
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode previous(Context& refCtxt) const;

  /// Get previous iteration item from the event loop context, but skip jump elements
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * @param jump      [IN]      Number of events to be skipped
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode previous(Context& refCtxt,int jump) const;

  /// Rewind the dataset
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode rewind(Context& refCtxt) const;

  /// Create new Opaque address corresponding to the current record
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode createAddress(const Context& refCtxt, IOpaqueAddress*&) const;

  /// Release existing event iteration context
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * 
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode releaseContext(Context*& refCtxt) const;

  /** Will set a new criteria for the selection of the next list of events and will change
      the state of the context in a way to point to the new list.
      
      @param cr The new criteria string.
      @param c Reference pointer to the Context object.
   */
  virtual StatusCode resetCriteria(const std::string& cr,Context& c)const;

  /// Service Constructor
  DbEvtSelector( const std::string& name, ISvcLocator* svcloc );

  /// Standard destructor
  virtual ~DbEvtSelector();

protected:
  // Data Members
  std::string                   m_selSvcName;
  std::string                   m_cnvSvcName;
  ConversionSvc*                m_cnvSvc;
  CLID                          m_rootCLID;
  std::string                   m_rootName;
  std::string                   m_criteria;
};


//static const SvcFactory<DbEvtSelector> s_DbEvtSelectorFactory;
//const ISvcFactory& DbEvtSelectorFactory = s_DbEvtSelectorFactory;

DECLARE_SERVICE_FACTORY(DbEvtSelector);


/** @class DbContext
 *
 *  @author  M.Frank
 *  @version 1.0
 */
class DbContext : public IEvtSelector::Context {
public:
  typedef std::list<std::string> ListName;
private:
  const DbEvtSelector*          m_pSelector;
  //dbIter<DbObject>*             m_iterator;
  ListName                      m_files;
  std::string                   m_criteria;
  ListName::const_iterator      m_fileIterator;

public:
  /// Standard constructor
  DbContext(const DbEvtSelector* pSelector);
  /// Standard destructor 
  virtual ~DbContext();
  ListName& files() {
    return m_files;
  }
  virtual void* identifier() const {
    return (void*)m_pSelector;
  }
  void setCriteria(const std::string& crit) {
    m_criteria = crit;
  }
  ListName::const_iterator fileIterator()        {
    return m_fileIterator;
  }
  void setFileIterator(ListName::const_iterator new_iter)  {
    m_fileIterator = new_iter;   
  }
//  dbIter<DbObject>* iterator() const {
//    return m_iterator;
//  }
//  void setiterator(dbIter<DbObject>* new_iter)  {
//    if ( m_iterator ) delete m_iterator;
//    m_iterator = new_iter;
//  }
};

DbContext::~DbContext() {
}

DbContext::DbContext(const DbEvtSelector* pSelector)
  : m_pSelector(pSelector) //, m_iterator(0)
{
}

DbEvtSelector::DbEvtSelector(const std::string& name, ISvcLocator* svcloc )  : Service( name, svcloc),
    m_rootCLID(CLID_NULL)
{
  declareProperty("CnvService",    m_selSvcName="DbEvtSelector");
  declareProperty("ConversionSvc", m_cnvSvcName="EventCnvSvc");
}

DbEvtSelector::~DbEvtSelector() {
}

// IInterface::queryInterface
StatusCode DbEvtSelector::queryInterface(const InterfaceID& riid, void** ppvIf) {
  if ( riid == IEvtSelector::interfaceID() )  {
    *ppvIf = (IEvtSelector*)this;
    addRef();
    return SUCCESS;
  }
  return Service::queryInterface( riid, ppvIf );
}

/// IService implementation: Db event selector override
StatusCode DbEvtSelector::initialize()    {
  // Initialize base class
  StatusCode status = Service::initialize();
  MsgStream log(messageService(), name());
  if ( !status.isSuccess() )    {
    log << MSG::ERROR << "Error initializing base class Service!" << endreq;
    return status;
  }
  // Retrieve conversion service handling event iteration
  IService* tmpSvc;
  status = serviceLocator()->service(m_cnvSvcName, tmpSvc, true);
  if( !status.isSuccess() ) {
    log << MSG::ERROR 
        << "Unable to localize the IID from service:" 
        << m_cnvSvcName << endreq;
    return status;
  }
  m_cnvSvc = dynamic_cast<ConversionSvc*>(tmpSvc);

  // Get DataSvc
  IDataManagerSvc* eds = 0;
  status = serviceLocator()->service("EventDataSvc", eds);
  if( !status.isSuccess() ) {
    log << MSG::ERROR 
      << "Unable to localize interface IID_IDataManagerSvc "
      << "from service EventDataSvc"
	    << endreq;
    return status;
  }
  m_rootCLID = eds->rootCLID();
  m_rootName = eds->rootName();
  eds->release();
  log << MSG::INFO << "Selection root:" << m_rootName 
                   << " CLID:" << m_rootCLID << endreq;
  return status;
}

/// IService implementation: Service finalization
StatusCode DbEvtSelector::finalize()    {
  // Initialize base class
  MsgStream log(messageService(), name());
  if ( m_cnvSvc )  {
    m_cnvSvc->release();
    m_cnvSvc = 0;
  }
  StatusCode status = Service::finalize();
  if ( !status.isSuccess() )    {
    log << MSG::ERROR << "Error initializing base class Service!" << endreq;
  }
  return status;
}

StatusCode
DbEvtSelector::createContext(Context*& refpCtxt) const 
{
  refpCtxt = new DbContext(this);
  return StatusCode::SUCCESS;
}

/// Access last item in the iteration
StatusCode DbEvtSelector::last(Context& /* refContext */) const  {
  return StatusCode::FAILURE;
}

StatusCode
DbEvtSelector::next(Context& ctxt) const 
{
  DbContext* pCtxt = dynamic_cast<DbContext*>(&ctxt);
  if ( pCtxt != 0 ) {
      /*
    dbIter<DbObject>* it = pCtxt->iterator();
    if ( !it ) {
      DbContext::ListName::const_iterator fileit = pCtxt->fileIterator();
      if ( fileit != pCtxt->files().end() ) {
        std::string input = *fileit;
        pCtxt->setFileIterator(++fileit);
        it = m_dbMgr->createIterator( input, m_rootName);
        pCtxt->setiterator(it);
        return StatusCode::SUCCESS;
      }
    }
    if ( it ) {
      if ( !it->next() )  {
        dbHandle<DbDataBase> db = it->container().containedIn();
        pCtxt->setiterator(0);
        db.close();
        return next(ctxt);
      }
      return StatusCode::SUCCESS;
    }
  }
  if ( pCtxt ) {
    dbIter<DbObject>* it = pCtxt->iterator();
    if ( it ) {
      dbHandle<DbDataBase> db = it->container().containedIn();
      pCtxt->setiterator(0);
      db.close();
    }
    pCtxt->setFileIterator(pCtxt->files().end());
    */
      return StatusCode::SUCCESS;
  }  
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::next(Context& ctxt,int jump) const 
{
  if ( jump > 0 ) {
    for ( int i = 0; i < jump; ++i ) {
      StatusCode status = next(ctxt);
      if ( !status.isSuccess() ) {
        return status;
      }
    }
    return StatusCode::SUCCESS;
  }
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::previous(Context& /* ctxt */) const 
{
  MsgStream log(msgSvc(), name());
  log << MSG::FATAL 
    << " EventSelector Iterator, operator -- not supported " << endreq;
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::previous(Context& ctxt, int jump) const 
{
  if ( jump > 0 ) {
    for ( int i = 0; i < jump; ++i ) {
      StatusCode status = previous(ctxt);
      if ( !status.isSuccess() ) {
        return status;
      }
    }
    return StatusCode::SUCCESS;
  }
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::rewind(Context& ctxt) const 
{
  DbContext* pCtxt = dynamic_cast<DbContext*>(&ctxt);
  if ( pCtxt ) {
//    dbIter<DbObject>* it = pCtxt->iterator();
//    if ( it ) {
//      dbHandle<DbDataBase> db = it->container().containedIn();
//      pCtxt->setiterator(0);
//      db.close();
//    }
    pCtxt->setFileIterator(pCtxt->files().begin());
    return StatusCode::SUCCESS;
  }
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::createAddress(const Context& ctxt, IOpaqueAddress*& pAddr) const
{
  const DbContext* pctxt = dynamic_cast<const DbContext*>(&ctxt);
  pAddr = 0;
  if ( pctxt ) {
//    dbIter<DbObject>* iter = pctxt->iterator();
//    if ( iter )  {
//      dbHandle<DbObject>& objH = iter->operator->();
//      return m_dbMgr->createAddress(objH, m_rootCLID, pAddr);
//    }
      std::string str("");
      unsigned long temp = 0;
      IAddressCreator* cnvAddr = dynamic_cast<IAddressCreator*>(m_cnvSvc);
      return cnvAddr->createAddress(SICB_StorageType, m_rootCLID, &str, &temp, pAddr);
  }
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::releaseContext(Context*& ctxt) const
{
  DbContext* pCtxt = dynamic_cast<DbContext*>(ctxt);
  if ( pCtxt ) {
//    dbIter<DbObject>* it = pCtxt->iterator();
//    if ( it ) {
//      dbHandle<DbDataBase> db = it->container().containedIn();
//      pCtxt->setiterator(0);
//      db.close();
//    }
    delete pCtxt;
    return StatusCode::SUCCESS;
  }
  return StatusCode::FAILURE;
}

StatusCode 
DbEvtSelector::resetCriteria(const std::string& criteria,
                             Context& context)  const
{
  MsgStream log(messageService(), name());
  DbContext* ctxt = dynamic_cast<DbContext*>(&context);
  std::string db, typ, item, sel, stmt, aut, addr;
  if ( ctxt )  {
    if ( criteria.substr(0,5) == "FILE " )  {
      // The format for the criteria is:
      //        FILE  filename1, filename2 ...
      db = criteria.substr(5);
    }
    else  {
      Tokenizer tok(true);
      tok.analyse(criteria," ","","","=","'","'");
      for(Tokenizer::Items::iterator i=tok.items().begin(); i!=tok.items().end();i++) {
        std::string tmp = (*i).tag().substr(0,3);
        if(tmp=="DAT")  {
          db = (*i).value();
        }
        if(tmp=="OPT")   {
          if((*i).value() != "REA")   {
            log << MSG::ERROR << "Option:\"" << (*i).value() << "\" not valid" << endreq;
            return StatusCode::FAILURE;
          }
        }
        if (tmp=="TYP") {
          typ = (*i).value();
        }
        if(tmp=="ADD")  {
          item = (*i).value();
        }
        if(tmp=="SEL")  {
          sel = (*i).value();
        }
        if(tmp=="FUN")  {
          stmt = (*i).value();
        }
        if(tmp=="AUT")  {
          aut = (*i).value();
        }
        if(tmp=="COL")  {
          addr = (*i).value();
        }
      }
    }
    // It's now time to parse the criteria for the event selection
    // The format for the criteria is:
    //        FILE  filename1, filename2 ...
    //        JOBID number1-number2, number3, ...
    std::string rest = db;
    ctxt->files().clear();
    while(true)  {
      int ipos = rest.find_first_not_of(" ,");
      if (ipos == -1 ) break;
      rest = rest.substr(ipos, std::string::npos);// remove blanks before
      int lpos  = rest.find_first_of(" ,");       // locate next blank
      ctxt->files().push_back( rest.substr(0,lpos )); // insert in list
      if (lpos == -1 ) break;
      rest = rest.substr(lpos, std::string::npos);// get the rest
    }
    ctxt->setFileIterator(ctxt->files().begin());
    ctxt->setCriteria(criteria);
    return StatusCode::SUCCESS;
  }
  log << MSG::ERROR << "Invalid iteration context." << endreq;
  return StatusCode::FAILURE;
}

