#ifndef LdfEventSelector_H
#define LdfEventSelector_H 1

#include "GaudiKernel/IEvtSelector.h"
#include "GaudiKernel/IProperty.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/Service.h"
//#include "src/LdfEvtIterator.h"

#include "ldfReader/EbfParser.h"

#include <list>
#include <vector>
#include <string>
#include <fstream>

class IAddressCreator;
class IDataProviderSvc;

class LdfSelectorContext;

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

 /// HMK Access "real" selector  try this out HMK :)
  const IEvtSelector* selector()  const  {
    return this;
  }

/// NEW from GaudiSvc's EventSelector.h that was required by IEvtSelector

  void printEvtInfo(const LdfSelectorContext* iter) const;

  /// Create a new event loop context
  /** @param refpCtxt   [IN/OUT]  Reference to pointer to store the context
    *
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode createContext(Context*& refpCtxt) const;

  /// Get next iteration item from the event loop context
  /** @param refCtxt   [IN/OUT]  Reference to the context
    *
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode next(Context& refCtxt) const;

  /// Get next iteration item from the event loop context, but skip jump 
  /// elements
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

  /// Get previous iteration item from the event loop context, but skip jump 
  /// elements
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
    * @param refpAddr  [OUT]     Reference to address pointer
    *
    * @return StatusCode indicating success or failure
    */

  /// Create new Opaque address corresponding to the current record
  /** @param refCtxt   [IN/OUT]  Reference to the context
    * @param refpAddr  [OUT]     Reference to address pointer
    *
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode createAddress(const Context& refCtxt, IOpaqueAddress*& refpAddr) const;

  /// Release existing event iteration context
  /** @param refCtxt   [IN/OUT]  Reference to the context
    *
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode releaseContext(Context*& refCtxt) const;

  /** Will set a new criteria for the selection of the next list of events and 
    *  will change
    * the state of the context in a way to point to the new list.
    *
    * @param cr The new criteria string.
    * @param c Reference pointer to the Context object.
    *
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode resetCriteria(const std::string& cr,Context& c)const;

  /** Access last item in the iteration
    * @param c Reference to the Context object.
    *
    * @return StatusCode indicating success or failure
    */

  /** Access last item in the iteration
    * @param c Reference to the Context object.
    *
    * @return StatusCode indicating success or failure
    */
  virtual StatusCode last(Context& c) const;


/// OLD
//  static IService* createClassObject( const std::string& svcname, ISvcLocator* svcloc );


  // IInterface implementation
  virtual StatusCode queryInterface(const InterfaceID& riid, void** ppvInterface);
    
  LdfEventSelector(const std::string& name, ISvcLocator* svcloc );
  virtual ~LdfEventSelector();

  // IEvtSelector implementation
  virtual StatusCode initialize();

  // All the information is in instance variables so don't really need
  // to pass it in, but have to have a setCriteria with string arg.
  //! The string from Input parameter gets passed in here
  virtual StatusCode setCriteria(const std::string& storageType);

  // Also needed for IEvtSelector interface
  //virtual StatusCode setCriteria( const SelectionCriteria& criteria ); 
    
  //virtual IEvtSelector::Iterator* begin() const;
 // virtual IEvtSelector::Iterator* end() const;
 // virtual IEvtSelector::Iterator& next(IEvtSelector::Iterator& it) const;
 // virtual IEvtSelector::Iterator& previous(IEvtSelector::Iterator& it) const; 
 // IOpaqueAddress* reference(const IEvtSelector::Iterator& it) const;
    
// From gaudi v18r1 EventSelector in GaudiSvc
 //HMK typedef std::vector<EventSelectorDataStream*>  Streams;
  //HMKtypedef std::vector<std::string>               StreamSpecs;
  //HMKtypedef std::vector<StringProperty>            Properties;
    
private:

  void checkForSkippedEvents(unsigned long long& counter, bool &lastEventFlag);

  enum CriteriaType { 
    EBFFILE,
    EBFFITS,
    LDFFILE, //When the input is a raw file
    LDFFITS, //When the input is fitsified raw ldf
    LSFFILE, // New LSF file 
    CCSDSFILE, // Access to packet archive setup through XML file
    NONE     //When there is no input
  };
    
 // typedef std::list<std::string>   ListName;
    
  //! Used internally to get the MaxEvent from the ApplicationMgr
  StatusCode LdfEventSelector::getMaxEvent();

  /// Job options parameter. Options are LDFFILE, LDFFITS, LSFFILE
  StringProperty        m_storageType;

  /// Job options parameter. List of input data files
  //StringArrayProperty   m_inputList;

  /// same info as m_storageType, but as enum
  CriteriaType          m_criteriaType;

  //LdfEvtIterator        m_evtEnd;
  //ListName*             m_inputDataList;
  IDataProviderSvc*     m_eventDataSvc;
  //LdfEvtIterator*       m_it; 
//  IntegerProperty       m_evtMax;  // see below
  std::string           m_instrument;  // for now should be EM or LAT
  std::string           m_fileName;
  ldfReader::EbfParser* m_ebfParser;
  IntegerProperty       m_ebfDebugLevel;  // If 0 (the default) no debug.
  IntegerProperty       m_sweepSearch; // If 0, skip search
  IAddressCreator*      m_addrCreator;
  CLID                  m_rootCLID;
  IntegerProperty       m_gemCheck; // If 0 (default) no GEM condition summary
                                    // check, other check for zero and print
                                    // log message

  IntegerProperty       m_startEventIndex; // used if we want to skip to 
                                          // Event N in processing
  IntegerProperty       m_startEventNumber; // used if we want to skip to event
                                            // with event number N as stored
                                            // in LDF

  std::string           m_acdRemap; // For beamtest 2006

  IntegerProperty       m_ignoreSegFault;


  /// Default 0 (false)
  /// Flag that determines what quantity from the evt file to use for the
  /// runId
  IntegerProperty       m_oldStyleRunId;

  IntegerProperty       m_server;  // For socket commections
  IntegerProperty       m_socket; // Boolean JO parameter to denote use of 
                                  // sockets rather than file input

// From gaudi v18r1 EventSelector in GaudiSvc
  //IIncidentSvc*          m_incidentSvc;
  /// Reconfigure occurred
  bool                  m_reconfigure;
  /// Input stream specifiers (for job options)
  //HMKStreamSpecs           m_streamSpecs;
  /// Input stream specifiers (last used)
  //HMKStreamSpecs           m_streamSpecsLast;
  /// Input streams
  //HMKKStreams               m_streams;
  /// Input stream counter (0..oo, monotonely increasing)
  //HMKint                   m_streamCount;
  /// First event to be processed
  int                   m_firstEvent;
  /// Maximum number of events to be processed
  int                   m_evtMax; // Obsolete parameter
  /// Printout frequency
  int                   m_evtPrintFrequency;
  /// Potential list of gem ids to skip
  UnsignedLongLongArrayProperty m_gemIdSkipList;
  std::vector<unsigned long long> m_gemIdSkipVec;
  /// Potential list of event indices to skip
  UnsignedLongLongArrayProperty m_eventIndexSkipList;
  std::vector<unsigned long long> m_eventIndexSkipVec;
  /// JO parameter denotes whether to insert an empty event if we skip events
  bool m_insertEmptyEvent;
  std::string m_eventIndicesToSkipEnv;
  bool m_printCounter;
};

#endif  // LdfEventSelector_H
