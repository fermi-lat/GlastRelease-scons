/**  @file InterleaveMapTool.cxx
    @brief implementation of class InterleaveMapTool
    
  $Header$  
*/

#include "IInterleaveMap.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IDataProviderSvc.h"

// access to the tuple
#include "ntupleWriterSvc/INTupleWriterSvc.h"

#include <stdexcept>
#include <sstream>

/** @class InterleaveMapTool
    @brief Provides the connection between input interleave events and the output ntuple
    @author Dan Flath

This provides a mechanisn for relating back to the original intput ntuple row for interleaving with 
generated events. 

*/
//HMK do I take this out?? Aug,2010
static const InterfaceID IID_InterleaveMapTool("InterleaveMapTool", 1 , 0);

class InterleaveMapTool : public AlgTool, virtual public IInterleaveMap
{
public:

    // Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_InterleaveMapTool; }

    // Standard Gaudi Tool constructor
    InterleaveMapTool(const std::string& type, const std::string& name, const IInterface* parent);

    // After building it, destroy it
    ~InterleaveMapTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// Defined in IInterleaveMap
    /// Return a pointer to the map variables
    InterleaveMap* getInterleaveMap() {return m_interleaveMap;}

    /// Save a row
    void saveInterleaveMapRow();

private:

    /// The interleave map
    InterleaveMap*     m_interleaveMap;

    /// access the ntupleWriter service to write out to ROOT ntuples
    INTupleWriterSvc * m_rootTupleSvc;

    StringProperty     m_mapName;     ///< name of the map tree
};

//static ToolFactory<InterleaveMapTool> s_factory;
//const IToolFactory& InterleaveMapToolFactory = s_factory;
DECLARE_TOOL_FACTORY(InterleaveMapTool);

//------------------------------------------------------------------------
InterleaveMapTool::InterleaveMapTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent) :
                                 AlgTool(type, name, parent)
                                 , m_rootTupleSvc(0)
                                 , m_interleaveMap(0)
{
    //Declare the additional interface
    declareInterface<InterleaveMapTool>(this);

    // declare properties with setProperties calls
    declareProperty("MapName", m_mapName="interleave_map");
}
//------------------------------------------------------------------------
InterleaveMapTool::~InterleaveMapTool()
{}

StatusCode InterleaveMapTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Set the properties
    setProperties();

    // Retrieve tuple info an set up
    // get a pointer to RootTupleSvc
    sc = service("RootTupleSvc", m_rootTupleSvc);       
    if( sc.isFailure() ) 
    {
        log << MSG::ERROR << "failed to get the RootTupleSvc" << endreq;
        return sc;
    }

    // Get an InterfaceMap
    m_interleaveMap = new InterleaveMap();

    // Setup a mapping tree between interleave input and output
    m_rootTupleSvc->addItem(m_mapName.value(), "run",    &m_interleaveMap->run);
    m_rootTupleSvc->addItem(m_mapName.value(), "event",  &m_interleaveMap->curEvent);
    m_rootTupleSvc->addItem(m_mapName.value(), "irun",   &m_interleaveMap->irun);
    m_rootTupleSvc->addItem(m_mapName.value(), "ievent", &m_interleaveMap->ievent);
    m_rootTupleSvc->addItem(m_mapName.value(), "itype",   m_interleaveMap->type);
    m_rootTupleSvc->addItem(m_mapName.value(), "value",  &m_interleaveMap->value);

    return sc;
}

void InterleaveMapTool::saveInterleaveMapRow()
{
    // Save this row
    m_rootTupleSvc->saveRow(m_mapName.value());
}

