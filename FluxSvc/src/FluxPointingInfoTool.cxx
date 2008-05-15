/**  @file FluxPointingInfoTool.cxx
    @brief implementation of class FluxPointingInfoTool
    
  $Header$  
*/

#include "FluxSvc/IPointingInfo.h"

#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/IDataProviderSvc.h"

#include "ntupleWriterSvc/INTupleWriterSvc.h"
#include "FluxSvc/PointingInfo.h"

/** @class FluxPointingInfoTool
    @brief Manages the Pointing information
    @author Tracy Usher
*/
class FluxPointingInfoTool : public AlgTool, virtual public IPointingInfo
{
public:

    // Standard Gaudi Tool constructor
    FluxPointingInfoTool(const std::string& type, const std::string& name, const IInterface* parent);

    // After building it, destroy it
    ~FluxPointingInfoTool();

    /// @brief Intialization of the tool
    StatusCode initialize();

    /// @brief Finalize method for the tool
    StatusCode finalize();

    //!@brief fill the pointing info for the current orbital status
    virtual void set() {m_pointingInfo.set();}

    //!@brief finish it. 
    virtual void finish(double stop_time, double live) {m_pointingInfo.finish(stop_time, live);}

    //!@brief accessor for time
    virtual double start_time() const {return m_pointingInfo.start_time();}

    //!@brief return TDS object for old scheme
    virtual Event::Exposure* forTDS() const {return m_pointingInfo.forTDS();}

    //!@brief  associate it with the the FT2 tuple
    virtual void setFT2Tuple(const std::string& tname);

    //!@brief  associate it with the the Pt part of the "merit" tuple
    virtual void setPtTuple(const std::string& tname);

    //!@brief Provide ability to read all of PointingInfo data members
    virtual const double get_start()       const {return m_pointingInfo.start;}
    virtual const double get_stop()        const {return m_pointingInfo.stop;}
    virtual const float* get_sc_position() const {return m_pointingInfo.sc_position;}
    virtual const float  get_lat_geo()     const {return m_pointingInfo.lat_geo;}
    virtual const float  get_lon_geo()     const {return m_pointingInfo.lon_geo;}
    virtual const float  get_lat_mag()     const {return m_pointingInfo.lat_mag;}
    virtual const float  get_rad_geo()     const {return m_pointingInfo.rad_geo;}
    virtual const float  get_ra_zenith()   const {return m_pointingInfo.ra_zenith;}
    virtual const float  get_dec_zenith()  const {return m_pointingInfo.dec_zenith;}
    virtual const float  get_ra_scz()      const {return m_pointingInfo.ra_scz;}
    virtual const float  get_dec_scz()     const {return m_pointingInfo.dec_scz;}
    virtual const float  get_ra_scx()      const {return m_pointingInfo.ra_scx;}
    virtual const float  get_dec_scx()     const {return m_pointingInfo.dec_scx;}
    virtual const float  get_in_saa()      const {return m_pointingInfo.in_saa;}
    virtual const float  get_livetime()    const {return m_pointingInfo.livetime;}
    virtual const float  get_L()           const {return m_pointingInfo.L;}
    virtual const float  get_B()           const {return m_pointingInfo.B;}
    virtual const float  get_zenith_scz()  const {return m_pointingInfo.zenith_scz;}

private:
    /// Pointer to the Gaudi data provider service
    INTupleWriterSvc* m_rootTupleSvc;

    /// The root tree for the pointing info
    StringProperty    m_info_tree;
    StringProperty    m_history_tree;

    /// The all important PointingInfo class
    PointingInfo      m_pointingInfo;
};

static ToolFactory<FluxPointingInfoTool> s_factory;
const IToolFactory& FluxPointingInfoToolFactory = s_factory;
//------------------------------------------------------------------------

FluxPointingInfoTool::FluxPointingInfoTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent) :
                                 AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<IPointingInfo>(this);

    // declare properties with setProperties calls
    declareProperty("Pointing_Info_Tree_name",  m_info_tree    = "MeritTuple");
    declareProperty("Pointing_History_Tree",    m_history_tree = "pointing_history"); //doesn't work???

    return;
}
//------------------------------------------------------------------------
FluxPointingInfoTool::~FluxPointingInfoTool()
{
}

StatusCode FluxPointingInfoTool::initialize()
{
    StatusCode sc   = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    // Set the properties
    setProperties();

    // get a pointer to RootTupleSvc, use only if available 
    if( (service("RootTupleSvc", m_rootTupleSvc, true) ). isFailure() ) 
    {
        log << MSG::WARNING << " RootTupleSvc is not available, will not write Pt tuple" << endreq;
        m_rootTupleSvc=0;
    }

    // I can't imagine why the following would ever be the case but we preserve it for posterity...
    if( !m_info_tree.value().empty() ) 
    {     
        m_pointingInfo.setPtTuple(m_rootTupleSvc, m_info_tree.value());
    }

    return sc;
}

StatusCode FluxPointingInfoTool::finalize ()
{
    StatusCode  status = StatusCode::SUCCESS;
    
    return status;
}

void FluxPointingInfoTool::setFT2Tuple(const std::string& tname) 
{
    if (tname != m_history_tree.value() && !m_history_tree.value().empty())
    {
        MsgStream log(msgSvc(), name());

        log << MSG::WARNING << "Input history tree name (" << tname << ") does not match default name " 
            << m_history_tree.value() << "!" << endreq;
    }

    m_pointingInfo.setFT2Tuple(m_rootTupleSvc, tname);
}

void FluxPointingInfoTool::setPtTuple(const std::string& tname) 
{
    if (tname != m_info_tree.value() && !m_info_tree.value().empty())
    {
        MsgStream log(msgSvc(), name());

        log << MSG::WARNING << "Input info tree name (" << tname << ") does not match default name " 
            << m_info_tree.value() << "!" << endreq;
    }

    m_pointingInfo.setPtTuple(m_rootTupleSvc, tname);
}
