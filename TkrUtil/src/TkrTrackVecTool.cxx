/**
* @class TkrTrackVecTool
*
* @brief Tool for setting the ghost attributes of hits and tracks
*
* @author The Tracking Software Group
*
* $Header$
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "TkrUtil/ITkrTrackVecTool.h"


class TkrTrackVecTool : public AlgTool, virtual public ITkrTrackVecTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrTrackVecTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrTrackVecTool() {}

    /// @brief 

    StatusCode initialize();

    std::vector<Event::TkrTrack*> getTrackVec();

private:

    DataSvc*               m_dataSvc;
};

//static ToolFactory<TkrTrackVecTool> s_factory;
//const IToolFactory& TkrTrackVecToolFactory = s_factory;
DECLARE_TOOL_FACTORY(TkrTrackVecTool);

TkrTrackVecTool::TkrTrackVecTool(const std::string& type, 
                                 const std::string& name, 
                                 const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrTrackVecTool>(this);

    return;
}

StatusCode TkrTrackVecTool::initialize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        return sc;
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    return sc;
}

std::vector<Event::TkrTrack*> TkrTrackVecTool::getTrackVec()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    std::vector<Event::TkrTrack*> trackVec(0);

    SmartDataPtr<Event::TkrTrackCol> 
        trackCol(m_dataSvc, EventModel::TkrRecon::TkrTrackCol);
    int size;

    SmartDataPtr<Event::TkrTrackCol> 
        cRTrackCol(m_dataSvc, EventModel::TkrRecon::TkrCRTrackCol);

    Event::TkrTrackColConPtr tcolIter;

    int trackCount = 0;
    if(trackCol) {
        size = trackCol->size();
        tcolIter = trackCol->begin();
        for(; tcolIter!=trackCol->end(); ++tcolIter,++trackCount) {
            Event::TkrTrack* track = *tcolIter;
            trackVec.push_back(track);
        }
    }

    //std::cout << trackCol->size() << " normal tracks, trackCount = " << trackCount<< std::endl;

    if(cRTrackCol) {
        size = cRTrackCol->size();
        tcolIter = cRTrackCol->begin();
        for(; tcolIter!=cRTrackCol->end(); ++tcolIter,++trackCount) {
            Event::TkrTrack* track = *tcolIter;
            trackVec.push_back(track);
        }
    //std::cout << cRTrackCol->size() << " CR tracks, trackCount = " << trackCount<< std::endl;
    }

    //std::cout << "final size: " << trackVec.size() << std::endl;
    size = trackVec.size();

    return trackVec;
}
