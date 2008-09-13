/**
* @class TkrGhostTool
*
* @brief Implements a Gaudi Tool for setting the candidate track energies before 
*        the track fit
*
* @author The Tracking Software Group
*
* $Header$
*/

#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/TkrRecon/TkrTrack.h"
#include "Event/Recon/TkrRecon/TkrEventParams.h"
#include "Event/TopLevel/DigiEvent.h"
#include "Event/Digi/TkrDigi.h"

#include "LdfEvent/Gem.h"
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "src/Track/TkrGhostTool.h"
#include "src/Track/TkrControl.h"
#include "TkrUtil/ITkrGeometrySvc.h"

#include <iomanip>
#include <map>

class TkrGhostTool : public AlgTool, virtual public ITkrGhostTool
{
public:
    /// Standard Gaudi Tool interface constructor
    TkrGhostTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~TkrGhostTool() {}

    /// @brief Method to fit a single candidate track. Will retrieve any extra info 
    ///        needed from the TDS, then create and use a new KalFitTrack object to 
    ///        fit the track via a Kalman Filter. Successfully fit tracks are then 
    ///        added to the collection in the TDS.
    StatusCode initialize();

    StatusCode flagSingles();
    StatusCode flagEarlyHits();
    StatusCode flagEarlyTracks();

    /// @brief Tool for identifying and flagging ghost clusters

private:

    TkrControl*            m_control;
    Event::DigiEvent*      m_digiEvent;

    /// Pointer to the Gaudi data provider service
    DataSvc*               m_dataSvc;
    IGlastDetSvc*          m_pDetSvc;
    ITkrGeometrySvc*       m_tkrGeom;

    int m_numTowers;
    int m_numLayers;
};

static ToolFactory<TkrGhostTool> s_factory;
const IToolFactory& TkrGhostToolFactory = s_factory;

TkrGhostTool::TkrGhostTool(const std::string& type, 
                           const std::string& name, 
                           const IInterface* parent) :
AlgTool(type, name, parent)
{
    //Declare the additional interface
    declareInterface<ITkrGhostTool>(this);

    return;
}

StatusCode TkrGhostTool::initialize()
{
    // Purpose and Method: finds the "Global Event Energy" and constrains the
    //                     first two track energies to sum to it.
    // Inputs:  Calorimeter Energy
    // Outputs: Sets the "constrained" energy for all Candidates
    //          Note: This is the energy that will be used for the 
    //          final track fit. 
    // Dependencies: None
    // Restrictions and Caveats:  None.

    //Always believe in success
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());

    m_control = TkrControl::getPtr();

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    sc = service("GlastDetSvc", m_pDetSvc);
    if(sc.isFailure()) {
        log << MSG::ERROR << "GlastDetSvc not found!" << endreq;
        return sc;
    }
    sc = service("TkrGeometrySvc", m_tkrGeom);
    if(sc.isFailure()) {
        log << MSG::ERROR << "TkrGeometrySvc not found!" << endreq;
        return sc;
    }

    int numX, numY;
    m_pDetSvc->getNumericConstByName("xNum", &numX);
    m_pDetSvc->getNumericConstByName("yNum", &numY);  
    m_numTowers = numX*numY;

    m_numLayers = m_tkrGeom->numLayers();

    return sc;
}

StatusCode TkrGhostTool::flagSingles() 
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;


    std::map<int, int> hC[3];

    //get the clusters
    SmartDataPtr<Event::TkrClusterCol> 
        clusterCol(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);

    int clusSize = clusterCol->size();
    int i;
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        idents::TkrId tkrid = clus->getTkrId();
        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
        int tower = twrid.id();
        int layer = clus->getLayer();
        int view  = tkrid.getView();
        int end   = clus->getEnd();

        int index = 1000*tower + 2*layer + view ;

        if(hC[end].find(index)==hC[end].end()) hC[end][index] = 0;
        hC[end][index]++;
    }

    // Now mark lone hits

    for(i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        int tower = clus->tower();
        int layer = clus->getLayer();
        int end   = clus->getEnd();
        idents::TkrId id = clus->getTkrId();
        int view = id.getView();
        int index = 1000*tower + 2*layer + view;

        int totCount;
        int endCount[2] = {0,0};
        if(hC[0].find(index)!=hC[0].end()) {
            endCount[0] = hC[0][index];
        }
        if(hC[1].find(index)!=hC[1].end()) {
            endCount[1] = hC[1][index];
        }
        totCount = endCount[0] + endCount[1];

        if(end==2) { // this might never happen... 
            if(totCount==0) {
                clus->setMask(Event::TkrCluster::maskALONE);
                clus->setMask(Event::TkrCluster::maskALONEEND);
            }
        } else {
            if(totCount==1) clus->setMask(Event::TkrCluster::maskALONE);
            if(endCount[end]==1) clus->setMask(Event::TkrCluster::maskALONEEND);
        }
    }
    return sc;
}

StatusCode TkrGhostTool::flagEarlyHits()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    //get the tkrVector
    // Retrieve the Event Summary data for this event
    SmartDataPtr<LdfEvent::Gem> gem(m_dataSvc, "/Event/Gem");

    if (!gem) {
        log << MSG::DEBUG << "No GEM found on TDS" << endreq;
        return sc;
    }

    unsigned int tkrVector = gem->tkrVector();

    //get the clusters
    SmartDataPtr<Event::TkrClusterCol> 
        clusterCol(m_dataSvc, EventModel::TkrRecon::TkrClusterCol);

    //set up the tower bits vector
    towerVec clusterTrigger(m_numTowers);

    int i;
    for(i=0;i<m_numTowers;++i) {
        clusterTrigger[i] = new TkrTowerBits();
    }

    //fill the bits vector

    int clusSize = clusterCol->size();
    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        int tower = clus->tower();
        clusterTrigger[tower]->setBit(clus);
    }

    //generate the tower trigger word

    unsigned int trigBits = 0;
    std::vector<unsigned int> towerBits(m_numTowers,0);
    for(i=0;i<m_numTowers;++i) {
        towerBits[i] = clusterTrigger[i]->getTriggeredBits();
        if(towerBits[i]>0) {trigBits |= (1<<i);}
    }

    if((tkrVector&trigBits)!=tkrVector) {
        log << MSG::DEBUG;
        if (log.isActive()) {
            log <<"tkrVector and calculated tower trigger disagree:" << endreq 
                << "tkrVector = " << std::hex << tkrVector 
            << ", calculation = " << trigBits << std::dec << endreq;
        }
    }

    // Flag the hits that would have made the trigger
    // This is an assumption... for a three-in-a-row, any one
    //   missing hit would disable the trigger
    // Also, flag ToT==255

    for (i=0;i<clusSize;++i) {
        Event::TkrCluster* clus = (*clusterCol)[i];
        idents::TkrId tkrid = clus->getTkrId();
        idents::TowerId twrid = idents::TowerId(tkrid.getTowerX(), tkrid.getTowerY());
        int tower = twrid.id();
        // we only look at towers with software trigger but no hardware trigger
        if((tkrVector&(1<<tower))!=0) continue;
        if((trigBits&(1<<tower))==0) continue;

        int layer = clus->getLayer();
        if(clus->getRawToT()==255) clus->setMask(Event::TkrCluster::mask255);
        if(towerBits[tower]&(1<<layer)) {
            clus->setMask(Event::TkrCluster::maskGHOST);

            log << MSG::DEBUG << "Ghost bit set for cluster " 
                << i << ", t/l "  << tower << ", " << layer  << ", isGhost "  
                << clus->isSet(Event::TkrCluster::maskGHOST) <<endreq;
            log << "tower bits " << std::hex << towerBits[tower] 
            << " "  << (1<<layer) << std::dec << endreq;
        }
    }
    return sc;
}

StatusCode TkrGhostTool::flagEarlyTracks()
{

    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;

    // Flag the hits on tracks containing ghosts or 255's
    //get the tracks
    SmartDataPtr<Event::TkrTrackCol> 
        trackCol(m_dataSvc, EventModel::TkrRecon::TkrTrackCol);

    int trackCount = 0;
    Event::TkrTrackColConPtr tcolIter = trackCol->begin();
    for(tcolIter; tcolIter!=trackCol->end();++tcolIter,++trackCount) {
        Event::TkrTrack* track = *tcolIter;
        Event::TkrTrackHitVecItr pHit = track->begin();

        int ghostCount = 0;
        int _255Count = 0;
        int count = 0;
        while(pHit != track->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            Event::TkrClusterPtr pClus = hit->getClusterPtr();
            if(!pClus) continue;
            bool is255      = pClus->isSet(Event::TkrCluster::mask255);
            bool isGhost    = pClus->isSet(Event::TkrCluster::maskGHOST);
            bool isAlone    = pClus->isSet(Event::TkrCluster::maskALONE);
            bool isAloneEnd = pClus->isSet(Event::TkrCluster::maskALONEEND);
            if(is255&&isAloneEnd) _255Count++;
            if(isGhost)           ghostCount++;
            count++;
        }

        if(_255Count==0&&ghostCount==0) continue;

        count = 0;
        pHit = track->begin();
        while(pHit!=track->end()) {
            const Event::TkrTrackHit* hit = *pHit++;
            Event::TkrClusterPtr pClus = hit->getClusterPtr();
            if(!pClus) continue;
            if(_255Count>0||ghostCount>0) {
                bool is255   = pClus->isSet(Event::TkrCluster::mask255);
                bool isGhost = pClus->isSet(Event::TkrCluster::maskGHOST);
                bool isAlone = pClus->isSet(Event::TkrCluster::maskALONE);
                bool isAloneEnd = pClus->isSet(Event::TkrCluster::maskALONEEND);
                if(!is255&&!isGhost) {
                    pClus->setMask(Event::TkrCluster::maskSAMETRACK);
                }
            } 
            count++;
        }
    }

    return sc;
}
