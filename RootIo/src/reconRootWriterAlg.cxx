#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"

//#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "Event/TopLevel/Event.h"
#include "Event/TopLevel/EventModel.h"

#include "Event/Recon/TkrRecon/TkrClusterCol.h"
#include "Event/Recon/TkrRecon/TkrPatCand.h"
#include "Event/Recon/TkrRecon/TkrFitTrack.h"
#include "Event/Recon/TkrRecon/TkrKalFitTrack.h"
#include "Event/Recon/TkrRecon/TkrVertex.h"

#include "Event/Recon/CalRecon/CalCluster.h"   
#include "Event/Recon/CalRecon/CalXtalRecData.h"   

#include "Event/Recon/AcdRecon/AcdRecon.h"

#include "LdfEvent/EventSummaryData.h"

#include "idents/CalXtalId.h"

#include "facilities/Util.h"
#include "commonData.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TVector3.h"

#include "reconRootData/ReconEvent.h"
#include "RootIo/IRootIoSvc.h"

// ADDED FOR THE FILE HEADERS DEMO
#include "RootIo/FhTool.h"
#include <cstdlib>

/** @class reconRootWriterAlg
* @brief Writes Recon TDS data to a persistent ROOT file.
*
* @author Heather Kelly and Tracy Usher
* $Header$
*/

class reconRootWriterAlg : public Algorithm
{	
public:
    
    reconRootWriterAlg(const std::string& name, ISvcLocator* pSvcLocator);
    
    /// Handles setup by opening ROOT file in write mode and creating a new TTree
    StatusCode initialize();
    
    /// Orchastrates reading from TDS and writing to ROOT for each event
    StatusCode execute();
    
    /// Closes the ROOT file and cleans up
    StatusCode finalize();
    
private:
    
    /// Retrieves event Id and run Id from TDS and fills the Recon ROOT object
    StatusCode writeReconEvent();
    
    /// Retrieves the TKR reconstruction data from the TDS and fills the TkrRecon
    /// ROOT object
    StatusCode writeTkrRecon();
    
    /// These are the methods specific to filling the pieces of the TkrRecon stuff
    void fillTkrClusterCol(TkrRecon* recon, Event::TkrClusterCol* clusterColTds);
    void       fillCandidateTracks(TkrRecon* recon, Event::TkrPatCandCol*  candidatesTds);
    void       fillFitTracks(      TkrRecon* recon, Event::TkrFitTrackCol* tracksTds);
    void       fillVertices(       TkrRecon* recon, Event::TkrVertexCol*   verticesTds, Event::TkrFitTrackCol* tracksTds);
    
    TkrHitPlane convertTkrHitPlane(Event::TkrFitPlane& planeTds);
    TObject*    convertTkrFitTrack(Event::TkrFitTrack* trackTds, int trkId);
    TObject*    convertTkrKalFitTrack(Event::TkrKalFitTrack* trackTds, int trkId);
    
    /// Retrieves the CAL reconstruction data from the TDS and fills the CalRecon
    /// ROOT object
    StatusCode writeCalRecon();
    
    /// Retrieves the CalCluster collection from the TDS and fills the ROOT    
    /// collection   
    void fillCalCluster(CalRecon *calRec, Event::CalClusterCol* clusterColTds);   
    
    /// Retrieves the CalXtalRecData collection from the TDS and fills the ROOT   
    /// xtal collection   
    void fillCalXtalRec(CalRecon *calRec, Event::CalXtalRecCol* xtalColTds); 
    
    StatusCode writeAcdRecon();
    
    /// Calls TTree::Fill for each event and clears m_mcEvt
    void writeEvent();
    
    /// Performs the final write to the ROOT file and closes
    void close();
    
    /// ROOT file pointer
    TFile *m_reconFile;
    /// ROOT tree pointer
    TTree *m_reconTree;
    /// Top-level Monte Carlo ROOT object
    ReconEvent *m_reconEvt;
    /// name of the output ROOT file
    std::string m_fileName;
    /// name of the TTree in the ROOT file
    std::string m_treeName;
    /// ROOT split mode
    int m_splitMode;
    /// Buffer Size for the ROOT file
    int m_bufSize;
    /// Compression level for the ROOT file
    int m_compressionLevel;

    /// Keep track of relation between TDS objs and ROOT counterparts
    commonData m_common;
    IRootIoSvc* m_rootIoSvc;

    // ADDED FOR THE FILE HEADERS DEMO
    IFhTool * m_headersTool ;
};


static const AlgFactory<reconRootWriterAlg>  Factory;
const IAlgFactory& reconRootWriterAlgFactory = Factory;

reconRootWriterAlg::reconRootWriterAlg(const std::string& name, 
                                       ISvcLocator* pSvcLocator) : 
Algorithm(name, pSvcLocator)
{
    // Input parameters available to be set via the jobOptions file
    declareProperty("reconRootFile",m_fileName="recon.root");
    declareProperty("splitMode", m_splitMode=1);
    declareProperty("bufferSize", m_bufSize=64000);
    // ROOT default compression
    declareProperty("compressionLevel", m_compressionLevel=1);
    declareProperty("treeName", m_treeName="Recon");    
    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool = 0 ;
}

StatusCode reconRootWriterAlg::initialize()
{
    // Purpose and Method:  Called once before the run begins.  This method
    //    opens a new ROOT file and prepares for writing.
    
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    StatusCode headersSc = toolSvc()->retrieveTool("FhTool",m_headersTool) ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to retreive headers tool" << endreq;
    }
    headersSc = m_headersTool->newReconHeader() ;
    if (headersSc.isFailure()) {
        log<<MSG::WARNING << "Failed to create a new Recon FileHeader" << endreq;
    }
    
    // Use the Job options service to set the Algorithm's parameters
    // This will retrieve parameters set in the job options file
    setProperties();
    
    m_rootIoSvc = 0 ;
    if ( service("RootIoSvc", m_rootIoSvc, true).isFailure() ){
        log << MSG::INFO << "Couldn't find the RootIoSvc!" << endreq;
        log << MSG::INFO << "No Auto Saving" << endreq;
        m_rootIoSvc = 0;
    } 

    facilities::Util::expandEnvVar(&m_fileName);
    
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;   
    // Create the new ROOT file
    m_reconFile = new TFile(m_fileName.c_str(), "RECREATE");
    if (!m_reconFile->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    m_reconFile->cd();
    m_reconFile->SetCompressionLevel(m_compressionLevel);
    m_reconTree = new TTree(m_treeName.c_str(), "GLAST Reconstruction Data");
    m_reconEvt = new ReconEvent();
    m_reconTree->Branch("ReconEvent","ReconEvent", &m_reconEvt, m_bufSize, m_splitMode);
    m_common.m_reconEvt = m_reconEvt;
    
    saveDir->cd();
    return sc;
    
}

StatusCode reconRootWriterAlg::execute()
{
    // Purpose and Method:  Called once per event.  This method calls
    //   the appropriate methods to read data from the TDS and write data
    //   to the ROOT file.
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    if(!m_reconTree->GetCurrentFile()->IsOpen()) {
        log << MSG::ERROR << "ROOT file " << m_fileName 
            << " could not be opened for writing." << endreq;
        return StatusCode::FAILURE;
    }
    
    if (m_reconEvt) m_reconEvt->Clear();

    sc = writeReconEvent();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write ReconEvent" << endreq;
        return sc;
    }
    
    sc = writeTkrRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Tkr Recon data" << endreq;
        return sc;
    }
    
    sc = writeCalRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Cal Recon Data" << endreq;
        return sc;
    }
    
    sc = writeAcdRecon();
    if (sc.isFailure()) {
        log << MSG::ERROR << "Failed to write Acd Recon Data" << endreq;
        return sc;
    }
    
    writeEvent();
    return sc;
}


StatusCode reconRootWriterAlg::writeReconEvent() {
    // Purpose and Method:  Retrieve the Event object from the TDS and set the
    //    event and run numbers in the ReconEvent ROOT object
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Retrieve the Event data for this event
    SmartDataPtr<Event::EventHeader> evtTds(eventSvc(), EventModel::EventHeader);
    
    if (!evtTds) return sc;
    
    UInt_t evtId = evtTds->event();
    UInt_t runId = evtTds->run();
    
    log << MSG::DEBUG;
    if( log.isActive())evtTds->fillStream(log.stream());
    log << endreq;
    
    m_reconEvt->initialize(evtId, runId, new TkrRecon, new CalRecon, new AcdRecon);

    // For simulated data - this may not exist on the TDS and that is ok
    // no need to fail for that
    SmartDataPtr<LdfEvent::EventSummaryData> summaryTds(eventSvc(), "/Event/EventSummary");
    if (summaryTds) 
        m_reconEvt->initEventFlags(summaryTds->eventFlags());
        
    
    return sc;
}

StatusCode reconRootWriterAlg::writeTkrRecon() {
    // Purpose and Method:  Retrieve the Tkr Recon data from the TDS
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Get pointer to TkrRecon part of ReconEvent
    TkrRecon* recon = m_reconEvt->getTkrRecon();
    if (!recon) return StatusCode::FAILURE;
    recon->initialize();
    
    SmartDataPtr<Event::TkrClusterCol> clusterColTds(eventSvc(), EventModel::TkrRecon::TkrClusterCol);
    if (clusterColTds) fillTkrClusterCol(recon, clusterColTds);
    
    // Retrieve the information on candidate tracks
    SmartDataPtr<Event::TkrPatCandCol> candidatesTds(eventSvc(), EventModel::TkrRecon::TkrPatCandCol);
    
    // Fill the candidate tracks
    if (candidatesTds) fillCandidateTracks(recon, candidatesTds);
    
    // Retrieve the information on fit tracks
    SmartDataPtr<Event::TkrFitTrackCol>    tracksTds(eventSvc(), EventModel::TkrRecon::TkrFitTrackCol);
    
    // Fill the fit tracks
    if (tracksTds) fillFitTracks(recon, tracksTds);
    
    // Retrieve the information on vertices
    SmartDataPtr<Event::TkrVertexCol> verticesTds(eventSvc(), EventModel::TkrRecon::TkrVertexCol);
    
    // Fill the vertices
    if (verticesTds && tracksTds) fillVertices(recon, verticesTds, tracksTds);
    
    return sc;
}

void reconRootWriterAlg::fillTkrClusterCol(TkrRecon* recon, Event::TkrClusterCol* clusterColTds) {
    
    // Purpose and Method:  Reads collection of clusters from TDS and copies their
    //   contents to a ROOT version.
    
    int numClusters = clusterColTds->nHits();
    
    int icluster;
    for (icluster = 0; icluster < numClusters; icluster++) {
        Event::TkrCluster *clusterTds = clusterColTds->getHit(icluster);
        Event::TkrCluster::view viewTds = clusterTds->v();
        TkrCluster::view viewRoot;
        if (viewTds == Event::TkrCluster::X) viewRoot = TkrCluster::X;
        else viewRoot = TkrCluster::Y;
        Point posTds = clusterTds->position();
        TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());
       
        TkrCluster *clusterRoot = new TkrCluster(clusterTds->id(), 
            clusterTds->plane(), viewRoot, clusterTds->firstStrip(), 
            clusterTds->lastStrip(), posRoot, clusterTds->ToT(), 
            clusterTds->hitFlagged(), clusterTds->tower());
        
        recon->addCluster(clusterRoot);
    }
}

void reconRootWriterAlg::fillCandidateTracks(TkrRecon* recon, Event::TkrPatCandCol* candidatesTds)
{
    // Purpose and Method:  This creates root candidate tracks from tds candidate tracks 
    //                      and adds them to the list kept in TkrRecon
    
    // Loop over the candidate tracks in the TDS
    unsigned int                     candId  = 0;
    Event::TkrPatCandColPtr cands   = candidatesTds->begin();
    while(candId < candidatesTds->size())
    {
        Event::TkrPatCand* candTds = *cands++;
        Double_t           x = candTds->getPosition().x();
        Double_t           y = candTds->getPosition().y();
        Double_t           z = candTds->getPosition().z();
        TVector3           pos(x,y,z);
        TVector3           dir(candTds->getDirection().x(),candTds->getDirection().z(),candTds->getDirection().z());
        TkrCandTrack *cand = new TkrCandTrack(candId, candTds->getLayer(), candTds->getTower(),
            candTds->getQuality(), candTds->getEnergy(), pos, dir);

        // Keep relation between Event and Root candidate tracks
        TRef ref = cand;
        m_common.m_tkrCandMap[candTds] = ref;

        // Now we go through the candidate hits (if they are included) and fill those objects
        int                nHits   = 0;
        while(nHits < candTds->numPatCandHits())
        {
            Event::TkrPatCandHit* candHitTds = candTds->getCandHit(nHits);
            TVector3              pos(candHitTds->Position().x(),candHitTds->Position().y(),candHitTds->Position().z());
            TkrCandHit::AXIS      view = candHitTds->View() == Event::TkrCluster::X ? TkrCandHit::X : TkrCandHit::Y;
            TkrCandHit candHit;
            
            candHit.initialize(pos, candHitTds->HitIndex(), candHitTds->TowerIndex(), candHitTds->PlaneIndex(), view); 
            
            cand->addHit(candHit);
            nHits++;
        }
        
        // Add the candidate to the list
        recon->addTrackCand(cand);
        candId++;
    }
}

void reconRootWriterAlg::fillFitTracks(TkrRecon* recon, Event::TkrFitTrackCol* tracksTds)
{
    // Purpose and Method:  This creates root tracks from tds tracks 
    //                      and adds them to the list kept in TkrRecon

    // Iterate over the tracks in the TDS
    int                 trkId  = 0;
    Event::TkrFitColPtr trkPtr = tracksTds->begin();
    while(trkPtr != tracksTds->end())
    {
        TObject*                trackRoot   = 0;
        Event::TkrFitTrackBase* trackBase   = *trkPtr++;       // The TDS track
        
        // Use dynamic casting to determine what type of fit track is in the TDS
        if (Event::TkrFitTrack* trackFitTds = dynamic_cast<Event::TkrFitTrack*>(trackBase))
            trackRoot = convertTkrFitTrack(trackFitTds, trkId++);
        
        else if (Event::TkrKalFitTrack* trackFitTds = dynamic_cast<Event::TkrKalFitTrack*>(trackBase))
            trackRoot = convertTkrKalFitTrack(trackFitTds, trkId++);

        // Keep relation between Event and Root candidate tracks
        TRef ref = trackRoot;
        m_common.m_tkrTrackMap[trackBase] = ref;

        // Ok, now add the track to the list!
        recon->addTrack(trackRoot);
    }
    
    return;
}

void reconRootWriterAlg::fillVertices(TkrRecon* recon, Event::TkrVertexCol* verticesTds, Event::TkrFitTrackCol* tracksTds)
{
    // Purpose and Method:  This creates root vertex output objects from tds vertices 
    //                      and adds them to the list kept in TkrRecon
    
    // Loop over the candidate tracks in the TDS
    int                     vtxId  = 0;
    Event::TkrVertexConPtr  vtxPtr = verticesTds->begin();
    while(vtxPtr != verticesTds->end())
    {
        Event::TkrVertex*   vtxTds = *vtxPtr++;
        TkrVertex*          vtx    = new TkrVertex();
        TVector3            pos(vtxTds->getPosition().x(), vtxTds->getPosition().y(), vtxTds->getPosition().z());
        TVector3            dir(vtxTds->getDirection().x(),vtxTds->getDirection().y(),vtxTds->getDirection().z());
        Event::TkrFitPar    parTds = vtxTds->getTrackPar();
        Event::TkrFitMatrix covTds = vtxTds->getTrackCov();
        
        // Keep relation between Event and Root vertices
        TRef ref = vtx;
        m_common.m_tkrVertexMap[vtxTds] = ref;

        vtx->initializeVals(TkrParams(parTds.getXPosition(), parTds.getXSlope(),
            parTds.getYPosition(), parTds.getYSlope()),
            TkrCovMat(covTds.getcovX0X0(), covTds.getcovX0Sx(), covTds.getcovX0Y0(), covTds.getcovX0Sy(),
            covTds.getcovSxX0(), covTds.getcovSxSx(), covTds.getcovSxY0(), covTds.getcovSxSy(),
            covTds.getcovY0X0(), covTds.getcovY0Sx(), covTds.getcovY0Y0(), covTds.getcovY0Sy(),
            covTds.getcovSyX0(), covTds.getcovSySx(), covTds.getcovSyY0(), covTds.getcovSySy()),
            pos,
            dir );
        vtx->initializeInfo(vtxId, 
            vtxTds->getLayer(), 
            vtxTds->getTower(), 
            vtxTds->getQuality(), 
            vtxTds->getEnergy());
        
        // Now add the track ids 
        // This is pretty ugly because we don't store track ids in the TDS classes
        SmartRefVector<Event::TkrFitTrackBase>::const_iterator vtxTrkIter = vtxTds->getTrackIterBegin();
        while(vtxTrkIter != vtxTds->getTrackIterEnd())
        {
            // Basically... take the pointer to the track from the vertex list and loop
            // over the track pointers in the track list looking for the match. Assign 
            // the id according to the loop variable value.
            // This ALWAYS succeeds (and I'm also selling valuable swampland in Florida!)
            int                 trkId  = 0;
            Event::TkrFitColPtr trkPtr = tracksTds->begin();
            SmartRef<Event::TkrFitTrackBase> vtxTrk = *vtxTrkIter++;
            
            while(trkPtr != tracksTds->end())
            {
                SmartRef<Event::TkrFitTrackBase> fitTrk = *trkPtr++;
                
                if (vtxTrk == fitTrk) break;
                
                trkId++;
            }
            
            vtx->addTrackId(trkId);
        }
        
        // Add the candidate to the list
        recon->addVertex(vtx);
        vtxId++;
    }
}

TObject* reconRootWriterAlg::convertTkrFitTrack(Event::TkrFitTrack* trackTds, int trkId)
{
    // Purpose and Method:  This converts Event::TkrFitTrack objects into root TkrTrack objects
    
    TkrTrack* track = new TkrTrack();      // Create a new Root Track
    
    // Initialize the track
    track->initializeInfo(trkId,
        trackTds->getNumXGaps(),
        trackTds->getNumYGaps(),
        trackTds->getNumXFirstGaps(),
        trackTds->getNumYFirstGaps() );
    
    track->initializeQual(trackTds->getChiSquare(),
        trackTds->getChiSquareSmooth(),
        trackTds->getScatter(),
        trackTds->getQuality(),
        trackTds->getKalEnergy(),
        trackTds->getKalThetaMS() );
    
    // Now loop over the hit planes and fill that information
    Event::TkrFitPlaneConPtr planePtr = trackTds->begin();
    while(planePtr != trackTds->end())
    {
        Event::TkrFitPlane  planeTds = *planePtr++;
        TkrHitPlane         plane    = convertTkrHitPlane(planeTds);
        
        track->addHit(plane);
    }
    
    return track;
}

TObject* reconRootWriterAlg::convertTkrKalFitTrack(Event::TkrKalFitTrack* trackTds, int trkId)
{
    // Purpose and Method:  This converts Event::TkrKalFitTrack objects into root TkrKalFitTrack objects
    
    TkrKalFitTrack* track    = new TkrKalFitTrack();
    
    // Need to extract status, initial position and direction first
    TkrKalFitTrack::Status status   = TkrKalFitTrack::EMPTY;
    TVector3               iniPos   = TVector3(trackTds->getInitialPosition().x(),
        trackTds->getInitialPosition().y(),
        trackTds->getInitialPosition().z());
    TVector3               iniDir   = TVector3(trackTds->getInitialDirection().x(),
        trackTds->getInitialDirection().y(),
        trackTds->getInitialDirection().z());
    
    if      (trackTds->status() == Event::TkrKalFitTrack::FOUND) status = TkrKalFitTrack::FOUND;
    else if (trackTds->status() == Event::TkrKalFitTrack::CRACK) status = TkrKalFitTrack::CRACK;
    
    // Initialize the track
    track->initializeBase(status, 
        trackTds->getType(),
        trkId,
        trackTds->getStartEnergy(),
        iniPos,
        iniDir);
    
    track->initializeQual(trackTds->getChiSquare(),
        trackTds->getChiSquareSmooth(),
        trackTds->getScatter(),
        trackTds->getQuality(),
        trackTds->getKalEnergy(),
        trackTds->getKalEnergyError(),
        trackTds->getKalThetaMS() );
    
    track->initializeGaps(trackTds->getNumXGaps(),
        trackTds->getNumYGaps(),
        trackTds->getNumXFirstGaps(),
        trackTds->getNumYFirstGaps());
    
    track->initializeKal( trackTds->getNumSegmentPoints(),
        trackTds->chiSquareSegment(),
        trackTds->getNumXHits(),
        trackTds->getNumYHits(),
        trackTds->getTkrCalRadlen());
    
    // Now loop over the hit planes and fill that information
    Event::TkrFitPlaneConPtr planePtr = trackTds->begin();
    while(planePtr != trackTds->end())
    {
        Event::TkrFitPlane  planeTds = *planePtr++;
        TkrHitPlane         plane    = convertTkrHitPlane(planeTds);
        
        track->addHit(plane);
    }
    
    return track;
}

TkrHitPlane reconRootWriterAlg::convertTkrHitPlane(Event::TkrFitPlane& planeTds)
{
    TkrHitPlane plane;
    
    // Wouldn't it be great if there was a STANDARD definition for these!
    TkrHitPlane::AXIS   proj     = planeTds.getProjection() == Event::TkrCluster::X ? TkrHitPlane::X : TkrHitPlane::Y;
    TkrHitPlane::AXIS   projPlus = planeTds.getNextProj()   == Event::TkrCluster::X ? TkrHitPlane::X : TkrHitPlane::Y;
    
    plane.initializeInfo(planeTds.getIDHit(),
        planeTds.getIDTower(),
        planeTds.getIDPlane(),
        proj,
        projPlus,
        planeTds.getZPlane(),
        planeTds.getEnergy(),
        planeTds.getRadLen(),
        planeTds.getActiveDist() );
    
    // Here we build the hit info (one at a time) starting with measured
    TkrParams           params;
    TkrCovMat           covMat;
    Event::TkrFitHit    hitTds = planeTds.getHit(Event::TkrFitHit::MEAS);
    Event::TkrFitPar    parTds = hitTds.getPar();
    Event::TkrFitMatrix covTds = hitTds.getCov();
    params.initialize(parTds.getXPosition(),
        parTds.getXSlope(),
        parTds.getYPosition(),
        parTds.getYSlope() );
    covMat.initialize(covTds.getcovX0X0(), covTds.getcovX0Sx(), covTds.getcovX0Y0(), covTds.getcovX0Sy(),
        covTds.getcovSxX0(), covTds.getcovSxSx(), covTds.getcovSxY0(), covTds.getcovSxSy(),
        covTds.getcovY0X0(), covTds.getcovY0Sx(), covTds.getcovY0Y0(), covTds.getcovY0Sy(),
        covTds.getcovSyX0(), covTds.getcovSySx(), covTds.getcovSyY0(), covTds.getcovSySy());
    
    TkrFitHit           measHit(TkrFitHit::MEAS, params, covMat);
    
    // Now the predicted hit
    hitTds = planeTds.getHit(Event::TkrFitHit::PRED);
    parTds = hitTds.getPar();
    covTds = hitTds.getCov();
    params.initialize(parTds.getXPosition(),
        parTds.getXSlope(),
        parTds.getYPosition(),
        parTds.getYSlope() );
    covMat.initialize(covTds.getcovX0X0(), covTds.getcovX0Sx(), covTds.getcovX0Y0(), covTds.getcovX0Sy(),
        covTds.getcovSxX0(), covTds.getcovSxSx(), covTds.getcovSxY0(), covTds.getcovSxSy(),
        covTds.getcovY0X0(), covTds.getcovY0Sx(), covTds.getcovY0Y0(), covTds.getcovY0Sy(),
        covTds.getcovSyX0(), covTds.getcovSySx(), covTds.getcovSyY0(), covTds.getcovSySy());
    
    TkrFitHit           predHit(TkrFitHit::PRED, params, covMat);
    
    // Now the fit (filtered) hit
    hitTds = planeTds.getHit(Event::TkrFitHit::FIT);
    parTds = hitTds.getPar();
    covTds = hitTds.getCov();
    params.initialize(parTds.getXPosition(),
        parTds.getXSlope(),
        parTds.getYPosition(),
        parTds.getYSlope() );
    covMat.initialize(covTds.getcovX0X0(), covTds.getcovX0Sx(), covTds.getcovX0Y0(), covTds.getcovX0Sy(),
        covTds.getcovSxX0(), covTds.getcovSxSx(), covTds.getcovSxY0(), covTds.getcovSxSy(),
        covTds.getcovY0X0(), covTds.getcovY0Sx(), covTds.getcovY0Y0(), covTds.getcovY0Sy(),
        covTds.getcovSyX0(), covTds.getcovSySx(), covTds.getcovSyY0(), covTds.getcovSySy());
    
    TkrFitHit           filtHit(TkrFitHit::FIT, params, covMat);
    
    // Now the smoothed hit
    hitTds = planeTds.getHit(Event::TkrFitHit::SMOOTH);
    parTds = hitTds.getPar();
    covTds = hitTds.getCov();
    params.initialize(parTds.getXPosition(),
        parTds.getXSlope(),
        parTds.getYPosition(),
        parTds.getYSlope() );
    covMat.initialize(covTds.getcovX0X0(), covTds.getcovX0Sx(), covTds.getcovX0Y0(), covTds.getcovX0Sy(),
        covTds.getcovSxX0(), covTds.getcovSxSx(), covTds.getcovSxY0(), covTds.getcovSxSy(),
        covTds.getcovY0X0(), covTds.getcovY0Sx(), covTds.getcovY0Y0(), covTds.getcovY0Sy(),
        covTds.getcovSyX0(), covTds.getcovSySx(), covTds.getcovSyY0(), covTds.getcovSySy());
    
    TkrFitHit           smooHit(TkrFitHit::SMOOTH, params, covMat);
    
    // Here retrieve the scattering matrix
    Event::TkrFitMatrix scatTds = planeTds.getQmaterial();
    TkrCovMat         scatCov(scatTds.getcovX0X0(), scatTds.getcovX0Sx(), scatTds.getcovX0Y0(), scatTds.getcovX0Sy(),
        scatTds.getcovSxX0(), scatTds.getcovSxSx(), scatTds.getcovSxY0(), scatTds.getcovSxSy(),
        scatTds.getcovY0X0(), scatTds.getcovY0Sx(), scatTds.getcovY0Y0(), scatTds.getcovY0Sy(),
        scatTds.getcovSyX0(), scatTds.getcovSySx(), scatTds.getcovSyY0(), scatTds.getcovSySy());
    
    plane.initializeHits(measHit, predHit, filtHit, smooHit, scatCov);
    
    return plane;
}


StatusCode reconRootWriterAlg::writeCalRecon() {
    // Purpose and Method:  Retrieve the Cal Recon data from the TDS 
    //  calls the two helper methods that do the work of filling the ROOT    
    //  version of the CalRecon object for this event. 
    
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Get pointer to CalRecon part of ReconEvent   
    CalRecon* calRec = m_reconEvt->getCalRecon();   
    if (!calRec) return StatusCode::FAILURE;   
    calRec->initialize();   
    
    // Retrieve the cal cluster collection   
    SmartDataPtr<Event::CalClusterCol> clusterColTds(eventSvc(), EventModel::CalRecon::CalClusterCol);   
    if (clusterColTds) fillCalCluster(calRec, clusterColTds);   
    
    // Retrieve the cal xtal collection   
    SmartDataPtr<Event::CalXtalRecCol> xtalColTds(eventSvc(), EventModel::CalRecon::CalXtalRecCol);   
    if (xtalColTds) fillCalXtalRec(calRec, xtalColTds); 
    
    
    return sc;
}

void reconRootWriterAlg::fillCalCluster(CalRecon *calRec, Event::CalClusterCol* clusterColTds) {   
    // Purpose and Method:  Given the CalCluster collection from the TDS, we fill the ROOT   
    //  CalCluster collection.   
    
    unsigned int numClusters = clusterColTds->num();   
    unsigned int iCluster;   
    for (iCluster = 0; iCluster < numClusters; iCluster++) {   
        Event::CalCluster *clusterTds = clusterColTds->getCluster(iCluster);   
        Point posTds = clusterTds->getPosition();   
        TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());   
        
        CalCluster *clusterRoot = new CalCluster(   
            clusterTds->getEnergySum(), posRoot);   
        
        Vector dirTds = clusterTds->getDirection();   
        TVector3 dirRoot(dirTds.x(), dirTds.y(), dirTds.z());   
        
        std::vector<Vector> posLayerTds = clusterTds->getPosLayer();   
        std::vector<Vector> rmsLayerTds = clusterTds->getRmsLayer();   
        std::vector<TVector3> posLayerRoot;   
        std::vector<Vector>::const_iterator tdsIt;   
        for (tdsIt = posLayerTds.begin(); tdsIt != posLayerTds.end(); tdsIt++) {   
            TVector3 curVec(tdsIt->x(), tdsIt->y(), tdsIt->z());   
            posLayerRoot.push_back(curVec);   
        }   
        std::vector<TVector3> rmsLayerRoot;   
        for (tdsIt = rmsLayerTds.begin(); tdsIt != rmsLayerTds.end(); tdsIt++) {   
            TVector3 curVec(tdsIt->x(), tdsIt->y(), tdsIt->z());   
            rmsLayerRoot.push_back(curVec);   
        }   
        
        clusterRoot->initialize(clusterTds->getEnergyLeak(), 
            clusterTds->getEnergyCorrected(),
            clusterTds->getEneLayer(),   
            posLayerRoot, rmsLayerRoot,    
            clusterTds->getRmsLong(), clusterTds->getRmsTrans(),    
            dirRoot, clusterTds->getTransvOffset());   
        
        clusterRoot->initProfile(clusterTds->getFitEnergy(),    
            clusterTds->getProfChisq(), clusterTds->getCsiStart(),   
            clusterTds->getCsiAlpha(), clusterTds->getCsiLambda());   
        
        calRec->addCalCluster(clusterRoot);   
    }   
    
    return;   
}   

void reconRootWriterAlg::fillCalXtalRec(CalRecon *calRec, Event::CalXtalRecCol* xtalColTds) {   
    // Purpose and Method:  Given the CalXtalRecData collection from the TDS,   
    //   this method fills the ROOT CalXtalRecData collection.   
    
    Event::CalXtalRecCol::const_iterator xtalTds;   
    
    for (xtalTds = xtalColTds->begin(); xtalTds != xtalColTds->end(); xtalTds++) {   
        CalXtalRecData *xtalRoot = new CalXtalRecData();   
        idents::CalXtalId::CalTrigMode modeTds = (*xtalTds)->getMode();   
        idents::CalXtalId idTds = (*xtalTds)->getPackedId();   
        CalXtalId idRoot;   
        idRoot.init(idTds.getTower(), idTds.getLayer(), idTds.getColumn());   
        if (modeTds == idents::CalXtalId::BESTRANGE) {   
            xtalRoot->initialize(CalXtalId::BESTRANGE, idRoot);   
            Event::CalXtalRecData::CalRangeRecData *xtalRangeTds =    
                (*xtalTds)->getRangeRecData(0);   
            CalRangeRecData recRoot(   
                xtalRangeTds->getRange(idents::CalXtalId::POS),    
                xtalRangeTds->getEnergy(idents::CalXtalId::POS),    
                xtalRangeTds->getRange(idents::CalXtalId::NEG),   
                xtalRangeTds->getEnergy(idents::CalXtalId::NEG) );   
            Point posTds = xtalRangeTds->getPosition();   
            TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());   
            recRoot.initialize(posRoot);   
            xtalRoot->addRangeRecData(recRoot);   
        } else {   
            xtalRoot->initialize(CalXtalId::ALLRANGE, idRoot);   
            int range;   
            for (range = idents::CalXtalId::LEX8; range <= idents::CalXtalId::HEX1; range++) {   
                Event::CalXtalRecData::CalRangeRecData *xtalRangeTds =    
                    (*xtalTds)->getRangeRecData(range);   
                CalRangeRecData recRoot(   
                    xtalRangeTds->getRange(idents::CalXtalId::POS),    
                    xtalRangeTds->getEnergy(idents::CalXtalId::POS),    
                    xtalRangeTds->getRange(idents::CalXtalId::NEG),   
                    xtalRangeTds->getEnergy(idents::CalXtalId::NEG) );   
                Point posTds = xtalRangeTds->getPosition();   
                TVector3 posRoot(posTds.x(), posTds.y(), posTds.z());   
                recRoot.initialize(posRoot);   
                xtalRoot->addRangeRecData(recRoot);   
            }   
        }   
        
        calRec->addXtalRecData(xtalRoot);   
    }   
    
    return;   
} 


StatusCode reconRootWriterAlg::writeAcdRecon() {
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    
    // Get pointer to AcdRecon part of ReconEvent   
    AcdRecon* acdRec = m_reconEvt->getAcdRecon();   
    if (!acdRec) return StatusCode::FAILURE;
    SmartDataPtr<Event::AcdRecon> acdRecTds(eventSvc(), EventModel::AcdRecon::Event);  
    if (!acdRecTds) return StatusCode::SUCCESS;
    idents::AcdId acdIdTds = acdRecTds->getMinDocaId();
    std::vector<AcdId> idRootCol;
    std::vector<idents::AcdId>::const_iterator idTdsIt;
    for (idTdsIt = acdRecTds->getIdCol().begin(); idTdsIt != acdRecTds->getIdCol().end(); idTdsIt++) {
        idRootCol.push_back(AcdId(idTdsIt->layer(), idTdsIt->face(), 
            idTdsIt->row(), idTdsIt->column()));
    }
    AcdId acdIdRoot(acdIdTds.layer(), acdIdTds.face(), acdIdTds.row(), acdIdTds.column());
    acdRec->initialize(acdRecTds->getEnergy(), acdRecTds->getTileCount(),
        acdRecTds->getGammaDoca(), acdRecTds->getDoca(), acdRecTds->getActiveDist(), acdIdRoot, 
        acdRecTds->getRowDocaCol(), acdRecTds->getRowActDistCol(),
        idRootCol, acdRecTds->getEnergyCol());
    
    return sc;
}

void reconRootWriterAlg::writeEvent() 
{
    // Purpose and Method:  Stores the DigiEvent data for this event in the ROOT
    //    tree.  The m_digiEvt object is cleared for the next event.
    
    static int eventCounter = 0;
    TDirectory *saveDir = gDirectory;
    m_reconTree->GetCurrentFile()->cd();
    //m_reconFile->cd();
    m_reconTree->Fill();
    ++eventCounter;
    if (m_rootIoSvc)
        if (eventCounter % m_rootIoSvc->getAutoSaveInterval() == 0) m_reconTree->AutoSave();
    saveDir->cd();

    saveDir->cd();
    return;
}

void reconRootWriterAlg::close() 
{
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kWriteDelete parameter is used in the Write method
   //    Used rather than TObject::kOverwrite - supposed to be safer but slower
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite.
    
    TDirectory *saveDir = gDirectory;
    TFile *f = m_reconTree->GetCurrentFile();
    f->cd();
    //m_reconFile->cd();
    m_reconTree->BuildIndex("m_runId", "m_eventId");
    f->Write(0, TObject::kWriteDelete);
    f->Close();
    saveDir->cd();
    return;
}

StatusCode reconRootWriterAlg::finalize()
{
    MsgStream log(msgSvc(), name());
    
    // ADDED FOR THE FILE HEADERS DEMO
    m_headersTool->writeReconHeader(m_reconTree->GetCurrentFile()) ;
    
    close();
    
    StatusCode sc = StatusCode::SUCCESS;
    setFinalized();

    log << MSG::DEBUG << "Finalized" << endreq;

    return sc;
}

