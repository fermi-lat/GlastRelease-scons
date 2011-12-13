#include "IMipFindingTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/IToolSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalClusterTab.h"
#include "Event/Recon/CalRecon/CalMipClasses.h"

#include "GlastSvc/Reco/IPropagatorTool.h"
#include "GlastSvc/Reco/IPropagatorSvc.h"
#include "GlastSvc/Reco/IPropagator.h" 
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#include "TkrUtil/ITkrGeometrySvc.h"

#include "idents/TowerId.h" 
#include "idents/VolumeIdentifier.h"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"

#include "geometry/Vector.h"
#include "geometry/Ray.h"

#include "TMath.h"

/**   
* @class StdMipFindingTool
*
* $Header$
*/

//-----------------------------------------------------------------------------------------------------------------
class StdMipFindingTool : public IMipFindingTool,  public AlgTool 
{
public:    
    StdMipFindingTool(const std::string & type, const std::string & name, const IInterface * parent );
    virtual ~StdMipFindingTool() {};
  
    /// @brief Intialization of the tool
    virtual StatusCode initialize();
    /// @brief Default cluster finding framework
    virtual StatusCode findMIPCandidates();
  
private:
    /// Private methods
    StatusCode               readGlastDet();
    void                     getSingleClusterCentroid();
    int                      findMipXtals();
    int                      findMipTracks();
    void                     trackProperties();
    bool                     findC0();
    bool                     findC1();
    bool                     findC2();
    bool                     findCn();
    int                      propagate(Point x0, Vector dir);
    bool                     leastSquares();
    double                   D2Edge             (Event::CalMipXtal     *calMipXtal);
    void                     clearCalMipTrackVec(Event::CalMipTrackVec *calMipTrackVec);
    void                     clearCalMipXtalVec (Event::CalMipXtalVec  *calMipXtalVec);
    void                     readCalMipXtals    (Event::CalMipXtalVec   calMipXtalVec);
    void                     readCalMipTracks   (Event::CalMipTrackVec  calMipTrackVec);
    void                     readCalMipTrackCol();
    StatusCode               storeCalMipTracks();

    /// Private data members
    /// Cut values for MIP energy
    double             m_e1;
    double             m_e2;
    double             m_ecor1;
    double             m_ecor2;
    double             m_pi;
    double             m_radToDeg;

    /// MsgStream member variable to speed up execution
    MsgStream          m_log;

    /// Pointer to the Gaudi data provider service
    DataSvc*           m_dataSvc;

    IPropagator * m_G4PropTool; 

    /// the GlastDetSvc used for access to detector info
    IGlastDetSvc*      m_detSvc;

    /// TkrGeometrySvc used for access to tracker geometry info
    ITkrGeometrySvc*   m_geoSvc;
    //@@@FP 07/10/05

    double             m_CsILength;
    double             m_CsIWidth;
    double             m_CsIHeight;

    int                m_hitId[16][8][12];

    double             m_calZTop;
    double             m_calZBot;
    double             m_calXLo;
    double             m_calXHi;
    double             m_calYLo;
    double             m_calYHi;

    //  Mip Xtals candidate vector from recon hit list;
    Event::CalMipXtalVec     m_calMipXtalVec;
    // Mip track candidate vector from Mip Xtals candidate vector
    Event::CalMipTrackVec    m_calMipTrackVec;
    Event::CalMipTrackCol*   m_calMipTrackCol;

    Point              m_singleClusterCentroid;

    int                m_hid;
    int                m_hid0;
    int                m_hid1;
    int                m_hid2;
    int                m_hidn;

    Vector             m_uu;
    Vector             m_vv;
    Vector             m_dir;
    Point              m_refP;
} ;

//-----------------------------------------------------------------------------------------------------------------
//static ToolFactory<StdMipFindingTool> s_factory;
//const IToolFactory& StdMipFindingToolFactory = s_factory;
DECLARE_TOOL_FACTORY(StdMipFindingTool);

//-----------------------------------------------------------------------------------------------------------------
StdMipFindingTool::StdMipFindingTool(const std::string & type, 
                                     const std::string & nameIn,
                                     const IInterface * parent ) : AlgTool( type, nameIn, parent ),
                                                                   m_log(msgSvc(), name())
{ 
    declareInterface<IMipFindingTool>(this) ; 

    declareProperty("MIPEne1",       m_e1         =  7.0      );
    declareProperty("MIPEne2",       m_e2         = 30.0      );
    //50% of Landau area
    //     declareProperty("Ecor1",       m_ecor1         = 11.05      );
    //     declareProperty("Ecor2",       m_ecor2         = 13.05     );
    //20%
    //     declareProperty("Ecor1",       m_ecor1         = 10.50      );
    //     declareProperty("Ecor2",       m_ecor2         = 14.45     );
    //10%
    //declareProperty("Ecor1",       m_ecor1         = 10.25      );
    //declareProperty("Ecor2",       m_ecor2         = 15.50     );

    declareProperty("Ecor1",       m_ecor1         =  7.0      );
    declareProperty("Ecor2",       m_ecor2         = 30.0      );

    return;
}
    
//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::initialize()
{
    StatusCode sc = StatusCode::SUCCESS;

    m_log << MSG::INFO << "BEGIN initialize()" << endreq ;

    m_pi = TMath::Pi();
    m_radToDeg = 180. / m_pi;

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
      throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);
    
    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
      m_log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
      return StatusCode::FAILURE;
    }

    // find TkrGeometrySvc service
    if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
      m_log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
      return StatusCode::FAILURE;
    }

    // find G4 propagation tool
    if(!toolSvc()->retrieveTool("G4PropagationTool", m_G4PropTool)) {
      m_log << MSG::ERROR << "Couldn't find the G4PropationTool!" << endreq;
      return StatusCode::FAILURE;
    }

    // read geometry
    readGlastDet();

    return StatusCode::SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::readGlastDet()
{
    StatusCode sc = StatusCode::SUCCESS;
    //    m_log << MSG::DEBUG << "readGlastDet in StdMipFindingTool" << endreq;
    
    m_detSvc->getNumericConstByName("CsILength", &m_CsILength);
    m_detSvc->getNumericConstByName("CsIWidth",  &m_CsIWidth);
    m_detSvc->getNumericConstByName("CsIHeight", &m_CsIHeight);
 
    m_calZTop = m_geoSvc->calZTop();
    m_calZBot = m_geoSvc->calZBot();

    double towerPitch = m_geoSvc->towerPitch();
    int xNum = m_geoSvc->numXTowers();
    int yNum = m_geoSvc->numYTowers();
    double calXWidth = m_geoSvc->calXWidth();
    double calYWidth = m_geoSvc->calYWidth();
    double deltaX = 0.5*(xNum*towerPitch - calXWidth);
    double deltaY = 0.5*(yNum*towerPitch - calYWidth);

    m_calXLo = m_geoSvc->getLATLimit(0,LOW)  + deltaX;
    m_calXHi = m_geoSvc->getLATLimit(0,HIGH) - deltaX;
    m_calYLo = m_geoSvc->getLATLimit(1,LOW)  + deltaY;
    m_calYHi = m_geoSvc->getLATLimit(1,HIGH) - deltaY;

//     m_log << MSG::DEBUG << "readGlastDet -end-" << endreq;
    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::getSingleClusterCentroid()
{
    //    m_log << MSG::DEBUG << "findSingleClusterCentroid in StdMipFindingTool" << endreq;

    Event::CalClusterCol* pCals = SmartDataPtr<Event::CalClusterCol>(m_dataSvc, EventModel::CalRecon::CalClusterCol);

    if (pCals != 0 && !pCals->empty())
    {
        Event::CalCluster* calCluster = pCals->front();

        m_singleClusterCentroid=calCluster->getPosition();
        //    m_log << MSG::DEBUG << "findSingleClusterCentroid -end- x=" << m_singleClusterCentroid.x() << " y=" << m_singleClusterCentroid.y() << " z=" << m_singleClusterCentroid.z() << endreq;
    }
    else m_singleClusterCentroid = Point(0.,0.,0.);

    return;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipXtals()
{
    //    m_log << MSG::DEBUG << "findMipXtals in StdMipFindingTool" << endreq;

    double emip1= 2.;
    double emip2=50.;

    int numMipXtals=0;

    for (int itow=0; itow<16; itow++)
        for (int ilay=0; ilay<8; ilay++)
            for (int icol=0; icol<12; icol++)
                m_hitId[itow][ilay][icol]=-1;

    Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(m_dataSvc, EventModel::CalRecon::CalXtalRecCol); 

    // loop over CalXtalRecdata
    for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
    {
        Event::CalXtalRecData* xTalData = *xTalIter;
        double                 xTalE    =  xTalData->getEnergy();
        if (xTalE>emip1 && xTalE<emip2)
        {
            m_calMipXtalVec.push_back(Event::CalMipXtal());
        
            int itow=xTalData->getPackedId().getTower();
            int ilay=xTalData->getPackedId().getLayer();
            int icol=xTalData->getPackedId().getColumn();
            m_hitId[itow][ilay][icol]=numMipXtals;

            numMipXtals++;

            Event::CalMipXtal& mipXtal = m_calMipXtalVec.back();
            mipXtal.setFree(true);
            mipXtal.setFreeC0(true);
            mipXtal.setXtal(xTalData);
            m_uu = m_singleClusterCentroid - xTalData->getPosition();
            double d2C=sqrt(m_uu*m_uu);
            mipXtal.setD2C(d2C);
            mipXtal.setEcor(-1);
        }
    }

    //    m_log << MSG::DEBUG << "findMipXtals -end- numMipXtals=" << numMipXtals << endreq;
    return numMipXtals;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readCalMipXtals(Event::CalMipXtalVec calMipXtalVec)
{
    //    m_log << MSG::DEBUG << "readCalMipXtals in StdMipFindingTool" << endreq;
    int counter=0;

    //Loop over crystals in the calMipXtalVec
    for(Event::CalMipXtalVec::iterator xTalIter=calMipXtalVec.begin(); xTalIter != calMipXtalVec.end(); xTalIter++)
    {
        Event::CalMipXtal calMipXtal = *xTalIter;

//         m_log << MSG::DEBUG << "readCalMipXtals - MipXtal No " << counter << "-------" << endreq;
//         m_log << MSG::DEBUG << "readCalMipXtals - Free=" << calMipXtal.getFree()  << endreq;
//         m_log << MSG::DEBUG << "readCalMipXtals - D2C =" << calMipXtal.getD2C()   << endreq;
//         m_log << MSG::DEBUG << "readCalMipXtals - Ener=" << calMipXtal.getXtal()->getEnergy()         << endreq;
//         m_log << MSG::DEBUG << "readCalMipXtals - PosX=" << calMipXtal.getXtal()->getPosition().x()   << endreq;
//         m_log << MSG::DEBUG << "readCalMipXtals - PosY=" << calMipXtal.getXtal()->getPosition().y()   << endreq;
//         m_log << MSG::DEBUG << "readCalMipXtals - PosZ=" << calMipXtal.getXtal()->getPosition().z()   << endreq;
        
//calMipXtal.print();
        counter++;
    }
        
    //    m_log << MSG::DEBUG << "readCalMipXtals -end-" << endreq;
    return;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readCalMipTracks(Event::CalMipTrackVec calMipTrackVec)
{
    m_log << MSG::DEBUG << "readCalMipTracks in StdMipFindingTool" << endreq;
    int counter=0;
    //Loop over tracks in the calMipTrackVec
    for(Event::CalMipTrackVec::iterator trackIter=calMipTrackVec.begin(); trackIter != calMipTrackVec.end(); trackIter++)
    {
        Event::CalMipTrack calMipTrack = *trackIter;

        m_log << MSG::DEBUG << "readCalMipTrack - MipTrack No " << counter << "-------" << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - Nh  =" << calMipTrack.getNh()        << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - PosX=" << calMipTrack.getPoint().x() << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - PosY=" << calMipTrack.getPoint().y() << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - PosZ=" << calMipTrack.getPoint().z() << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - DirX=" << calMipTrack.getDir().x()   << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - DirY=" << calMipTrack.getDir().y()   << endreq;
        m_log << MSG::DEBUG << "readCalMipTrack - DirZ=" << calMipTrack.getDir().z()   << endreq;

        Event::CalMipXtalVec calMipXtalVec = calMipTrack;
        readCalMipXtals(calMipXtalVec);

        counter++;
    }
        
    m_log << MSG::DEBUG << "readCalMipTracks -end-" << endreq;
    return;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::clearCalMipXtalVec(Event::CalMipXtalVec *calMipXtalVec)
{
    //    m_log << MSG::DEBUG << "clearCalMipXtalVec in StdMipFindingTool" << endreq;
    
    calMipXtalVec->clear();
    
    //    m_log << MSG::DEBUG << "clearCalMipXtalVec -end-" << endreq;
    return;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::clearCalMipTrackVec(Event::CalMipTrackVec* calMipTrackVec)
{
    //    m_log << MSG::DEBUG << "clearCalMipTrackVec in StdMipFindingTool" << endreq;
    
    calMipTrackVec->clear();

    //    m_log << MSG::DEBUG << "clearCalMipXtalVec -end-" << endreq;
    return;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipTracks()
{
    m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    m_log << MSG::DEBUG << "findMipTracks in StdMipFindingTool" << endreq;
    
    while (findC0())
    {
        if (!findC1())
        {
            m_calMipXtalVec[m_hid0].setFree(true);
            continue;
        }
    

        if (!findC2())
        {
            m_calMipXtalVec[m_hid0].setFree(true);
            m_calMipXtalVec[m_hid1].setFree(true);
            continue;
        }

        if (!leastSquares())
        {
            m_calMipXtalVec[m_hid0].setFree(true);
            m_calMipXtalVec[m_hid1].setFree(true);
            m_calMipXtalVec[m_hid2].setFree(true);
            m_calMipTrackVec.pop_back();
            continue;
        }

        while (findCn()) leastSquares();
    
        if (m_calMipTrackVec.back().getNh()<4)
        {
            for (int ih=0; ih<m_calMipTrackVec.back().getNh(); ih++) m_calMipXtalVec[ih].setFree(true);
            m_calMipTrackVec.pop_back();
        }
    }//end of track finding
    
    int numTracks=m_calMipTrackVec.size();
    if (numTracks>0) trackProperties();

    m_log << MSG::DEBUG << "findMipTracks -end- numTracks=" << m_calMipTrackVec.size() << endreq;
    return numTracks;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::trackProperties()
{
    //   m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    //   m_log << MSG::DEBUG << "trackProperties in StdMipFindingTool" << endreq;

    //@@@@@@@@@@@@@@@@@@@@@@ loop over tracks to compute their properties
    for(Event::CalMipTrackVec::iterator trackIter=m_calMipTrackVec.begin(); trackIter != m_calMipTrackVec.end(); trackIter++)
    {
        Event::CalMipTrack& calMipTrack = *trackIter;

        // distance to centroid
        m_dir  = calMipTrack.getDir();
        m_refP = calMipTrack.getPoint();
        m_uu   = m_refP-m_singleClusterCentroid;
        m_vv   = m_uu.cross(m_dir);
        double d2C = sqrt(m_vv*m_vv);
        calMipTrack.setD2C(d2C);

        // energy, distance to closest edge and arcLen
        double Edep[16][8];
        double ArcL[16][8];
        for (int itow=0; itow<16; itow++)
            for (int ilay=0; ilay<8; ilay++)
            {
                Edep[itow][ilay]=0;
                ArcL[itow][ilay]=0;
            }

        //-- get arclen per layer
        Point exitP1;
        Point exitP2;
        if (m_dir.z()!=0)
        {
            exitP1=m_refP+((m_calZBot-m_refP.z())/m_dir.z())*m_dir;
            exitP2=m_refP+((m_calZTop-m_refP.z())/m_dir.z())*m_dir;
        }
        else if (m_dir.x()!=0)
        {
            exitP1=m_refP+((m_calXLo-m_refP.x())/m_dir.x())*m_dir;
            exitP2=m_refP+((m_calXHi-m_refP.x())/m_dir.x())*m_dir;
        }
        else if (m_dir.y()!=0)
        {
            exitP1=m_refP+((m_calYLo-m_refP.y())/m_dir.y())*m_dir;
            exitP2=m_refP+((m_calYHi-m_refP.y())/m_dir.y())*m_dir;
        }
        m_uu=exitP2-exitP1;
        m_vv=m_uu.unit();
        m_G4PropTool->setStepStart(exitP1,m_vv);
        m_G4PropTool->step(sqrt(m_uu*m_uu));
        // Now loop over the steps to extract the materials
        int numSteps = m_G4PropTool->getNumberSteps();
        idents::VolumeIdentifier volId;
        idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
        for(int istep=0; istep<numSteps; ++istep)
        {
            volId=m_G4PropTool->getStepVolumeId(istep);
            volId.prepend(prefix);
            if(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0)// in Xtal ? 
            {
                int itow=4*volId[1]+volId[2];
                int ilay=volId[4];
                ArcL[itow][ilay]+=m_G4PropTool->getStepArcLen(istep);
            }
        }

        //-- get sum of track hit energy per layer - and d2edge at the same time
        double d2edge=99999.;
        int Nh=calMipTrack.getNh();
        for (int ih=0; ih<Nh; ih++)
        {
            Event::CalMipXtal* calMipXtal0=new Event::CalMipXtal();
            *calMipXtal0=calMipTrack.at(ih);    

            double d=D2Edge(calMipXtal0);
            if (d<d2edge) d2edge=d;

            int itow = calMipXtal0->getXtal()->getPackedId().getTower();
            int ilay = calMipXtal0->getXtal()->getPackedId().getLayer();
            Edep[itow][ilay]+=calMipXtal0->getXtal()->getEnergy();
            delete calMipXtal0;
        }
        calMipTrack.setD2Edge(d2edge);
        // to be modified later
        calMipTrack.setCalEdge(-1);
        
        //-- get trackArcLen, trackEcor and trackEcorRms
        double trackEcor=0;
        double trackEcorRms=0;
        double trackArcLen=0.;
        int nlay=0;
        for (int itow=0; itow<16; itow++)
            for (int ilay=0; ilay<8; ilay++)
            {
                double edep=Edep[itow][ilay];
                if (edep>0)
                {
                    double arcl=ArcL[itow][ilay];
                    if (arcl>0)
                    {
                        nlay++;
                        trackArcLen+=arcl;
                        trackEcor+=edep/arcl;
                        trackEcorRms+=(edep*edep)/(arcl*arcl);
                    }
                }
            }
        if (nlay>0)
        {
            trackEcor/=nlay;
            trackEcorRms/=nlay;
            trackEcorRms=trackEcorRms-trackEcor*trackEcor;
            if (trackEcorRms>=0)
              trackEcorRms=sqrt(trackEcorRms);
            else
              trackEcorRms=0.;
        }
        calMipTrack.setEcor(trackEcor*m_CsIHeight);
        calMipTrack.setEcorRms(trackEcorRms*m_CsIHeight);
        calMipTrack.setArcLen(trackArcLen);      

        // total energy (including track energy) in a cylinder of radius 1Rm=3.8cm around track
        double Rm=38;
        double trackErm=0;

        //-- find start and end of track segment
        Point startPos;
        Point endPos;
        double dmax=-1.;
        for(Event::CalMipXtalVec::iterator xTalIter1  = calMipTrack.begin(); 
                                           xTalIter1 != calMipTrack.end(); 
                                           xTalIter1++)
        {
            Event::CalMipXtal& xTal1 = *xTalIter1;
            //++ get the vector of reconstructed position
            Point H1=xTal1.getXtal()->getPosition();
            //++ projection on the track
            Point P1=m_refP+((H1-m_refP)*m_dir)*m_dir;
            for(Event::CalMipXtalVec::iterator xTalIter2  = calMipTrack.begin(); 
                                               xTalIter2 != calMipTrack.end(); 
                                               xTalIter2++)
            {
                Event::CalMipXtal& xTal2 = *xTalIter2;
                //++ get the vector of reconstructed position
                Point H2=xTal2.getXtal()->getPosition();
                //++ projection on the track
                Point P2=m_refP+((H2-m_refP)*m_dir)*m_dir;
                m_uu=P2-P1;
                double dist=sqrt(m_uu*m_uu);
                if (dist>dmax)
                {
                    startPos=P1;
                    endPos=P2;
                    dmax=dist;
                }
            }
        }

        //-- loop over CalXtalRecdata
        Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(m_dataSvc, EventModel::CalRecon::CalXtalRecCol); 
        for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
        {
            Event::CalXtalRecData* xTalData = *xTalIter;
            Point H = xTalData->getPosition();
            Point P = m_refP+((H-m_refP)*m_dir)*m_dir;
            //++ add energy if projection on track is between start and end of track segment
            if ((P-startPos)*(P-endPos)<0)
            {
                m_uu=H-P;
                double d2H=sqrt(m_uu*m_uu);
                if (d2H<Rm) trackErm+=xTalData->getEnergy();
            }
        }
        calMipTrack.setErm(trackErm);
    }
    
    //   m_log << MSG::DEBUG << "trackProperties -end-" << endreq;
    return;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findC0()
{
  //m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  //m_log << MSG::DEBUG << "findC0 in StdMipFindingTool" << endreq;

    m_hid     = -1;
    m_hid0    = -1;
    double d2Cmax = -1.;

    for(Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;
        double ene=calMipXtal.getXtal()->getEnergy();
          
        //int itow0 = calMipXtal.getXtal()->getPackedId().getTower();
        //int ilay0 = calMipXtal.getXtal()->getPackedId().getLayer();
        //int icol0 = calMipXtal.getXtal()->getPackedId().getColumn();

        //double x=calMipXtal.getXtal()->getPosition().x();
        //double y=calMipXtal.getXtal()->getPosition().y();
        //double z=calMipXtal.getXtal()->getPosition().z();
        //double d=calMipXtal.getD2C();
      
//        m_log << MSG::DEBUG << "find C0 " << itow0 << " " << ilay0 << " " << icol0  << " " << x << " " << y << " " << z  << " " << d << " " << ene << endreq;

        if (calMipXtal.getFree() && calMipXtal.getFreeC0() && ene>m_e1 && ene<m_e2 && calMipXtal.getD2C()>d2Cmax)
          //        if (calMipXtal.getFree() && calMipXtal.getFreeC0() && calMipXtal.getD2C()>d2Cmax)
        {
            d2Cmax = calMipXtal.getD2C();
            m_hid0 = m_hid;
        }
    }

    if (m_hid0>=0)
    {
        Event::CalMipXtal& calMipXtalRef = m_calMipXtalVec[m_hid0];
        calMipXtalRef.setFree(false);
        calMipXtalRef.setFreeC0(false);

        m_refP = m_calMipXtalVec[m_hid0].getXtal()->getPosition();

        //int itow = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
        //int ilay = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();
        //int icol = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getColumn();
        //m_log << MSG::DEBUG << "findC0 -end- true " << itow << " " << ilay << " " << icol  << endreq;
        return true;
    }

    m_log << MSG::DEBUG << "findC0 -end- false" << endreq;
    return false;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findC1()
{
  //m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  //m_log << MSG::DEBUG << "findC1 in StdMipFindingTool" << endreq;

    m_hid  = -1;
    m_hid1 = -1;
    double dmin = 99999.;
    double d01;

    int itow0   = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
    int ilay0   = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();
    double ene0 = m_calMipXtalVec[m_hid0].getXtal()->getEnergy();

    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;

        if (!calMipXtal.getFree())
          continue;

        int itow1 = calMipXtal.getXtal()->getPackedId().getTower();
        int ilay1 = calMipXtal.getXtal()->getPackedId().getLayer();
        int icol1 = calMipXtal.getXtal()->getPackedId().getColumn();

        //don't want initial direction to be horizontal
        if (itow1==itow0 && ilay1==ilay0)
          continue;

        m_uu=calMipXtal.getXtal()->getPosition()-m_refP;
        m_vv=m_uu.unit();

        Point xStart=m_refP-m_CsILength*m_vv;//go backward to step forward and cross C0 (to compute also its arcLen)
        m_G4PropTool->setStepStart(xStart,m_vv);

        d01=sqrt(m_uu*m_uu);
        //m_log << MSG::DEBUG << "findC1 current C1 " << itow1 << " " << ilay1 << " " << icol1 << " d01=" << d01 << endreq;
        m_G4PropTool->step(d01+3*m_CsILength);

        // Now loop over the steps to extract the materials
        int numSteps = m_G4PropTool->getNumberSteps();
        //m_log << MSG::DEBUG << "findC1 G4prop numSteps=" << numSteps << endreq;
        idents::VolumeIdentifier volId;
        idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
        int itow=-1;
        int ilay=-1;
        int icol=-1;
        bool foundOneCalMipXtalNotAlreadyUsed=false;
        double arcLen0=0;
        double arcLen1=0;
        int hidp=-1;
        bool lnext=true;
        bool crossedC0=false;
        for(int istep=0; istep<numSteps && lnext; ++istep)
        {
            volId=m_G4PropTool->getStepVolumeId(istep);
            volId.prepend(prefix);
            //        m_log << MSG::DEBUG << "findC1 istep=" << istep << " volId.name=" << volId.name() << endreq;

            if(!(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0))// in Xtal ? 
                continue;

            int iitow=4*volId[1]+volId[2];
            int iilay=volId[4];
            int iicol=volId[6];
            int iiseg=volId[8];
            int hid=m_hitId[iitow][iilay][iicol];
            //m_log << MSG::DEBUG << "findC1 seg in Xtal itow/ilay/icol/hid/iseg " << iitow << "  " << iilay << "  " << iicol << " " << hid << " " << iiseg << endreq;

            // continue if C0 still not crossed
            if (!crossedC0 && hid!=m_hid0)
                continue;

            double arcLen_step = m_G4PropTool->getStepArcLen(istep);
            if (hid==m_hid0)
            {
                arcLen0+=arcLen_step;
                crossedC0=true;
                //m_log << MSG::DEBUG << "findC1 seg in Xtal C0 itow/ilay/icol/hid/iseg " << iitow << "  " << iilay << "  " << iicol << " " << hid << " " << iiseg << " arcLen0 " << arcLen0 << endreq;
                continue;
            }
            // at this point C0 has been crossed and the volume does not belong to C0 crystal

            // continue if still in same layer as C0
            if (iitow==itow0 && iilay==ilay0)
                continue;
        
            // stop if one log has been found and current log is different from it
            if (foundOneCalMipXtalNotAlreadyUsed && !(iitow==itow && iilay==ilay && iicol==icol))
            {
                lnext=false;
                continue;
            }

            // stop if in different layer from current C1 layer
            if (!(iitow==itow1 && iilay==ilay1))
            {
              lnext=false;
              continue;
            }

            // stop if not in the hit map (CalMipXtalVec)
            if (hid<0)
            {
              lnext=false;
              continue;
            }

            // stop if not free
            if (!m_calMipXtalVec[hid].getFree())
            {
              lnext=false;
              continue;
            }

            // at this point a free crystal has been found, in same layer as current C1
            // continue if different from C1
            if (!(iitow==itow1 && iilay==ilay1 && iicol==icol1))
            {
              continue;
            }
            else
            {
                itow=iitow;
                ilay=iilay;
                icol=iicol;
                hidp=hid;
                foundOneCalMipXtalNotAlreadyUsed=true;
                arcLen1+=arcLen_step;
                //m_log << MSG::DEBUG << "findC1 seg in Xtal itow/ilay/icol/hid/iseg " << iitow << "  " << iilay << "  " << iicol << " " << hid << " " << iiseg << " arc " << arcLen_step << endreq;
            }
        }
        
        if (hidp==m_hid && arcLen0>0 && arcLen1>0)
        {
            double ecor0=ene0*m_CsIHeight/arcLen0;
            double ene1=m_calMipXtalVec[hidp].getXtal()->getEnergy();
            double ecor1=ene1*m_CsIHeight/arcLen1;
            //m_log << MSG::DEBUG << "------> findC1 candidate in Xtal itow/ilay/icol " << itow << "  " << ilay << "  " << icol << " arc0 " << arcLen0  << " ecor0 " << ecor0 << " arc1 " << arcLen1  << " ecor1 " << ecor1 << endreq;
            
            if (ecor0>m_ecor1 && ecor0<m_ecor2 && ecor1>m_ecor1 && ecor1<m_ecor2)
            {
                if (d01<dmin)
                {
                    dmin   = d01;
                    m_hid1 = hidp;
                    // impose downwards direction
                    if (m_vv.z()>0)  m_vv=-1.*m_vv;
                    m_dir = m_vv;
                }
            }
        }
    }

    if (m_hid1>=0)
    {
        Event::CalMipXtal& calMipXtalRef = m_calMipXtalVec[m_hid1];
        calMipXtalRef.setFree(false);

        //int itow = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getTower();
        //int ilay = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getLayer();
        //int icol = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getColumn();
        //        m_log << MSG::DEBUG << "findC1 -end- true " << itow << " " << ilay << " " << icol  << endreq;
        return true;
    }

    m_log << MSG::DEBUG << "findC1 -end- false" << endreq;
    return false;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findC2()
{
  //m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
  //m_log << MSG::DEBUG << "findC2 in StdMipFindingTool" << endreq;

    Point  xStart = m_refP;
    Vector dir    = m_dir;
    m_hid2=propagate(xStart,dir);

    if (m_hid2>=0)
    {
        Event::CalMipXtal& calMipXtalRef_m_hid2 = m_calMipXtalVec[m_hid2];
        calMipXtalRef_m_hid2.setFree(false);

        m_calMipTrackVec.push_back(Event::CalMipTrack());
        
        Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

        Event::CalMipXtal& calMipXtalRef0 = m_calMipXtalVec[m_hid0];
        calMipTrack.push_back(calMipXtalRef0);
        calMipXtalRef0.writeOut(m_log);

        Event::CalMipXtal& calMipXtalRef1 = m_calMipXtalVec[m_hid1];
        calMipTrack.push_back(calMipXtalRef1);
        calMipXtalRef1.writeOut(m_log);

        Event::CalMipXtal& calMipXtalRef2 = m_calMipXtalVec[m_hid2];
        calMipTrack.push_back(calMipXtalRef2);
        calMipXtalRef2.writeOut(m_log);

        calMipTrack.setPoint(m_refP);
        calMipTrack.setDir(m_dir);

        //m_log << MSG::DEBUG << "findC2 -end- true" << endreq;
        //m_log << MSG::DEBUG << "findC2 -end- refP =" << m_refP.x() << " " << m_refP.y() << "  " << m_refP.z() << endreq;  
        //m_log << MSG::DEBUG << "findC2 -end- dir  =" << m_dir.x() << " " << m_dir.y() << "  " << m_dir.z() << endreq;  
        return true;
    }

    m_log << MSG::DEBUG << "findC2 -end- false" << endreq;
    return false;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findCn()
{
  //m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
  //m_log << MSG::DEBUG << "findCn in StdMipFindingTool" << endreq;

    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

    m_refP = calMipTrack.getPoint();
    m_dir  = calMipTrack.getDir();

    Point xStart = m_refP;
    Vector dir   = m_dir;
    m_hidn=propagate(xStart,dir);

    if (m_hidn>=0)
    {
        Event::CalMipXtal& calMipXtalRef_m_hidn = m_calMipXtalVec[m_hidn];
        calMipXtalRef_m_hidn.setFree(false);
        calMipTrack.push_back(calMipXtalRef_m_hidn);

        //        m_log << MSG::DEBUG << "findCn -end- true" << endreq;
        //int itow = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getTower();
        //int ilay = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getLayer();
        //int icol = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getColumn();
        //        m_log << MSG::DEBUG << "findCn -end- " << itow << " " << ilay << " " << icol  << endreq;
        return true;
    }

    m_log << MSG::DEBUG << "findCn -end- false" << endreq;
    return false;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::propagate(Point xStart, Vector dir)
{
  //m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
  //m_log << MSG::DEBUG << "propagate in StdMipFindingTool" << endreq;
//   m_log << MSG::DEBUG << "propagate xStart=" << xStart.x() << " " << xStart.y() << "  " << xStart.z() << endreq;  
//   m_log << MSG::DEBUG << "propagate dir   =" << dir.x() << " " << dir.y() << "  " << dir.z() << endreq;  
    int hidn=-1;
    int radn=-1;
    double dmin = 99999;
    // get one unit vector perpendicular to track direction
    double costheta = dir.z();
    double sintheta = sqrt(1.-costheta*costheta);
    double cosphi;
    if (sintheta!=0) cosphi=dir.x()/sintheta;
    else             cosphi=1;
    Vector p(costheta/cosphi, 0., -sintheta);
    p=p.unit();
//     int nRad=5;
//     int nAng1=3;
//     double deltaRadius=4;
    int nRad=2;
    int nAng1=3;
    double deltaRadius=10;
    
    for (int iRad=1; iRad<=nRad; iRad++)
    {
        double radius = iRad * deltaRadius;
        int      nAng = iRad * nAng1;

        for(int iAng=0; iAng<=nAng; iAng++)
        {
            if (iRad>1 && iAng==0) continue;
            Point x0=xStart;
            if(iAng>0)
            {
                // Compute new starting point onto of Cal
                double rotAng=2.*m_pi*(iAng-1)/nAng; 
                //m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ propagate G4prop iRad=" << iRad << " radius=" << radius << " iAng=" << iAng << " rotAng " << rotAng << endreq;
                CLHEP::HepRotation rot(dir,rotAng);
                // get unit vector delta perpendicular to track direction with variable cylindrical phi angle 
                Vector delta=rot*p;
                // get starting point on cylinder surface
                Point xI = xStart + radius*delta;
                double s;
                if (costheta!=0) s=(xStart.z() - xI.z())/costheta;
                else             s=0;
                Ray segmt(xI,dir); 
                x0=segmt.position(s);
            }
      
            // search for another hit forward and backward
            for(int iforward=0;iforward<2;iforward++)
            {
                Point exitP;
                if (dir.z()!=0)
                {
                    if (iforward==0)
                    {
                        m_uu=dir;// m_dir always downwards
                        exitP=x0+((m_calZBot-x0.z())/m_uu.z())*m_uu;        
                    }
                    else
                    {
                        m_uu=-dir;// m_dir always downwards
                        exitP=x0+((m_calZTop-x0.z())/m_uu.z())*m_uu;
                    }
                }
                else
                {
                    if (dir.x()!=0)
                    {
                        if (iforward==0) m_uu=dir;
                        else             m_uu=-dir;
                        if (m_uu.x()>0)
                            exitP=x0+((m_calXHi-x0.x())/m_uu.x())*m_uu;
                        else
                            exitP=x0+((m_calXLo-x0.x())/m_uu.x())*m_uu;
                    }
                    else if (dir.y()!=0)
                    {
                        if (iforward==0)
                            m_uu=dir;
                        else
                            m_uu=-dir;
                        if (m_uu.y()>0)
                            exitP=x0+((m_calYHi-x0.y())/m_uu.y())*m_uu;
                        else
                            exitP=x0+((m_calYLo-x0.y())/m_uu.y())*m_uu;
                    }
                    else
                        continue;
                }
                //m_log << MSG::DEBUG << "propagate " << " " << iRad << " " << iAng << " x0= " << x0.x() << " " << x0.y() << "  " << x0.z() << " dir= " << m_uu.x() << " " << m_uu.y() << "  " << m_uu.z() << endreq;  
                m_G4PropTool->setStepStart(x0,m_uu);
                m_uu=x0-exitP;
                m_G4PropTool->step(sqrt(m_uu*m_uu));          
                // Now loop over the steps to extract the materials
                int numSteps = m_G4PropTool->getNumberSteps();
                //m_log << MSG::DEBUG << "propagate G4prop iforward=" << iforward << " numSteps=" << numSteps << endreq;
                idents::VolumeIdentifier volId;
                idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
                int itow=-1;
                int ilay=-1;
                int icol=-1;
                bool foundOneCalMipXtalNotAlreadyUsed=false;
                double arcLen=0;
                int hidp=-1;
                bool lnext=true;
                for(int istep=0; istep<numSteps && lnext; ++istep)
                {
                    volId=m_G4PropTool->getStepVolumeId(istep);
                    volId.prepend(prefix);
                    //      m_log << MSG::DEBUG << "propagate istep "<<  istep << " volId.name " << volId.name() << endreq;
                    if(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0)// in Xtal ? 
                    {
                        int iitow=4*volId[1]+volId[2];
                        int iilay=volId[4];
                        int iicol=volId[6];
                        //          int iiseg=volId[8];
                        //          m_log << MSG::DEBUG << "propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << endreq;
                        // stop if current log is different from the log previously found
                        if (foundOneCalMipXtalNotAlreadyUsed && !(iitow==itow && iilay==ilay && iicol==icol))
                        {
                            //m_log << MSG::DEBUG << "propagate end of vol search 1" << endreq;
                            lnext=false;
                            continue;
                        }
                        int hid=m_hitId[iitow][iilay][iicol];
                        if (hid<0) // stop if current log is not in the hit map (CalMipXtalVec)
                        {
                            //m_log << MSG::DEBUG << "propagate end of vol search 2" << endreq;
                            lnext=false;
                            continue;
                        }
                        double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
            
                        if (m_calMipXtalVec[hid].getFree())// can have crossing tracks....
                        {
                            itow=iitow;
                            ilay=iilay;
                            icol=iicol;
                            hidp=hid;
                            foundOneCalMipXtalNotAlreadyUsed=true;
                            arcLen+=arcLen_step;
                            //m_log << MSG::DEBUG << "propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << " arc " << arcLen_step << endreq;
                        }
                        // else
                        //  m_log << MSG::DEBUG << "propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " hit already in track" << endreq;
                    }
                }

                if (hidp>=0 && arcLen>0)
                {
                    double ecor=m_calMipXtalVec[hidp].getXtal()->getEnergy()*m_CsIHeight/arcLen;
                    //m_log << MSG::DEBUG << "------> propagate " << iRad << " " << iAng << " x0= " << x0.x() << " " << x0.y() << "  " << x0.z() << " dir= " << m_uu.x() << " " << m_uu.y() << "  " << m_uu.z() << endreq;  
                    //m_log << MSG::DEBUG << "------> propagate candidate in Xtal itow/ilay/icol " << itow << "  " << ilay << "  " << icol << " arc " << arcLen << " ecor " << ecor << endreq;
                    if (ecor>m_ecor1 && ecor<m_ecor2)
                    {
                        m_uu = m_calMipXtalVec[hidp].getXtal()->getPosition()-xStart;
                        m_vv = dir.cross(m_uu);
                        double d2dir = sqrt(m_vv*m_vv);
                        if (d2dir<dmin)
                        {
                            dmin = d2dir;
                            hidn = hidp;
                            radn=radius;
                        }
                    }
                }
            }
        }
    }

    //m_log << MSG::DEBUG << "propagate -end- hidn=" << hidn << " radius=" << radn << endreq;
    return hidn;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::leastSquares()
{
//   m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
//   m_log << MSG::DEBUG << "leastSquares in StdMipFindingTool" << endreq;

    int icoord;
    double xmeas,ymeas,zmeas,err2;

    // sigma2.z() is useless here no matter its value.
    double cov_xz = 0;  // covariance x,z
    double cov_yz = 0;  // covariance y,z
    double mx=0;        // mean x
    double my=0;        // mean y
    double mz=0;        // mean z
    double mz1=0;       // mean z for x pos
    double mz2=0;       // mean z for y pos
    double norm=0;      // sum of weights
    double norm1=0;     // sum of weights for odd layers
    double norm2=0;     // sum of weights for even layers
    double var_z1=0;    // variance of z for odd layers 
    double var_z2=0;    // variance of z for even layers
    double var_x=0;
    double var_y=0;

    // number of hits with non-zero energy deposition in X and Y direction, respectively
    int nlx=0,nly=0;
      
    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();
    Event::CalMipXtal   calMipXtal;
    //  m_log << MSG::DEBUG << "leastSquares Nh=" << calMipTrack.getNh() << endreq;
    for (int ih=0; ih<calMipTrack.getNh(); ih++)
    {
        calMipXtal = calMipTrack.at(ih);

        Point x0   = calMipXtal.getXtal()->getPosition();
        double ene = calMipXtal.getXtal()->getEnergy();
        int ilay   = calMipXtal.getXtal()->getPackedId().getLayer();

        for(int k=0;k<=1;k++)
        {
            if (k==0)// transverse info
            {
                icoord=1-ilay%2;
                // crude estimation of lateral position error
                err2=pow(m_CsIWidth,2)/12.;
            }
            else if (k==1)// longitudinal info
            {
                icoord=ilay%2;        
                // take err=10mm*(11.2/E_MeV) (MPV=11.2MeV for vertical muons)
                //err2=12544/(ene*ene);
                err2=100.;
            }
            zmeas=x0.z();
            // calculate weighting coefficient
            // double w = ene*ene/err2;
            double w = 1./err2;
            mz+=w*zmeas;
            norm+=w;
            if(icoord==0)
            {
                xmeas=x0.x();
                nlx++;
                // calculate sums for least square linear fit in XZ plane
                cov_xz += w*xmeas*zmeas;
                var_z1 += w*zmeas*zmeas;
                mx     += w*xmeas;
                var_x  += w*xmeas*xmeas;
                mz1    += w*zmeas;
                norm1  += w;
            }
            else
            {
                ymeas=x0.y();
                nly++;
                // calculate sums for least square linear fit in YZ plane
                cov_yz += w*ymeas*zmeas;
                var_z2 += w*zmeas*zmeas;
                my     += w*ymeas;
                var_y  += w*ymeas*ymeas;
                mz2    += w*zmeas;
                norm2  += w;
            }
        }
    }

    // linear fit requires at least 3 hits in both XZ and YZ planes
    if(nlx<3 || nly<3) return false;
      
    mx /= norm1;
    var_x/=norm1;
    var_x-=mx*mx;
    mz1 /= norm1;
    cov_xz /= norm1;
    cov_xz -= mx*mz1;
    var_z1 /= norm1;
    var_z1 -= mz1*mz1;
      
    // protection against dividing by 0 in the next statement
    if(var_z1==0) return false;
      
    // Now we have cov(x,z) and var(z) we can
    // deduce slope in XZ plane
    double tgthx= cov_xz/var_z1;

    my /= norm2;
    var_y/=norm2;
    var_y-=my*my;
    mz2 /= norm2;
    cov_yz /= norm2;
    cov_yz -= my*mz2;
    var_z2 /= norm2;
    var_z2 -= mz2*mz2;
      
    // protection against dividing by 0 in the next statement
    if(var_z2==0) return false;
      
    // Now we have cov(y,z) and var(z) we can
    // deduce slope in YZ plane
    double tgthy = cov_yz/var_z2;
      
    mz/=norm;

    // combining slope in XZ and YZ planes to get normalized 3-vector
    // of particle direction
    double tgtheta_sqr = tgthx*tgthx+tgthy*tgthy;
    double costheta = 1/sqrt(1+tgtheta_sqr);

    // impose uz<0
    double ux=-costheta*tgthx;
    double uy=-costheta*tgthy;
    double uz=-costheta;
  
    m_refP.setX(mx);
    m_refP.setY(my);
    m_refP.setZ(mz);
    calMipTrack.setPoint(m_refP);

    m_dir.setX(ux);
    m_dir.setY(uy);
    m_dir.setZ(uz);
    if (m_dir.z()>0)
        m_dir=-1.*m_dir;
    calMipTrack.setDir(m_dir);
  
    double chi2x=(var_x+tgthx*tgthx*var_z1-2*tgthx*cov_xz)*norm1;
    double chi2y=(var_y+tgthy*tgthy*var_z2-2*tgthy*cov_yz)*norm2;
    double chi2=chi2x/(nlx-2)+chi2y/(nly-2);
    calMipTrack.setChi2(chi2);
    
    //m_log << MSG::DEBUG << "leastSquares -end- refP=" << m_refP.x() << " " << m_refP.y() << "  " << m_refP.z() << endreq;
    //m_log << MSG::DEBUG << "leastSquares -end- dir =" << m_dir.x() << " " << m_dir.y() << "  " << m_dir.z() << endreq;

    return true;
}

//-----------------------------------------------------------------------------------------------------------------
double StdMipFindingTool::D2Edge(Event::CalMipXtal *calMipXtal)
{
    double x = calMipXtal->getXtal()->getPosition().x();
    double y = calMipXtal->getXtal()->getPosition().y();
    double z = calMipXtal->getXtal()->getPosition().z();

    double dd;
    double d2Edge = 99999;
  
    dd = (x-m_calXLo)*(x-m_calXLo);
    if (dd<d2Edge) d2Edge = dd;

    dd = (x-m_calXHi)*(x-m_calXHi);
    if (dd<d2Edge) d2Edge = dd;

    dd = (y-m_calYLo)*(y-m_calYLo);
    if (dd<d2Edge) d2Edge = dd;

    dd = (y-m_calYHi)*(y-m_calYHi);
    if (dd<d2Edge) d2Edge = dd;

    dd = (z-m_calZBot)*(z-m_calZBot);
    if (dd<d2Edge) d2Edge = dd;

    dd = (z-m_calZTop)*(z-m_calZTop);
    if (dd<d2Edge) d2Edge = dd;

    return sqrt(d2Edge);
}

//-----------------------------------------------------------------------------------------------------------------
// Main method
StatusCode StdMipFindingTool::findMIPCandidates()
{
    // Presume success unless something bad happens
    StatusCode sc = StatusCode::SUCCESS;
        
    m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    m_log << MSG::DEBUG << "findMIPCandidates in StdMipFindingTool" << endreq;

    //
    // Task 0: CalMipTrackVec and CalMipXtalVec cleaners in case of interrupted simulation of previous event
    //
    clearCalMipXtalVec (&m_calMipXtalVec);
    clearCalMipTrackVec(&m_calMipTrackVec);

    //
    // Task 1: get single cluster centroid from TDS
    //
    getSingleClusterCentroid();

    //
    // Task 2: Selection of Mip Hits in CalMipXtalVec
    //
    int numMipXtals=findMipXtals();
    if (numMipXtals<4)//at least 4 mipXtals required
    {
        clearCalMipXtalVec(&m_calMipXtalVec);
        return sc;
    }
    m_log << MSG::DEBUG << "numMipXtals=" << numMipXtals  << endreq;
        
    //
    // Task 3: Mip Tracks Finder
    //
    int numTracks=findMipTracks();
    m_log << MSG::DEBUG << "numTracks=" << numTracks << endreq;

    //
    // Task 4: Store m_calMipTrackVec in TDS
    if (numTracks>0)
      sc=storeCalMipTracks();

    //
    // Task 5a: CalMipXtalVec  cleaner
    // Task 5b: CalMipTrackVec cleaner
    //
    //readCalMipXtals (m_calMipXtalVec);
    //readCalMipTracks(m_calMipTrackVec);
    clearCalMipXtalVec (&m_calMipXtalVec);
    clearCalMipTrackVec(&m_calMipTrackVec);
    //readCalMipXtals (m_calMipXtalVec);
    //readCalMipTracks(m_calMipTrackVec);

    //
    // Task 6: read CalMipTrackCol from TDS for debugging
    //
    // if(numTracks>0)
    //  readCalMipTrackCol();

    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::storeCalMipTracks()
{
    StatusCode sc = StatusCode::SUCCESS;

    m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    m_log << MSG::DEBUG << "storeCalMipTracks in StdMipFindingTool" << endreq;

    m_calMipTrackCol = SmartDataPtr<Event::CalMipTrackCol>(m_dataSvc,EventModel::CalRecon::CalMipTrackCol);

    // If no pointer then create it
    if (m_calMipTrackCol == 0)
    {
        m_calMipTrackCol = new Event::CalMipTrackCol();
        sc = m_dataSvc->registerObject(EventModel::CalRecon::CalMipTrackCol, m_calMipTrackCol);
        if (sc.isFailure()) throw GaudiException("Failed to create Cal Mip Track Collection!", name(), sc);
    }

    int nbStoredTracks=0;

    for(Event::CalMipTrackVec::iterator trackIter=m_calMipTrackVec.begin(); trackIter != m_calMipTrackVec.end(); trackIter++)
    {
        Event::CalMipTrack& calMipTrack = *trackIter;

        // Need to create a new CalMipTrack which will be "owned" by the TDS
        Event::CalMipTrack* newMipTrack = new Event::CalMipTrack();

        // Now copy to this new track. This should properly copy all elements (including inherited vector)
        *newMipTrack = calMipTrack;

        // Store in collection (and reliquish ownership of the object)
        m_calMipTrackCol->push_back(newMipTrack); 
        nbStoredTracks++;
    }

    m_log << MSG::DEBUG << "storeCalMipTracks -end- nbStoredTracks=" << nbStoredTracks << endreq;

    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readCalMipTrackCol()
{
    m_log << MSG::DEBUG << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    m_log << MSG::DEBUG << "readCalMipTrackCol in StdMipFindingTool" << endreq;

    Event::CalMipTrackCol* p_calMipTrackCol = SmartDataPtr<Event::CalMipTrackCol>(m_dataSvc, EventModel::CalRecon::CalMipTrackCol); 

    int counter=-1;

    if (p_calMipTrackCol!=0)
    {
        for(Event::CalMipTrackCol::const_iterator trackIter=p_calMipTrackCol->begin(); trackIter != p_calMipTrackCol->end(); trackIter++)
        {
            counter++;
            m_log << "counter col=" << counter << endreq;
            Event::CalMipTrack* p_calMipTrack = *trackIter;
            p_calMipTrack->writeOut(m_log);
        }
    }
    
    m_log << MSG::DEBUG << "readCalMipTrackCol -end- " << counter << endreq;
    return;
}
