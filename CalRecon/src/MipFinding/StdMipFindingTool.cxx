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
public:    StdMipFindingTool(const std::string & type, const std::string & name, const IInterface * parent );
  virtual ~StdMipFindingTool() {};
  
  /// @brief Intialization of the tool
  virtual StatusCode initialize();
  /// @brief Default cluster finding framework
  virtual StatusCode findMIPCandidates();
  
private:
  /// Private methods
  void                     getSingleClusterCentroid();
  int                      findMipXtals();
  int                      findMipTracks();
  void                     trackProperties();
  bool                     findC0();
  bool                     findC1();
  bool                     findC2();
  bool                     findCn();
  int                      propagate(Point x0, Vector dir);
  void                     leastSquares();
  double                   D2Edge             (Event::CalMipXtal     *calMipXtal);
  void                     readCalMipXtals    (Event::CalMipXtalVec   calMipXtalVec);
  void                     clearCalMipXtalVec (Event::CalMipXtalVec  *calMipXtalVec);
  void                     readCalMipTracks   (Event::CalMipTrackVec  calMipTrackVec);
  void                     clearCalMipTrackVec(Event::CalMipTrackVec *calMipTrackVec);
  void                     readCalMipTrackCol();
  StatusCode               readGlastDet();
  StatusCode               storeCalMipTracks(Event::CalMipTrackVec calMipTrackVec);

  /// Private data members
  /// Cut values for MIP energy
  double             m_e1;
  double             m_e2;
  double             m_ecor1;
  double             m_ecor2;
  double             m_pi;
  double             m_radToDeg;

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

  int                m_xNum;       ///< x tower number
  int                m_eLATTowers; ///< the value of fLATObjects field, defining LAT towers 
  int                m_eTowerCAL;  ///< value of fTowerObject field, defining cal. module 
  int                m_eXtal;      ///< the value of fCellCmp field defining CsI crystal
  int                m_nCsISeg;    ///< number of geometric segments per Xtal

  //  Mip Xtals candidate vector from recon hit list;
  Event::CalMipXtalVec     m_calMipXtalVec;
  // Mip track candidate vector from Mip Xtals candidate vector
  Event::CalMipTrackVec    m_calMipTrackVec;
  Event::CalMipTrackCol*   m_calMipTrackCol;

  Point              m_singleClusterCentroid;
  int                m_nbTracks;

  int                m_hid;
  int                m_hid0;
  int                m_hid1;
  int                m_hid2;
  int                m_hidn;

  Vector             m_uu;
  Vector             m_vv;
  Vector             m_dir;

  Point              m_refP;
  Point              m_C;

  bool               m_goodfit;
} ;

//-----------------------------------------------------------------------------------------------------------------
static ToolFactory<StdMipFindingTool> s_factory;
const IToolFactory& StdMipFindingToolFactory = s_factory;

//-----------------------------------------------------------------------------------------------------------------
StdMipFindingTool::StdMipFindingTool(const std::string & type, 
                                     const std::string & name,
                                     const IInterface * parent ) : AlgTool( type, name, parent )
{ 
    declareInterface<IMipFindingTool>(this) ; 

    declareProperty("MIPEne1",       m_e1         = 6.8      );
    declareProperty("MIPEne2",       m_e2         = 19.6     );
    //50% of Landau max
    //     declareProperty("Ecor1",       m_ecor1         = 10.3      );
    //     declareProperty("Ecor2",       m_ecor2         = 13.1     );
    //25%
    //     declareProperty("Ecor1",       m_ecor1         = 10.0      );
    //     declareProperty("Ecor2",       m_ecor2         = 14.2     );
    //10%
    declareProperty("Ecor1",       m_ecor1         = 9.7      );
    declareProperty("Ecor2",       m_ecor2         = 15.6     );
    return;
}
    
//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::initialize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;

    m_pi = TMath::Pi();
    m_radToDeg = 180. / m_pi;

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
      throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    m_dataSvc = dynamic_cast<DataSvc*>(iService);
    
    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
      log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
      return StatusCode::FAILURE;
    }

   // find TkrGeometrySvc service
    if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
      log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
      return StatusCode::FAILURE;
    }

    // find G4 propagation tool
    if(!toolSvc()->retrieveTool("G4PropagationTool", m_G4PropTool)) {
      log << MSG::ERROR << "Couldn't find the G4PropationTool!" << endreq;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS ;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::readGlastDet()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readGlastDet in StdMipFindingTool" << endreq;
    
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

    log << MSG::DEBUG << " FP : readGlastDet m_CsILength "<< m_CsILength << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_CsIWidth "<< m_CsIWidth << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_CsIHeight "<< m_CsIHeight << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calZTop "<< m_calZTop << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calZBot "<< m_calZBot << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calXLo "<< m_calXLo << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calXHi "<< m_calXHi << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calYLo "<< m_calYLo << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calYHi "<< m_calYHi << endreq;

    log << MSG::DEBUG << " SG : readGlastDet - End" << endreq;
    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::getSingleClusterCentroid()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findSingleClusterCentroid in StdMipFindingTool" << endreq;

    Event::CalClusterCol* pCals = SmartDataPtr<Event::CalClusterCol>(m_dataSvc, EventModel::CalRecon::CalClusterCol);
    Event::CalCluster* calCluster = pCals->front();

    m_singleClusterCentroid=calCluster->getPosition();
    log << MSG::DEBUG << " FP : findsingleClusterCentroid " << m_singleClusterCentroid.x() << " " << m_singleClusterCentroid.y() << " " << m_singleClusterCentroid.z() << " " << endreq;

    return;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipXtals()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findMipXtals in StdMipFindingTool" << endreq;

    double emip1    = 4.0;
    double emip2    = 40.0;
    int numMipXtals = 0;

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
            m_C  = xTalData->getPosition();
            m_uu = m_singleClusterCentroid - m_C;
            double d2C=sqrt(m_uu*m_uu);
            mipXtal.setD2C(d2C);
            mipXtal.setEcor(-1);
        }
    }
    log << MSG::DEBUG << " SG : findMipXtals " << numMipXtals << endreq;
    return numMipXtals;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readCalMipXtals(Event::CalMipXtalVec calMipXtalVec)
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readCalMipXtals in StdMipFindingTool" << endreq;
    int counter=0;

    //Loop over crystals in the calMipXtalVec
    for(Event::CalMipXtalVec::iterator xTalIter=calMipXtalVec.begin(); xTalIter != calMipXtalVec.end(); xTalIter++)
    {
        Event::CalMipXtal calMipXtal = *xTalIter;

//         log << MSG::DEBUG << "readCalMipXtals - MipXtal No" << counter << "-------" << endreq;
//         log << MSG::DEBUG << "readCalMipXtals - Free=" << calMipXtal.getFree()  << endreq;
//         log << MSG::DEBUG << "readCalMipXtals - D2C =" << calMipXtal.getD2C()   << endreq;
//         log << MSG::DEBUG << "readCalMipXtals - Ener=" << calMipXtal.getXtal()->getEnergy()   << endreq;
//         log << MSG::DEBUG << "readCalMipXtals - PosX=" << calMipXtal.getXtal()->getPosition().x()   << endreq;
//         log << MSG::DEBUG << "readCalMipXtals - PosY=" << calMipXtal.getXtal()->getPosition().y()   << endreq;
//         log << MSG::DEBUG << "readCalMipXtals - PosZ=" << calMipXtal.getXtal()->getPosition().z()   << endreq;
        
        //calMipXtal.print();
        counter++;
    }
        
    log << MSG::DEBUG << " SG : readCalMipXtals End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readCalMipTracks(Event::CalMipTrackVec calMipTrackVec)
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readCalMipTracks in StdMipFindingTool" << endreq;
    int counter=0;
    //Loop over tracks in the calMipTrackVec
    for(Event::CalMipTrackVec::iterator trackIter=calMipTrackVec.begin(); trackIter != calMipTrackVec.end(); trackIter++)
    {
        Event::CalMipTrack calMipTrack = *trackIter;

        log << MSG::DEBUG << "readCalMipTrack - MipTrack No" << counter << "-------" << endreq;
        log << MSG::DEBUG << "readCalMipTrack - Nh  =" << calMipTrack.getNh()        << endreq;
//         log << MSG::DEBUG << "readCalMipTrack - Ndof=" << calMipTrack.getNdof()      << endreq;
//         log << MSG::DEBUG << "readCalMipTrack - Ki2 =" << calMipTrack.getKi2()       << endreq;
        log << MSG::DEBUG << "readCalMipTrack - PosX=" << calMipTrack.getPoint().x() << endreq;
        log << MSG::DEBUG << "readCalMipTrack - PosY=" << calMipTrack.getPoint().y() << endreq;
        log << MSG::DEBUG << "readCalMipTrack - PosZ=" << calMipTrack.getPoint().z() << endreq;
        log << MSG::DEBUG << "readCalMipTrack - DirX=" << calMipTrack.getDir().x()   << endreq;
        log << MSG::DEBUG << "readCalMipTrack - DirY=" << calMipTrack.getDir().y()   << endreq;
        log << MSG::DEBUG << "readCalMipTrack - DirZ=" << calMipTrack.getDir().z()   << endreq;

        Event::CalMipXtalVec calMipXtalVec = calMipTrack;
        readCalMipXtals(calMipXtalVec);

        counter++;
    }
        
    log << MSG::DEBUG << " SG : readCalMipTracks - End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::clearCalMipXtalVec(Event::CalMipXtalVec *calMipXtalVec)
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : clearCalMipXtalVec in StdMipFindingTool" << endreq;
    
    calMipXtalVec->clear();

    log << MSG::DEBUG << " SG : clearCalMipXtalVec End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::clearCalMipTrackVec(Event::CalMipTrackVec* calMipTrackVec)
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : clearCalMipTrackVec in StdMipFindingTool" << endreq;
    
    calMipTrackVec->clear();

    log << MSG::DEBUG << " SG : clearCalMipXtalVec End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipTracks()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findMipTracks in StdMipFindingTool" << endreq;

    m_nbTracks=-1;

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
	
	leastSquares();
	
        if (!m_goodfit)
	  {
	    m_calMipXtalVec[m_hid0].setFree(true);
	    m_calMipXtalVec[m_hid1].setFree(true);
	    m_calMipXtalVec[m_hid2].setFree(true);
	    m_calMipTrackVec.pop_back();
	    m_nbTracks--;
	    continue;
	  }
	
        while (findCn())
	  leastSquares();
	
        if (m_calMipTrackVec.back().getNh()<4)
	  {
	    for (int ih=0; ih<m_calMipTrackVec.back().getNh(); ih++)
	      m_calMipXtalVec[ih].setFree(true);
            m_calMipTrackVec.pop_back();
	    m_nbTracks--;
	  }
      }//end of track finding
       
    if (m_nbTracks+1>0)
      trackProperties();

    log << MSG::DEBUG << " SG : findMipTracks - End" << endreq;
    return m_nbTracks=m_calMipTrackVec.size();
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::trackProperties()
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  log << MSG::DEBUG << " FP : trackProperties in StdMipFindingTool" << endreq;

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
    double trackEcor=0;
    double trackEcorRms=0;
    double trackArcLen=0.;
    double d2edge=99999.;
    int Nh=calMipTrack.getNh();
    for (int ih=0; ih<Nh; ih++)
    {
      Event::CalMipXtal* calMipXtal0=new Event::CalMipXtal();
      *calMipXtal0=calMipTrack.at(ih);	

      int itow = calMipXtal0->getXtal()->getPackedId().getTower();
      int ilay = calMipXtal0->getXtal()->getPackedId().getLayer();
      int icol = calMipXtal0->getXtal()->getPackedId().getColumn();
      int hid  = m_hitId[itow][ilay][icol];

      //log << MSG::DEBUG << " FP : trackProperties " << itow << " " << ilay << " " << icol << " " << hid << endreq;

      Event::CalMipXtal& calMipXtal = m_calMipXtalVec[hid];
      double s=100;
      Point x0=calMipXtal.getXtal()->getPosition()-s*m_dir;
            m_G4PropTool->setStepStart(x0,m_dir);
      m_G4PropTool->step(3*s);	      
      // Now loop over the steps to extract the materials
      int numSteps = m_G4PropTool->getNumberSteps();
      idents::VolumeIdentifier volId;
      idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
      double arcLen=0;
      bool lnext=true;
      bool found=false;
      for(int istep=0; istep<numSteps && lnext; ++istep)
      {
        volId=m_G4PropTool->getStepVolumeId(istep);
        volId.prepend(prefix);
        if(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0)// in Xtal ?	
        {
          int iitow=4*volId[1]+volId[2];
          int iilay=volId[4];
          int iicol=volId[6];
          int iiseg=volId[8];
          if (found && !(iitow==itow && iilay==ilay && iicol==icol))
          {
            lnext=false;
            continue;
          }
          if (iitow==itow && iilay==ilay && iicol==icol)
          {
            found=true;
            double arcLenStep=m_G4PropTool->getStepArcLen(istep); 
            arcLen+=arcLenStep;
            log << MSG::DEBUG << " FP : seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << " " << arcLenStep << endreq;
          }
        }
      }
      double ene =calMipXtal.getXtal()->getEnergy();
      double ecor=calMipXtal.getEcor();
      log << MSG::DEBUG << " FP : trackProperties ene / ecor / arcLen = " << ene << " " << ecor << " " << arcLen << endreq;        
     if (ecor==-1)
      {
        ecor=ene;
        if (arcLen!=0)
          ecor/=arcLen;
        
        //calMipXtal.setEcor(ecor);
        trackEcor+=ecor;
        trackEcorRms+=ecor*ecor;
        trackArcLen+=arcLen;
      }
      double d=D2Edge(calMipXtal0);
      if (d<d2edge)
        d2edge=d;
      delete calMipXtal0;
    }
    trackEcor/=Nh;
    trackEcorRms/=Nh;
    trackEcorRms=sqrt(trackEcorRms-trackEcor*trackEcor);
    calMipTrack.setEcor(trackEcor);
    calMipTrack.setEcorRms(trackEcorRms);
    calMipTrack.setArcLen(trackArcLen);      
    calMipTrack.setD2Edge(d2edge);

    //check
    calMipTrack.setCalEdge(-1);
    calMipTrack.setErm(-1);
  }
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findC0()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    log << MSG::DEBUG << " SG : findC0 in StdMipFindingTool" << endreq;

    m_hid     = -1;
    m_hid0    = -1;
    double d2Cmax = -1.;

    for(Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
      m_hid++;
      Event::CalMipXtal calMipXtal = *xTalIter;
      double ene=calMipXtal.getXtal()->getEnergy();
      
      int itow0 = calMipXtal.getXtal()->getPackedId().getTower();
      int ilay0 = calMipXtal.getXtal()->getPackedId().getLayer();
      int icol0 = calMipXtal.getXtal()->getPackedId().getColumn();

      double x=calMipXtal.getXtal()->getPosition().x();
      double y=calMipXtal.getXtal()->getPosition().y();
      double z=calMipXtal.getXtal()->getPosition().z();
      double d=calMipXtal.getD2C();
      
       log << MSG::DEBUG << " FP : find C0 " << itow0 << " " << ilay0 << " " << icol0  << " " << x << " " << y << " " << z  << " " << d << " " << ene << endreq;

      if (calMipXtal.getFree() && calMipXtal.getFreeC0() && ene>m_e1 && ene<m_e2 && calMipXtal.getD2C() > d2Cmax)
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

	m_refP  = m_calMipXtalVec[m_hid0].getXtal()->getPosition();

        log << MSG::DEBUG << " SG : findC0 End - True" << endreq;
	int itow = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
	int ilay = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();
	int icol = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getColumn();
	log << MSG::DEBUG << " FP : foundC0 " << itow << " " << ilay << " " << icol  << endreq;
        return true;
    }

    log << MSG::DEBUG << " SG : findC0 End - False" << endreq;
    return false;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findC1()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    log << MSG::DEBUG << " SG : findC1 in StdMipFindingTool" << endreq;

    m_hid  = -1;
    m_hid1 = -1;
    double dmin = 99999.;
    double d01;

    int itow0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
    int ilay0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();

    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
      {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;
	int itow1 = calMipXtal.getXtal()->getPackedId().getTower();
	int ilay1 = calMipXtal.getXtal()->getPackedId().getLayer();
	int icol1 = calMipXtal.getXtal()->getPackedId().getColumn();
        if (calMipXtal.getFree())
	  {
	    m_uu=calMipXtal.getXtal()->getPosition()-m_refP;
 	    if (m_uu.z()==0)// don't want initial direction to be horizontal
 	      continue;
	    m_vv=m_uu.unit();
	    m_G4PropTool->setStepStart(m_refP,m_vv);
	    d01=sqrt(m_uu*m_uu);
	    log << MSG::DEBUG << " FP : findC1 current C1 " << itow1 << " " << ilay1 << " " << icol1 << " d01=" << d01 << endreq;
	    m_G4PropTool->step(1.5*d01);
	    // Now loop over the steps to extract the materials
	    int numSteps = m_G4PropTool->getNumberSteps();
	    log << MSG::DEBUG << " FP : findC1 G4prop numSteps = " << numSteps << endreq;
	    idents::VolumeIdentifier volId;
	    idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
	    int itow=-1;
	    int ilay=-1;
	    int icol=-1;
	    bool foundOneCalMipXtalNotAlreadyUsed=false;
	    double arcLen=0;
	    double radLen=0;
	    int hidp=-1;
	    bool lnext=true;
	    for(int istep=0; istep<numSteps && lnext; ++istep)
	      {
		volId=m_G4PropTool->getStepVolumeId(istep);
		volId.prepend(prefix);
		//		log << MSG::DEBUG << " FP : findC1 istep "<<  istep << " volId.name " << volId.name() << endreq;
		if(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0)// in Xtal ?	
		  {
		    int iitow=4*volId[1]+volId[2];
		    int iilay=volId[4];
		    int iicol=volId[6];
                    //		    int iiseg=volId[8];
                    //		    log << MSG::DEBUG << " FP : findC1 seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << endreq;
		    // stop if current log is different from the log previously found
		    if (foundOneCalMipXtalNotAlreadyUsed && !(iitow==itow && iilay==ilay && iicol==icol))
		      {
			//log << MSG::DEBUG << " FP : findC1 end of vol search 1" << endreq;
			lnext=false;
			continue;
		      }
		    int hid=m_hitId[iitow][iilay][iicol];
		    if (hid<0) // not in the hit map (CalMipXtalVec)
		      {
			//log << MSG::DEBUG << " FP : findC1 end of vol search 2" << endreq;
			//			lnext=false;
			continue;
		      }
		    double radLen_step = m_G4PropTool->getStepRadLength(istep);
		    double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
		    //	    Point x_step       = m_G4PropTool->getStepPosition(istep);
		    // ignore logs found by G4prop in the same layer as C0
		    if (m_calMipXtalVec[hid].getFree() && !(iitow==itow0 && iilay==ilay0))
		      {
			itow=iitow;
			ilay=iilay;
			icol=iicol;
			hidp=hid;
			foundOneCalMipXtalNotAlreadyUsed=true;
			arcLen+=arcLen_step;
			radLen+=radLen_step;
                        //			log << MSG::DEBUG << " FP : findC1 seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << " arc " << arcLen_step << " rad " << radLen_step << endreq;
		      }
                    //		    else
                      //		      log << MSG::DEBUG << " FP : findC1 seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " hit already in track" << endreq;
		  }
	      }
	    
	    if (hidp==m_hid)
	      {
		double ene=m_calMipXtalVec[hidp].getXtal()->getEnergy();
		double ecor=ene*m_CsIHeight/arcLen;
		log << MSG::DEBUG << "------> FP : findC1 candidate in Xtal itow/ilay/icol " << itow << "  " << ilay << "  " << icol << " arc " << arcLen << " rad " << radLen << " ecor " << ecor << endreq;
		    
		//		if (ecor>m_ecor1 && ecor<m_ecor2)
		if (ene>m_e1 && ene<m_e2)
		  {
		    if (d01<dmin)
		      {
			dmin   = d01;
			m_hid1 = hidp;
			// impose downwards direction
			if (m_vv.z()>0)
			  m_vv=-1.*m_vv;
			m_dir = m_vv;
		      }
		  }
	      }
	  }
      }

    if (m_hid1>=0)
    {
        Event::CalMipXtal& calMipXtalRef = m_calMipXtalVec[m_hid1];
        calMipXtalRef.setFree(false);

	log << MSG::DEBUG << " SG : findC1 End - True" << endreq;
	int itow = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getTower();
	int ilay = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getLayer();
	int icol = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getColumn();
        log << MSG::DEBUG << " FP : foundC1 " << itow << " " << ilay << " " << icol  << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " SG : findC1 End - False" << endreq;
        return false;
    }

}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findC2()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " FP : findC2 in StdMipFindingTool" << endreq;

    Point xStart=m_refP;
    Vector dir=m_dir;
    // log << MSG::DEBUG << " FP : findC2 before propagate" << endreq;
    m_hid2=propagate(xStart,dir);
    // log << MSG::DEBUG << " FP : findC2 after propagate" << endreq;

    if (m_hid2>=0)
    {
        Event::CalMipXtal& calMipXtalRef_m_hid2 = m_calMipXtalVec[m_hid2];
        calMipXtalRef_m_hid2.setFree(false);

	m_nbTracks++;
	m_calMipTrackVec.push_back(Event::CalMipTrack());
        
//	Event::CalMipTrack& calMipTrackRef = m_calMipTrackVec.back();
        //check
// 	calMipTrackRef.setNdof(-1);
// 	calMipTrackRef.setKi2(-1);
	
        Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

        Event::CalMipXtal& calMipXtalRef0 = m_calMipXtalVec[m_hid0];
        calMipTrack.push_back(calMipXtalRef0);
	calMipXtalRef0.writeOut(log);

        Event::CalMipXtal& calMipXtalRef1 = m_calMipXtalVec[m_hid1];
        calMipTrack.push_back(calMipXtalRef1);
	calMipXtalRef1.writeOut(log);

        Event::CalMipXtal& calMipXtalRef2 = m_calMipXtalVec[m_hid2];
        calMipTrack.push_back(calMipXtalRef2);
	calMipXtalRef2.writeOut(log);

        calMipTrack.setPoint(m_refP);
        calMipTrack.setDir(m_dir);

        log << MSG::DEBUG << " FP : findC2 End - True" << endreq;
        log << MSG::DEBUG << " FP : findC2 End - True m_dir.x = " << m_dir.x() << endreq;
        log << MSG::DEBUG << " FP : findC2 End - True m_dir.y = " << m_dir.y() << endreq;
        log << MSG::DEBUG << " FP : findC2 End - True m_dir.z = " << m_dir.z() << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " FP : findC2 End - False" << endreq;
        return false;
    }
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findCn()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " FP : findCn in StdMipFindingTool" << endreq;

    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

    m_refP = calMipTrack.getPoint();
    m_dir  = calMipTrack.getDir();

    Point xStart=m_refP;
    Vector dir=m_dir;
    // log << MSG::DEBUG << " FP : findCn before propagate" << endreq;
    m_hidn=propagate(xStart,dir);
    // log << MSG::DEBUG << " FP : findCn after propagate" << endreq;

    if (m_hidn>=0)
      {
        Event::CalMipXtal& calMipXtalRef_m_hidn = m_calMipXtalVec[m_hidn];
        calMipXtalRef_m_hidn.setFree(false);
        calMipTrack.push_back(calMipXtalRef_m_hidn);
        log << MSG::DEBUG << " FP : findCn End - True" << endreq;
	int itow = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getTower();
	int ilay = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getLayer();
	int icol = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getColumn();
        log << MSG::DEBUG << " FP : foundCn " << itow << " " << ilay << " " << icol  << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " FP : findCn End - False" << endreq;
        return false;
    }
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::propagate(Point xStart, Vector dir)
{
  MsgStream log(msgSvc(), name());
//   log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
//   log << MSG::DEBUG << " FP : propagate in StdMipFindingTool" << endreq;
//   log << MSG::DEBUG << " FP : propagate dir.x = " << dir.x() << endreq;
//   log << MSG::DEBUG << " FP : propagate dir.y = " << dir.y() << endreq;
//   log << MSG::DEBUG << " FP : propagate dir.z = " << dir.z() << endreq;  
  int hidn=-1;
  double dmin = 99999;
  // get one unit vector perpendicular to track direction
  double costheta = dir.z();
  double sintheta = sqrt(1.-costheta*costheta);
  double cosphi;
  if (sintheta!=0)
    cosphi=dir.x()/sintheta;
  else
    cosphi=1;
  Vector p(costheta/cosphi, 0., -sintheta);
  p=p.unit();
  int numSamples=24;
  double deltaRadius=8;
  int numRadius=4;

  for (int numRad=1; numRad<=numRadius; numRad++)
    for(int is=numRad-1; is<=numSamples/numRad; is++)
      {
        //log << MSG::DEBUG << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FP : propagate G4prop numSample = " << is << endreq;
	Point x0=xStart;
	if(is>0)
	  {
	    // Compute new starting point onto of Cal
	    double rotAng = (is-1)*2.*m_pi/numSamples; 
	    HepRotation rot(dir, rotAng);
	    // get unit vector delta perpendicular to track direction with angle variable cylindrical phi angle 
	    Vector delta = rot*p;
	    // get starting point on cylinder surface
	    Point xI = xStart + (deltaRadius*numRad)*delta;
	    double s;
	    if (costheta!=0)
	      s=(xStart.z() - xI.z())/costheta;
	    else
	      s=0;
	    Ray segmt(xI,dir); 
	    x0=segmt.position(s);
	  }
	   
	// 	  log << MSG::DEBUG << " FP : propagate x0.x = " << x0.x() << endreq;
	// 	  log << MSG::DEBUG << " FP : propagate x0.y = " << x0.y() << endreq;
	// 	  log << MSG::DEBUG << " FP : propagate x0.z = " << x0.z() << endreq;

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
		    if (iforward==0)
		      m_uu=dir;
		    else
		      m_uu=-dir;
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
	      }
	    m_G4PropTool->setStepStart(x0,m_uu);
	    m_uu=x0-exitP;
	    m_G4PropTool->step(sqrt(m_uu*m_uu));	      
	    // Now loop over the steps to extract the materials
	    int numSteps = m_G4PropTool->getNumberSteps();
            //log << MSG::DEBUG << " FP : propagate G4prop iforward = " << iforward << " numSteps = " << numSteps << endreq;
	    idents::VolumeIdentifier volId;
	    idents::VolumeIdentifier prefix=m_detSvc->getIDPrefix();
	    int itow=-1;
	    int ilay=-1;
	    int icol=-1;
	    bool foundOneCalMipXtalNotAlreadyUsed=false;
	    double arcLen=0;
	    double radLen=0;
	    int hidp=-1;
	    bool lnext=true;
	    for(int istep=0; istep<numSteps && lnext; ++istep)
	      {
		volId=m_G4PropTool->getStepVolumeId(istep);
		volId.prepend(prefix);
		//		log << MSG::DEBUG << " FP : propagate istep "<<  istep << " volId.name " << volId.name() << endreq;
		if(volId.size()>7 && volId[0]==0 && volId[3]==0 && volId[7]==0)// in Xtal ?	
		  {
		    int iitow=4*volId[1]+volId[2];
		    int iilay=volId[4];
		    int iicol=volId[6];
                    //		    int iiseg=volId[8];
                    //		    log << MSG::DEBUG << " FP : propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << endreq;
		    // stop if current log is different from the log previously found
		    if (foundOneCalMipXtalNotAlreadyUsed && !(iitow==itow && iilay==ilay && iicol==icol))
		      {
                        //log << MSG::DEBUG << " FP : propagate end of vol search 1" << endreq;
			lnext=false;
			continue;
		      }
		    int hid=m_hitId[iitow][iilay][iicol];
		    if (hid<0) // stop if current log is not in the hit map (CalMipXtalVec)
		      {
                        //log << MSG::DEBUG << " FP : propagate end of vol search 2" << endreq;
			lnext=false;
			continue;
		      }
		    double radLen_step = m_G4PropTool->getStepRadLength(istep);
		    double arcLen_step = m_G4PropTool->getStepArcLen(istep); 
		    //	    Point x_step       = m_G4PropTool->getStepPosition(istep);
		    if (m_calMipXtalVec[hid].getFree())
		      {
			itow=iitow;
			ilay=iilay;
			icol=iicol;
			hidp=hid;
			foundOneCalMipXtalNotAlreadyUsed=true;
			arcLen+=arcLen_step;
			radLen+=radLen_step;
                        //log << MSG::DEBUG << " FP : propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << " arc " << arcLen_step << " rad " << radLen_step << endreq;
		      }
                    // else
                    // log << MSG::DEBUG << " FP : propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " hit already in track" << endreq;
		  }
	      }

	    if (hidp>=0)
	      {
		double ecor=m_calMipXtalVec[hidp].getXtal()->getEnergy()*m_CsIHeight/arcLen;
                //log << MSG::DEBUG << "------> FP : propagate candidate in Xtal itow/ilay/icol " << itow << "  " << ilay << "  " << icol << " arc " << arcLen << " rad " << radLen << " ecor " << ecor << endreq;
		
		if (ecor>m_ecor1 && ecor<m_ecor2)
		  {
		    m_uu = m_calMipXtalVec[hidp].getXtal()->getPosition()-xStart;
		    m_vv = dir.cross(m_uu);
		    double d2dir = sqrt(m_vv*m_vv);
		    if (d2dir<dmin)
		      {
			dmin = d2dir;
			hidn = hidp;
		      }
		  }
	      }
	  }
      }
  // log << MSG::DEBUG << " FP : end of propagate in StdMipFindingTool hidn= " << hidn << endreq;
  return hidn;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::leastSquares()
{
  MsgStream log(msgSvc(), name());
//   log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
//   log << MSG::DEBUG << " FP : leastSquares in StdMipFindingTool" << endreq;

  m_goodfit=false;

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
  log << MSG::DEBUG << " FP : leastSquares Nh = "<< calMipTrack.getNh() << endreq;
  for (int ih=0; ih<calMipTrack.getNh(); ih++)
    {
      calMipXtal = calMipTrack.at(ih);

      m_uu     = calMipXtal.getXtal()->getPosition();
      int ilay = calMipXtal.getXtal()->getPackedId().getLayer();

      for(int k=0;k<=1;k++)
	{
	  if (k==0)// transverse info
	    {
	      icoord=1-ilay%2;
	      err2=pow(m_CsIWidth,2)/12.;
	    }
	  else if (k==1)// longitudinal info
	    {
	      icoord=ilay%2;  	  
	      err2=289;
	    }
	  zmeas=m_uu.z();
	  // calculate weighting coefficient
	  double w = 1/err2;
	  mz+=w*zmeas;
	  norm+=w;
	  
	  if(icoord==0)
	    {
	      xmeas=m_uu.x();
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
	      ymeas=m_uu.y();
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
  if(nlx <3 || nly <3 ) return;
  
  mx /= norm1;
  var_x/=norm1;
  var_x-=mx*mx;
  mz1 /= norm1;
  cov_xz /= norm1;
  cov_xz -= mx*mz1;
  var_z1 /= norm1;
  var_z1 -= mz1*mz1;
  
  // protection against dividing by 0 in the next statment
  if(var_z1 == 0) return;
  
  // Now we have cov(x,z) and var(z) we can
  // deduce slope in XZ plane
  double tgthx = cov_xz/var_z1;
  my /= norm2;
  var_y/=norm2;
  var_y-=my*my;
  mz2 /= norm2;
  cov_yz /= norm2;
  cov_yz -= my*mz2;
  var_z2 /= norm2;
  var_z2 -= mz2*mz2;
  
  // protection against dividing by 0 in the next statment
  if(var_z2 == 0) return;
  
  m_goodfit=true;
  mz/=norm;

  // Now we have cov(y,z) and var(z) we can
  // deduce slope in YZ plane
  double tgthy = cov_yz/var_z2;
  
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
  
  log << MSG::DEBUG << " FP : leastSquares - End" << endreq;
  
  return;
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
    MsgStream log(msgSvc(), name());
        
    log << MSG::DEBUG << " SG : findMIPCandidates in StdMipFindingTool" << endreq;
    readGlastDet();

        
    //
    // Task 1: get single cluster centroid from TDS
    //
    getSingleClusterCentroid();

    //
    // Task 2: Selection of Mip Hits in CalMipXtalVec
    //
    int numMipXtals = findMipXtals();
    if (numMipXtals<4)//at least 4 mipXtals required
      {
	clearCalMipXtalVec(&m_calMipXtalVec);
	return sc;
      }
    log << MSG::DEBUG << " SG : numMipXtals=" << numMipXtals  << endreq;
        
    //
    // Task 3: Mip Tracks Finder
    //
    int numTracks = findMipTracks();
    log << MSG::DEBUG << " SG : numTracks=" << numTracks << endreq;

    //
    // Task 4: Store m_calMipTrackVec in TDS
    //
    sc=storeCalMipTracks(m_calMipTrackVec);
        
    //
    // Task 5: CalMipXtalVec Cleaner
    //
    clearCalMipXtalVec(&m_calMipXtalVec);
    //    readCalMipXtals(m_calMipXtalVec);
        
    //
    // Task 6: CalMipTrackVec Cleaner
    //
    clearCalMipTrackVec(&m_calMipTrackVec);
    //readCalMipTracks(m_calMipTrackVec);

    //
    // Task 7: check : read CalMipTrack & CalMIPsCol from TDS
    //
    readCalMipTrackCol();

    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::storeCalMipTracks(Event::CalMipTrackVec calMipTrackVec )
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : StoreCalMipTracks in StdMipFindingTool" << endreq;

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

    log << MSG::DEBUG << " SG : StoreCalMipTracks -  End" << endreq;

    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readCalMipTrackCol()
{
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readCalMipTrackCol in StdMipFindingTool" << endreq;

    Event::CalMipTrackCol* p_calMipTrackCol = SmartDataPtr<Event::CalMipTrackCol>(m_dataSvc, EventModel::CalRecon::CalMipTrackCol); 

    int counter=0;

    if (p_calMipTrackCol!=0)
    {
        for(Event::CalMipTrackCol::const_iterator trackIter=p_calMipTrackCol->begin(); trackIter != p_calMipTrackCol->end(); trackIter++)
        {
          log << "counter col=" << counter << endreq;
          counter++;
          Event::CalMipTrack* p_calMipTrack = *trackIter;
          p_calMipTrack->writeOut(log);
        }
    }
    log << MSG::DEBUG << " SG : readCalMipTrackCol End " << counter << endreq;
}
