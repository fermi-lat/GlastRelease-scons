#include "IMipFindingTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
//@@@FP 07/11/05
#include "GaudiKernel/IToolSvc.h"
//@@@FP 07/11/05
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/SmartRefVector.h"


#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalClusterTab.h"
#include "Event/Recon/CalRecon/CalMipClasses.h"
//@@@FP 07/09/05
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
//@@@FP 07/09/05

#include "TMath.h"
#include "TMinuit.h"

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
    int                      findMipXtals(const Event::CalXtalRecCol* calXtalRecCol, const Event::CalCluster* calCluster);
    int                      findMipTracks();
    void                     calculateTrackProperties();
    bool                     findC0();
    bool                     findC1();
    bool                     findC2();
    bool                     findCn();
    int                      propagate(Point x0, Vector dir);
    bool                     findOldC1();
    bool                     findOldC2();
    bool                     findOldCn();
    void                     fitTrack(double ki2Cut, int maxCall);
    void                     leastSquares();
    double                   D2Edge(Event::CalMipXtal *calMipXtal);
    void                     readMipXtals (Event::CalMipXtalVec calMipXtalVec);
    void                     clearCalMipXtalVec(Event::CalMipXtalVec *calMipXtalVec);
    void                     readMipTracks(Event::CalMipTrackVec calMipTrackVec);
    void                     clearCalMipTrackVec(Event::CalMipTrackVec *calMipTrackVec);
    void                     readCalMipTrackCol();
    StatusCode               readGlastDet();
    StatusCode               storeCalMipTracks(Event::CalMipTrackVec calMipTrackVec);

    /// These are used by MINUIT for fitting
    TMinuit*                 m_minuit;
    static void              fcn(int& npar, double* gin, double& f, double* par, int iflag);
    static double            Ecor(Event::CalMipXtal calMipXtal);
  
    /// Private data members
    /// Cut values for MIP energy
    double                   m_mipE1;
    double                   m_mipE2;
    double                   m_ecor1;
    double                   m_ecor2;
    double                   m_radToDeg;

    double                   m_ki2Cut1;
    int                      m_maxCall1;
    double                   m_ki2Cut2;
    int                      m_maxCall2;
    double                   m_distMaxBetweenHits1;
    double                   m_distMaxBetweenHits2;

    /// Pointer to the Gaudi data provider service
    DataSvc*                 m_dataSvc;

    IPropagator * m_G4PropTool; 

    /// the GlastDetSvc used for access to detector info
    static IGlastDetSvc*     m_detSvc;
  //@@@FP 07/10/05
  /// TkrGeometrySvc used for access to tracker geometry info
  ITkrGeometrySvc* m_geoSvc;
  //@@@FP 07/10/05

    static double            m_CsILength;
    static double            m_CsIWidth;
    static double            m_CsIHeight;

    static double            m_XtalXlo[16][8][12];
    static double            m_XtalXhi[16][8][12];
    static double            m_XtalYlo[16][8][12];
    static double            m_XtalYhi[16][8][12];
    static double            m_XtalZlo[16][8][12];
    static double            m_XtalZhi[16][8][12];

    static int               m_hitId[16][8][12];

    static double            m_calZTop;
    static double            m_calZBot;
    static double            m_xMinCalEdge;
    static double            m_xMaxCalEdge;
    static double            m_yMinCalEdge;
    static double            m_yMaxCalEdge;
    static double            m_zMinCalEdge;
    static double            m_zMaxCalEdge;

    static int               m_xNum;       ///< x tower number
    static int               m_eLATTowers; ///< the value of fLATObjects field, defining LAT towers 
    static int               m_eTowerCAL;  ///< value of fTowerObject field, defining cal. module 
    static int               m_eXtal;      ///< the value of fCellCmp field defining CsI crystal
    static int               m_nCsISeg;    ///< number of geometric segments per Xtal

    static double            m_pi;

    //  Mip Xtals candidate vector from recon hit list;
    static Event::CalMipXtalVec     m_calMipXtalVec;

    // Mip track candidate vector from Mip Xtals candidate vector
    static Event::CalMipTrackVec    m_calMipTrackVec;
    static Event::CalMipTrackCol*   m_calMipTrackCol;

    /// These will be accessed by minuit fit function
    static Point             m_simpCluCent;
    static int               m_nbTracks;

    static int               m_hid;
    static int               m_hid0;
    static int               m_hid1;
    static int               m_hid2;
    static int               m_hidn;

    static Vector            m_uu;
    static Vector            m_vv;
    static Vector            m_dir;

    static Point             m_refP;
    static double            m_refPz;
    static Point             m_C;

    static int               m_lfit;
    static int               m_ki2Type;
    static int               m_ki2p;
    static double            m_ki2;
    static bool              m_goodfit;
} ;

//-----------------------------------------------------------------------------------------------------------------
static ToolFactory<StdMipFindingTool> s_factory;
const IToolFactory& StdMipFindingToolFactory = s_factory;

//-----------------------------------------------------------------------------------------------------------------
// Define the static variables
IGlastDetSvc*         StdMipFindingTool::m_detSvc; 
double                StdMipFindingTool::m_CsILength;
double                StdMipFindingTool::m_CsIWidth;
double                StdMipFindingTool::m_CsIHeight;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalXlo[16][8][12];
double                StdMipFindingTool::m_XtalXhi[16][8][12];
double                StdMipFindingTool::m_XtalYlo[16][8][12];
double                StdMipFindingTool::m_XtalYhi[16][8][12];
double                StdMipFindingTool::m_XtalZlo[16][8][12];
double                StdMipFindingTool::m_XtalZhi[16][8][12];
int                   StdMipFindingTool::m_hitId[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_calZTop;
double                StdMipFindingTool::m_calZBot;
double                StdMipFindingTool::m_xMinCalEdge;
double                StdMipFindingTool::m_xMaxCalEdge;
double                StdMipFindingTool::m_yMinCalEdge;
double                StdMipFindingTool::m_yMaxCalEdge;
double                StdMipFindingTool::m_zMinCalEdge;
double                StdMipFindingTool::m_zMaxCalEdge;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_xNum;
int                   StdMipFindingTool::m_eLATTowers;
int                   StdMipFindingTool::m_eTowerCAL;
int                   StdMipFindingTool::m_eXtal;
int                   StdMipFindingTool::m_nCsISeg;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_pi;
//-----------------------------------------------------------------------------------------------------------------
Event::CalMipXtalVec  StdMipFindingTool::m_calMipXtalVec;
Event::CalMipTrackVec StdMipFindingTool::m_calMipTrackVec;
Event::CalMipTrackCol* StdMipFindingTool::m_calMipTrackCol;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_simpCluCent;
int                   StdMipFindingTool::m_nbTracks;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_hid;
int                   StdMipFindingTool::m_hid0;
int                   StdMipFindingTool::m_hid1;
int                   StdMipFindingTool::m_hid2;
int                   StdMipFindingTool::m_hidn;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_uu;
Vector                StdMipFindingTool::m_vv;
Vector                StdMipFindingTool::m_dir;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_refP;
double                StdMipFindingTool::m_refPz;
Point                 StdMipFindingTool::m_C;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_lfit;
int                   StdMipFindingTool::m_ki2Type;
int                   StdMipFindingTool::m_ki2p;
double                StdMipFindingTool::m_ki2;
bool                  StdMipFindingTool::m_goodfit;
//-----------------------------------------------------------------------------------------------------------------
StdMipFindingTool::StdMipFindingTool(const std::string & type, 
                                     const std::string & name,
                                     const IInterface * parent ) : AlgTool( type, name, parent )
{ 
    declareInterface<IMipFindingTool>(this) ; 

    declareProperty("MIPEne1",       m_mipE1         = 6.8      );
    declareProperty("MIPEne2",       m_mipE2         = 19.6     );
    //50%
    //     declareProperty("Ecor1",       m_ecor1         = 10.3      );
    //     declareProperty("Ecor2",       m_ecor2         = 13.1     );
    //25%
    //     declareProperty("Ecor1",       m_ecor1         = 10.0      );
    //     declareProperty("Ecor2",       m_ecor2         = 14.2     );
    //10%
    declareProperty("Ecor1",       m_ecor1         = 9.7      );
    declareProperty("Ecor2",       m_ecor2         = 15.6     );
    declareProperty("Ki2Cut1",       m_ki2Cut1       = 10.      );
    declareProperty("MaxCall1",      m_maxCall1      = 1000     );
    declareProperty("Ki2Cut2",       m_ki2Cut2       = 99999.   );
    declareProperty("MaxCall2",      m_maxCall2      = 2000     );
    declareProperty("DistMaxBetweenHits1",m_distMaxBetweenHits1 =   80.    );
    declareProperty("DistMaxBetweenHits2",m_distMaxBetweenHits2 =  160.    );
    declareProperty("KI2TYPE",       m_ki2Type       = 0        );
    declareProperty("LFIT",          m_lfit          = 0        );
    return;
}
    
// This function extracts geometry constants from xml
// file using GlastDetSvc
StatusCode StdMipFindingTool::initialize()
{
    MsgStream log(msgSvc(), name());
    StatusCode sc = StatusCode::SUCCESS;
    log<<MSG::INFO<<"BEGIN initialize()"<<endreq ;

    //Locate and store a pointer to the data service
    IService* iService = 0;
    if ((sc = serviceLocator()->getService("EventDataSvc", iService)).isFailure())
    {
        throw GaudiException("Service [EventDataSvc] not found", name(), sc);
    }
    m_dataSvc = dynamic_cast<DataSvc*>(iService);

    //Define Minuit
    if (m_lfit)
      {
	m_minuit = new TMinuit(4);
	m_minuit->SetFCN(fcn);
      }

    readGlastDet();

    return StatusCode::SUCCESS ;
}

//-----------------------------------------------------------------------------------------------------------------
StatusCode StdMipFindingTool::readGlastDet()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readGlastDet in StdMipFindingTool" << endreq;

    /// constants defining the position of the fields in VolumeIdentifier 
    enum {fLATObjects, fTowerY, fTowerX, fTowerObjects, fLayer,
      fMeasure, fCALXtal,fCellCmp, fSegment};     

//@@@FP 07/09/05
//     m_detSvc=0;
//     sc = service("GlastDetSvc", m_detSvc);
    //@@@FP 07/09/05
    
    // extracting int constants
    double value;  // intermediate variable for reading constants from
    // xml file as doubles and converting them to integer 
    typedef std::map<int*,std::string> PARAMAP;
    PARAMAP param; // map containing pointers to integer constants to be read
    // with their symbolic names from xml file used as a key 
    
    // filling the map with information on constants to be read 

    param[&m_xNum]      =  std::string("xNum");
    param[&m_eLATTowers]=  std::string("eLATTowers");
    param[&m_eTowerCAL] =  std::string("eTowerCAL");
    param[&m_eXtal]     =  std::string("eXtal");
    param[&m_nCsISeg]   =  std::string("nCsISeg");
    
    // find TkrGeometrySvc service
    if (service("TkrGeometrySvc", m_geoSvc, true).isFailure()){
      log << MSG::ERROR << "Couldn't find the TkrGeometrySvc!" << endreq;
      return StatusCode::FAILURE;
    }

    m_calZTop = m_geoSvc->calZTop();
    m_calZBot = m_geoSvc->calZBot();

    // find GlastDevSvc service
    if (service("GlastDetSvc", m_detSvc, true).isFailure()){
      log << MSG::ERROR << "Couldn't find the GlastDetSvc!" << endreq;
      return StatusCode::FAILURE;
    }

    // loop over all constants information contained in the map
    for(PARAMAP::iterator it=param.begin(); it!=param.end();it++){
        //  attempt to get the constant value via the method of GlastDetSvc
        if(!m_detSvc->getNumericConstByName((*it).second, &value)) {
            // if not successful - give the error message and return
            log << MSG::ERROR << " constant " <<(*it).second <<" not defined" << endreq;
            return StatusCode::FAILURE;
            //  if successful - fill the constant using the pointer from the map
        } else *((*it).first)= int(value);
    }
    
    m_detSvc->getNumericConstByName("CsILength", &m_CsILength);
    m_detSvc->getNumericConstByName("CsIWidth",  &m_CsIWidth);
    m_detSvc->getNumericConstByName("CsIHeight", &m_CsIHeight);
    
    m_xMinCalEdge =  99999;
    m_xMaxCalEdge = -99999;
    m_yMinCalEdge =  99999;
    m_yMaxCalEdge = -99999;
    m_zMinCalEdge =  99999;
    m_zMaxCalEdge = -99999;

    for (int itower=0;itower<16;itower++)
    {
        for (int ilayer=0;ilayer<8;ilayer++)
        {
            for (int icolumn=0;icolumn<12;icolumn++)
            {
                // create Volume Identifier for segment 0 of this crystal
	            idents::VolumeIdentifier segm0Id;
	            segm0Id.append(0);
	            segm0Id.append((int)(itower/m_xNum));
	            segm0Id.append(itower%m_xNum);
	            segm0Id.append(0);
	            segm0Id.append(ilayer);
	            segm0Id.append(ilayer%2); 
	            segm0Id.append(icolumn);
	            segm0Id.append(0);
	            segm0Id.append(0);
	
	            //      
	            HepTransform3D transf;
	            //get 3D transformation for segment 0 of this crystal
	            m_detSvc->getTransform3DByID(segm0Id,&transf);
	            //get position of the center of the segment 0
	            Vector vect0 = transf.getTranslation();
	            // create Volume Identifier for the last segment of this crystal
	            idents::VolumeIdentifier segm11Id;
	            // copy all fields from segm0Id, except segment number
	            for(int ifield = 0; ifield<fSegment; ifield++) segm11Id.append(segm0Id[ifield]);
	            segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
	            //get 3D transformation for the last segment of this crystal
	            m_detSvc->getTransform3DByID(segm11Id,&transf);
	            //get position of the center of the last segment
	            Vector vect11 = transf.getTranslation();
	            Point p0(0.,0.,0.);	      
	            // position of the crystal center
	            Point pCenter = p0+(vect0+vect11)*0.5; 

	            double delta_x=(ilayer%2==0) ? 0.5*m_CsILength : 0.5*m_CsIWidth;
	            double delta_y=(ilayer%2==1) ? 0.5*m_CsILength : 0.5*m_CsIWidth;

	            m_XtalXlo[itower][ilayer][icolumn]=pCenter.x()-delta_x;
	            m_XtalXhi[itower][ilayer][icolumn]=pCenter.x()+delta_x;
	            m_XtalYlo[itower][ilayer][icolumn]=pCenter.y()-delta_y;
	            m_XtalYhi[itower][ilayer][icolumn]=pCenter.y()+delta_y;
	            m_XtalZlo[itower][ilayer][icolumn]=pCenter.z()-0.5*m_CsIHeight;
	            m_XtalZhi[itower][ilayer][icolumn]=pCenter.z()+0.5*m_CsIHeight;
	
	            if (pCenter.x()<m_xMinCalEdge) m_xMinCalEdge=pCenter.x();
	            if (pCenter.x()>m_xMaxCalEdge) m_xMaxCalEdge=pCenter.x();
	            if (pCenter.y()<m_yMinCalEdge) m_yMinCalEdge=pCenter.y();
	            if (pCenter.y()>m_yMaxCalEdge) m_yMaxCalEdge=pCenter.y();
	            if (pCenter.z()<m_zMinCalEdge) m_zMinCalEdge=pCenter.z();
	            if (pCenter.z()>m_zMaxCalEdge) m_zMaxCalEdge=pCenter.z();
            }
        }
    }
  
    m_xMinCalEdge-=m_CsIWidth/2;
    m_yMinCalEdge-=m_CsIWidth/2;
    m_zMinCalEdge-=m_CsIHeight/2;
    m_xMaxCalEdge+=m_CsIWidth/2;
    m_yMaxCalEdge+=m_CsIWidth/2;
    m_zMaxCalEdge+=m_CsIHeight/2;
    
    log << MSG::DEBUG << " FP : readGlastDet m_calZTop "<< m_calZTop << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_calZBot "<< m_calZBot << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_xMinCalEdge "<< m_xMinCalEdge << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_xMaxCalEdge "<< m_xMaxCalEdge << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_yMinCalEdge "<< m_yMinCalEdge << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_yMaxCalEdge "<< m_yMaxCalEdge << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_zMinCalEdge "<< m_zMinCalEdge << endreq;
    log << MSG::DEBUG << " FP : readGlastDet m_zMaxCalEdge "<< m_zMaxCalEdge << endreq;


    IToolSvc* toolSvc = 0;
    if(service("ToolSvc", toolSvc, true).isFailure()) {
      log << MSG::ERROR << "Couldn't find the ToolSvc!" << endreq;
      return StatusCode::FAILURE;
    }
    if(!toolSvc->retrieveTool("G4PropagationTool", m_G4PropTool)) {
      log << MSG::ERROR << "Couldn't find the G4PropationTool!" << endreq;
      return StatusCode::FAILURE;
    }

    log << MSG::DEBUG << " SG : readGlastDet - End" << endreq;
    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
double StdMipFindingTool::Ecor(Event::CalMipXtal calMipXtal)
{
    int    itower  = calMipXtal.getXtal()->getPackedId().getTower();
    int    ilayer  = calMipXtal.getXtal()->getPackedId().getLayer();
    int    icolumn = calMipXtal.getXtal()->getPackedId().getColumn();

    double xlow    = m_XtalXlo[itower][ilayer][icolumn];
    double xhig    = m_XtalXhi[itower][ilayer][icolumn];
    double ylow    = m_XtalYlo[itower][ilayer][icolumn];
    double yhig    = m_XtalYhi[itower][ilayer][icolumn];
    double zlow    = m_XtalZlo[itower][ilayer][icolumn];
    double zhig    = m_XtalZhi[itower][ilayer][icolumn];

    bool l1        = false;
    bool l2        = false;
    double lambda;

    Point E1;
    Point E2;

    if (m_dir.z()!=0)
      {
        // zlo plane
        lambda=(zlow-m_refP.z())/m_dir.z();
        m_C=m_refP+lambda*m_dir;

        if (m_C.y()>ylow && m_C.y()<yhig && m_C.x()>xlow && m_C.x()<xhig)
	  if (!l1)
	    {
	      E1=m_C;
	      l1=true;
	    }
	  else
	    {
	      E2=m_C;
	      l2=true;
	    }

        // zhi plane
        lambda=(zhig-m_refP.z())/m_dir.z();
        m_C=m_refP+lambda*m_dir;

        if (m_C.y()>ylow && m_C.y()<yhig && m_C.x()>xlow && m_C.x()<xhig)
	  if (!l1)
	    {
	      E1=m_C;
	      l1=true;
	    }
	  else
	    {
	      E2=m_C;
	      l2=true;
	    }
      }

    if (m_dir.x()!=0 && !(l1 && l2))
      {
        //xlow plane
        lambda=(xlow-m_refP.x())/m_dir.x();
        m_C=m_refP+lambda*m_dir;
        if (m_C.y()>ylow && m_C.y()<yhig && m_C.z()>zlow && m_C.z()<zhig)
	  if (!l1)
	    {
	      E1=m_C;
	      l1=true;
	    }
	  else
	    {
	      E2=m_C;
	      l2=true;
	    }

        //xhig plane
        lambda=(xhig-m_refP.x())/m_dir.x();
        m_C=m_refP+lambda*m_dir;

        if (m_C.y()>ylow && m_C.y()<yhig && m_C.z()>zlow && m_C.z()<zhig)
	  if (!l1)
	    {
	      E1=m_C;
	      l1=true;
	    }
	  else
	    {
	      E2=m_C;
	      l2=true;
	    }
      }

    if (m_dir.y()!=0 && !(l1 && l2))
      {
        //ylow plane
        lambda=(ylow-m_refP.y())/m_dir.y();
        m_C=m_refP+lambda*m_dir;

        if (m_C.x()>xlow && m_C.x()<xhig && m_C.z()>zlow && m_C.z()<zhig)
	  if (!l1)
	    {
	      E1=m_C;
	      l1=true;
	    }
	  else
	    {
	      E2=m_C;
	      l2=true;
	    }
	// yhig plane
	lambda=(yhig-m_refP.y())/m_dir.y();
	m_C=m_refP+lambda*m_dir;
	
	if (m_C.x()>xlow && m_C.x()<xhig && m_C.z()>zlow && m_C.z()<zhig)
	  if (!l1)
	    {
	      E1=m_C;
	      l1=true;
	    }
	  else
	    {
	      E2=m_C;
	      l2=true;
	    }
      }

    double ec=0.0001;

    if (l1 && l2)
      {
        m_uu=E2-E1;
        ec=sqrt(m_uu*m_uu)/m_CsIHeight;
      }
    else
      ec=1.;

    return ec;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipXtals(const Event::CalXtalRecCol* calXtalRecCol, const Event::CalCluster* calCluster)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findMipXtals in StdMipFindingTool" << endreq;

    double emip1       = 4.0;
    double emip2       = 40.0;
    int    numMipXtals = 0;


    for (int itow=0; itow<16; itow++)
      for (int ilay=0; ilay<8; ilay++)
	for (int icol=0; icol<12; icol++)
	  m_hitId[itow][ilay][icol]=-1;

    m_simpCluCent=calCluster->getPosition();
    log << MSG::DEBUG << " FP : findMipXtals centroid " << m_simpCluCent.x() << " " << m_simpCluCent.y() << " " << m_simpCluCent.z() << " " << endreq;
    //Loop over crystals in the collection

    //    idents::CalXtalId      xTalId;

    for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
    {
        Event::CalXtalRecData* xTalData = *xTalIter;
        double                 xTalE    =  xTalData->getEnergy();
        if ( xTalE>emip1 && xTalE<emip2)
        {
            m_calMipXtalVec.push_back(Event::CalMipXtal());
	    
	    int itow=xTalData->getPackedId().getTower();
	    int ilay=xTalData->getPackedId().getLayer();
	    int icol=xTalData->getPackedId().getColumn();
	    m_hitId[itow][ilay][icol]=numMipXtals;

            numMipXtals++;

	    //            xTalId = xTalData->getPackedId();

            Event::CalMipXtal& mipXtal = m_calMipXtalVec.back();
            //Update the values for this MipXtal
            mipXtal.setFree(true);
            //@@@FP 07/11/05
	    mipXtal.setFreeC0(true);
            //@@@FP 07/11/05
	    mipXtal.setXtal(xTalData);
            //Update the distance between hit and single cluster centroid
            m_C  = xTalData->getPosition();
            m_uu = m_simpCluCent - m_C;
            double d2C=sqrt(m_uu*m_uu);
            mipXtal.setD2C(d2C);
        }
    }
    log << MSG::DEBUG << " SG : findMipXtals " << numMipXtals << endreq;
    return numMipXtals;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readMipXtals(Event::CalMipXtalVec calMipXtalVec)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readMipXtals in StdMipFindingTool" << endreq;
    int counter=0;

    //Loop over crystals in the calMipXtalVec
    for(Event::CalMipXtalVec::iterator xTalIter=calMipXtalVec.begin(); xTalIter != calMipXtalVec.end(); xTalIter++)
    {
        Event::CalMipXtal calMipXtal=*xTalIter;

//         log << MSG::DEBUG << "readMipXtals - MipXtal No" << counter << "-------" << endreq;
//         log << MSG::DEBUG << "readMipXtals - Free=" << calMipXtal.getFree()  << endreq;
//         log << MSG::DEBUG << "readMipXtals - D2C =" << calMipXtal.getD2C()   << endreq;
//         log << MSG::DEBUG << "readMipXtals - Ener=" << calMipXtal.getXtal()->getEnergy()   << endreq;
//         log << MSG::DEBUG << "readMipXtals - PosX=" << calMipXtal.getXtal()->getPosition().x()   << endreq;
//         log << MSG::DEBUG << "readMipXtals - PosY=" << calMipXtal.getXtal()->getPosition().y()   << endreq;
//         log << MSG::DEBUG << "readMipXtals - PosZ=" << calMipXtal.getXtal()->getPosition().z()   << endreq;
        
        //calMipXtal.print();
        counter++;
    }
        
    log << MSG::DEBUG << " SG : readMipXtals End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readMipTracks(Event::CalMipTrackVec calMipTrackVec)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readMipTracks in StdMipFindingTool" << endreq;
    int counter=0;
    //Loop over tracks in the calMipTrackVec
    for(Event::CalMipTrackVec::iterator xTalIter=calMipTrackVec.begin(); xTalIter != calMipTrackVec.end(); xTalIter++)
    {
        Event::CalMipTrack calMipTrack=*xTalIter;

        log << MSG::DEBUG << "readMipTrack - MipTrack No" << counter << "-------" << endreq;
        log << MSG::DEBUG << "readMipTrack - Nh  =" << calMipTrack.getNh()        << endreq;
        log << MSG::DEBUG << "readMipTrack - Ndof=" << calMipTrack.getNdof()      << endreq;
        log << MSG::DEBUG << "readMipTrack - Ki2 =" << calMipTrack.getKi2()       << endreq;
        log << MSG::DEBUG << "readMipTrack - PosX=" << calMipTrack.getPoint().x() << endreq;
        log << MSG::DEBUG << "readMipTrack - PosY=" << calMipTrack.getPoint().y() << endreq;
        log << MSG::DEBUG << "readMipTrack - PosZ=" << calMipTrack.getPoint().z() << endreq;
        log << MSG::DEBUG << "readMipTrack - DirX=" << calMipTrack.getDir().x()   << endreq;
        log << MSG::DEBUG << "readMipTrack - DirY=" << calMipTrack.getDir().y()   << endreq;
        log << MSG::DEBUG << "readMipTrack - DirZ=" << calMipTrack.getDir().z()   << endreq;

        Event::CalMipXtalVec calMipXtalVec = calMipTrack;
        readMipXtals(calMipXtalVec);

        counter++;
    }
        
    log << MSG::DEBUG << " SG : readMipTracks - End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::clearCalMipXtalVec(Event::CalMipXtalVec *calMipXtalVec)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : clearCalMipXtalVec in StdMipFindingTool" << endreq;
    
    // Why isn't this all that is needed?
    calMipXtalVec->clear();

    log << MSG::DEBUG << " SG : clearCalMipXtalVec End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::clearCalMipTrackVec(Event::CalMipTrackVec* calMipTrackVec)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : clearCalMipTrackVec in StdMipFindingTool" << endreq;
    
    calMipTrackVec->clear();

    log << MSG::DEBUG << " SG : clearCalMipXtalVec End" << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipTracks()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findMipTracks in StdMipFindingTool" << endreq;

    m_pi = TMath::Pi();
    m_radToDeg = 180. / m_pi;
    sc=readGlastDet();

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
	
//         if (!findCn())
// 	  {
// 	    m_calMipXtalVec[m_hid0].setFree(true);
// 	    m_calMipXtalVec[m_hid1].setFree(true);
// 	    m_calMipXtalVec[m_hid2].setFree(true);
// 	    m_calMipTrackVec.pop_back();
// 	    m_nbTracks--;
// 	    continue;
// 	  }
	
	if (m_lfit)
	  fitTrack(m_ki2Cut1, m_maxCall1);
	else
	  leastSquares();
	
        if (!m_goodfit)
	  {
	    m_calMipXtalVec[m_hid0].setFree(true);
	    m_calMipXtalVec[m_hid1].setFree(true);
	    m_calMipXtalVec[m_hid2].setFree(true);
	    //	    m_calMipXtalVec[m_hidn].setFree(true);
	    m_calMipTrackVec.pop_back();
	    m_nbTracks--;
	    continue;
	  }
	
	bool another=false;
        while (findCn())
	  {
	    if (!m_lfit)
	      leastSquares();
	    another=true;
	  }

         if (another && m_lfit)
	   fitTrack(m_ki2Cut2, m_maxCall2);
	
        if (m_calMipTrackVec.back().getNh()<4)
	  //if (m_calMipTrackVec.back().getNh()<4 && m_calMipTrackVec.back().getKi2()>=100)
	  {
	    for (int ih=0; ih<m_calMipTrackVec.back().getNh(); ih++)
	      m_calMipXtalVec[ih].setFree(true);
	    
            m_calMipTrackVec.pop_back();
	    m_nbTracks--;
	  }
      }//end of track finding
       
    if (m_nbTracks+1>0)
      calculateTrackProperties();

    log << MSG::DEBUG << " SG : findMipTracks - End" << endreq;
    return m_nbTracks = m_calMipTrackVec.size();
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::calculateTrackProperties()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " FP : calculateTrackProperties in StdMipFindingTool" << endreq;

    //@@@@@@@@@@@@@@@@@@@@@@ loop over tracks to compute their properties
    for (int itr=0; itr<=m_nbTracks; itr++)
      {
	Event::CalMipTrack& calMipTrack = m_calMipTrackVec[itr];
	
	// ki2
	// 	    double kidof=-5;
	// 	    if (calMipTrack.getNdof()>0)
	// 		kidof=calMipTrack.getKi2()/calMipTrack.getNdof();
	//test
	calMipTrack.setKi2(123);
	calMipTrack.setNdof(321);
	
	// distance to centroid
	m_dir = calMipTrack.getDir();
	m_C   = calMipTrack.getPoint();
	m_uu  = m_C-m_simpCluCent;
	m_vv  = m_uu.cross(m_dir);
	double d2C = sqrt(m_vv*m_vv);
	calMipTrack.setD2C(d2C);
	
	// energy, distance to closest edge and length
	calMipTrack.setD2Edge(99999);
	calMipTrack.setLength(-99999);
	double energy=0.;
	for (int ih=0; ih<calMipTrack.getNh(); ih++)
	  {
	    Event::CalMipXtal* calMipXtal_ih=new Event::CalMipXtal();
	    *calMipXtal_ih=calMipTrack.at(ih);		
	    Point H1=calMipXtal_ih->getXtal()->getPosition();
	    // projection on the track
	    Point P1=m_C+((H1-m_C)*m_dir)*m_dir;
	    for (int ihh=0; ihh<calMipTrack.getNh() && ihh!=ih; ihh++)
	      {
		Event::CalMipXtal* calMipXtal_ihh=new Event::CalMipXtal();
		*calMipXtal_ihh=calMipTrack.at(ihh);
		Point H2=calMipXtal_ihh->getXtal()->getPosition();
		// projection on the track
		Point P2=m_C+((H2-m_C)*m_dir)*m_dir;
		// distance between two projection points
		m_uu=P2-P1;
		double d=sqrt(m_uu*m_uu);
		if (d>calMipTrack.getLength())
		  calMipTrack.setLength(d);
		delete calMipXtal_ihh;
	      }
	    energy+=calMipXtal_ih->getXtal()->getEnergy()/Ecor(*calMipXtal_ih);
	    double dd=D2Edge(calMipXtal_ih);
	    if (dd<calMipTrack.getD2Edge())
	      calMipTrack.setD2Edge(dd);
	    delete calMipXtal_ih;
	  }
	calMipTrack.setEnergy(energy);	    
      }
    
    log << MSG::DEBUG << " FP : calculateTrackProperties - End" << endreq;
    return;
}

//-----------------------------------------------------------------------------------------------------------------
  bool StdMipFindingTool::findC0()
{
    StatusCode sc = StatusCode::SUCCESS;
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

      if (calMipXtal.getFree() && calMipXtal.getFreeC0() && ene>m_mipE1 && ene<m_mipE2 && calMipXtal.getD2C() > d2Cmax)
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
	m_refPz = m_refP.z();

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
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    log << MSG::DEBUG << " SG : findC1 in StdMipFindingTool" << endreq;

    m_hid  = -1;
    m_hid1 = -1;
    double dmin = 99999.;
    double d01;

    int itow0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
    int ilay0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();
    int icol0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getColumn();

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
// 	    log << MSG::DEBUG << " FP : findC1 current C1 u.x= " << m_uu.x() << endreq;
// 	    log << MSG::DEBUG << " FP : findC1 current C1 u.y= " << m_uu.y() << endreq;
// 	    log << MSG::DEBUG << " FP : findC1 current C1 u.z= " << m_uu.z() << endreq;
 	    if (m_uu.z()==0)// don't want initial direction to be horizontal
 	      continue;
	    m_vv=m_uu.unit();
// 	    log << MSG::DEBUG << " FP : findC1 current C1 v.x= " << m_vv.x() << endreq;
// 	    log << MSG::DEBUG << " FP : findC1 current C1 v.y= " << m_vv.y() << endreq;
// 	    log << MSG::DEBUG << " FP : findC1 current C1 v.z= " << m_vv.z() << endreq;
// 	    // to avoid "G4Propagator: in danger of stuck particle"
// 	    if (fabs(m_vv.x())<0.00001)
// 	      m_vv.setX(0.);
// 	    if (fabs(m_vv.y())<0.00001)
// 	      m_vv.setY(0.);
// 	    if (fabs(m_vv.z())<0.00001)
// 	      m_vv.setZ(0.);
	    m_G4PropTool->setStepStart(m_refP,m_vv);
	    d01=sqrt(m_uu*m_uu);
	    log << MSG::DEBUG << " FP : findC1 current C1 " << itow1 << " " << ilay1 << " " << icol1 << " d01=" << d01 << endreq;
	    m_G4PropTool->step(2*d01);
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
		    int iiseg=volId[8];
		    log << MSG::DEBUG << " FP : findC1 seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << endreq;
		    // stop if current log is different from the log previously found
		    if (foundOneCalMipXtalNotAlreadyUsed && !(iitow==itow && iilay==ilay && iicol==icol))
		      {
			log << MSG::DEBUG << " FP : findC1 end of vol search 1" << endreq;
			lnext=false;
			continue;
		      }
		    int hid=m_hitId[iitow][iilay][iicol];
		    if (hid<0) // not in the hit map (CalMipXtalVec)
		      {
			log << MSG::DEBUG << " FP : findC1 end of vol search 2" << endreq;
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
			log << MSG::DEBUG << " FP : findC1 seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << " arc " << arcLen_step << " rad " << radLen_step << endreq;
		      }
		    else
		      log << MSG::DEBUG << " FP : findC1 seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " hit already in track" << endreq;
		  }
	      }
	    
	    if (hidp==m_hid)
	      {
		double ene=m_calMipXtalVec[hidp].getXtal()->getEnergy();
		double ecor=ene*m_CsIHeight/arcLen;
		log << MSG::DEBUG << "------> FP : findC1 candidate in Xtal itow/ilay/icol " << itow << "  " << ilay << "  " << icol << " arc " << arcLen << " rad " << radLen << " ecor " << ecor << endreq;
		    
		//		if (ecor>m_ecor1 && ecor<m_ecor2)
		if (ene>m_mipE1 && ene<m_mipE2)
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
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " FP : findC2 in StdMipFindingTool" << endreq;

    Point xStart=m_refP;
    Vector dir=m_dir;
    log << MSG::DEBUG << " FP : findC2 before propagate" << endreq;
    m_hid2=propagate(xStart,dir);
    log << MSG::DEBUG << " FP : findC2 after propagate" << endreq;

    if (m_hid2>=0)
    {
        Event::CalMipXtal& calMipXtalRef_m_hid2 = m_calMipXtalVec[m_hid2];
        calMipXtalRef_m_hid2.setFree(false);

	m_nbTracks++;
	m_calMipTrackVec.push_back(Event::CalMipTrack());
	Event::CalMipTrack& calMipTrackRef = m_calMipTrackVec.back();
	calMipTrackRef.setNdof(-1);
	calMipTrackRef.setKi2(-1);
	
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
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " FP : findCn in StdMipFindingTool" << endreq;

    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

    m_refP = calMipTrack.getPoint();
    m_dir    = calMipTrack.getDir();

    Point xStart=m_refP;
    Vector dir=m_dir;
    log << MSG::DEBUG << " FP : findCn before propagate" << endreq;
    m_hidn=propagate(xStart,dir);
    log << MSG::DEBUG << " FP : findCn after propagate" << endreq;

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
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
  log << MSG::DEBUG << " FP : propagate in StdMipFindingTool" << endreq;
  log << MSG::DEBUG << " FP : propagate dir.x = " << dir.x() << endreq;
  log << MSG::DEBUG << " FP : propagate dir.y = " << dir.y() << endreq;
  log << MSG::DEBUG << " FP : propagate dir.z = " << dir.z() << endreq;  
  int hidn=-1;
  double dmin = 99999;
  // get one unit vector perpendicular to track direction
  double costheta = dir.z();
  double sintheta = sqrt(1.-costheta*costheta);
  double cosphi;
  if (sintheta!=0)
    cosphi= dir.x()/sintheta;
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
	log << MSG::DEBUG << " @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FP : propagate G4prop numSample = " << is << endreq;
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
	    double s = (xStart.z() - xI.z())/costheta;
	    Ray segmt(xI,dir); 
	    x0 = segmt.position(s);
	  }
	   
	// 	  log << MSG::DEBUG << " FP : propagate x0.x = " << x0.x() << endreq;
	// 	  log << MSG::DEBUG << " FP : propagate x0.y = " << x0.y() << endreq;
	// 	  log << MSG::DEBUG << " FP : propagate x0.z = " << x0.z() << endreq;

	// search for another hit forward and backward
	for(int iforward=0;iforward<2;iforward++)
	  {
	    double sstep;
	    Point exitP;
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
	    m_G4PropTool->setStepStart(x0,m_uu);
	    m_uu=x0-exitP;
	    m_G4PropTool->step(sqrt(m_uu*m_uu));	      
	    // Now loop over the steps to extract the materials
	    int numSteps = m_G4PropTool->getNumberSteps();
	    log << MSG::DEBUG << " FP : propagate G4prop iforward = " << iforward << " numSteps = " << numSteps << endreq;
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
		    int iiseg=volId[8];
		    log << MSG::DEBUG << " FP : propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << endreq;
		    // stop if current log is different from the log previously found
		    if (foundOneCalMipXtalNotAlreadyUsed && !(iitow==itow && iilay==ilay && iicol==icol))
		      {
			log << MSG::DEBUG << " FP : propagate end of vol search 1" << endreq;
			lnext=false;
			continue;
		      }
		    int hid=m_hitId[iitow][iilay][iicol];
		    if (hid<0) // stop if current log is not in the hit map (CalMipXtalVec)
		      {
			log << MSG::DEBUG << " FP : propagate end of vol search 2" << endreq;
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
			log << MSG::DEBUG << " FP : propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " " << iiseg << " arc " << arcLen_step << " rad " << radLen_step << endreq;
		      }
		    else
		      log << MSG::DEBUG << " FP : propagate seg in Xtal itow/ilay/icol/iseg " << iitow << "  " << iilay << "  " << iicol << " hit already in track" << endreq;
		  }
	      }

	    if (hidp>=0)
	      {
		double ecor=m_calMipXtalVec[hidp].getXtal()->getEnergy()*m_CsIHeight/arcLen;
		log << MSG::DEBUG << "------> FP : propagate candidate in Xtal itow/ilay/icol " << itow << "  " << ilay << "  " << icol << " arc " << arcLen << " rad " << radLen << " ecor " << ecor << endreq;
		
		if (ecor>m_ecor1 && ecor<m_ecor2)
		  {
		    m_uu = m_calMipXtalVec[hidp].getXtal()->getPosition()-xStart;
		    m_vv = m_dir.cross(m_uu);
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
  log << MSG::DEBUG << " FP : end of propagate in StdMipFindingTool " << endreq;
  return hidn;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::leastSquares()
{
  StatusCode sc = StatusCode::SUCCESS;
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
  log << MSG::DEBUG << " FP : leastSquares in StdMipFindingTool" << endreq;

  m_goodfit=false;

  int icoord;
  double xmeas,ymeas,zmeas,err2;
  int kmin,kmax;
  if (m_ki2Type==0 || m_ki2Type==1)
    {
      kmin=0;
      kmax=1;
    }
  else if (m_ki2Type==2)
    {
      kmin=1;
      kmax=1;
    }
  else if (m_ki2Type==3)
    {
      kmin=0;
      kmax=0;
    }
  else
    return;

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
  
  // number of hits with non-zero energy deposition in X and Y direction, respectively
  int nlx=0,nly=0;
  
  // "non-physical vector of direction, which is returned
  // if fit is imposible due to insufficient number of hits
  
  Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();
  Event::CalMipXtal calMipXtal;
  log << MSG::DEBUG << " FP : leastSquares Nh = "<< calMipTrack.getNh() << endreq;
  for (int ih=0; ih<calMipTrack.getNh(); ih++)
    {
      calMipXtal = calMipTrack.at(ih);

      m_uu     = calMipXtal.getXtal()->getPosition();
      int ilay = calMipXtal.getXtal()->getPackedId().getLayer();

      for(int k=kmin;k<=kmax;k++)
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
	      mx += w*xmeas;
	      mz1 += w*zmeas;
	      norm1 += w;
	    }
	  else
	    {
	      ymeas=m_uu.y();
	      nly++;
	      // calculate sums for least square linear fit in YZ plane
	      cov_yz += w*ymeas*zmeas;
	      var_z2 += w*zmeas*zmeas;
	      my += w*ymeas;
	      mz2 += w*zmeas;
	      norm2 += w;
	    }
	}
    }

  // linear fit requires at least 2 hits in both XZ and YZ planes
  if(nlx <2 || nly < 2 ) return;
  
  mx /= norm1;
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
  m_refPz=mz;
  calMipTrack.setPoint(m_refP);

  m_dir.setX(ux);
  m_dir.setY(uy);
  m_dir.setZ(uz);
  if (m_dir.z()>0)
    m_dir=-1.*m_dir;
  calMipTrack.setDir(m_dir);
  
//   calMipTrack.setNdof(m_ndof);
//   calMipTrack.setKi2(m_ki2);

  log << MSG::DEBUG << " FP : leastSquares - End" << endreq;
  
  return;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::fitTrack(double ki2Cut, int maxCall)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " SG : fitTrack in StdMipFindingTool" << endreq;

    m_ki2 = 0;
  
    // Initial direction
    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

    m_refP = calMipTrack.getPoint();
    m_dir    = calMipTrack.getDir();

    double theta = TMath::ACos(m_dir.z());
    double phi;
    if (m_dir.x() != 0) phi  = atan(m_dir.y()/m_dir.x());
    if (m_dir.x() < 0)  phi += m_pi;
    if (phi > m_pi)   phi -= 2 * m_pi;
  
    // Beginning of m_minuit
    m_minuit->mninit(5,5,7);
  
    double arglist[10];
    int    ierflg = 0;
  
    //@@@FP 06/30/05
    //     arglist[0] = -1;
    //     m_minuit->mnexcm("SET PRI", arglist, 1, ierflg);
    //     m_minuit->mnexcm("SET NOW", arglist, 1, ierflg);
    //@@@FP 06/30/05
    arglist[0] = 1;    // up=1 car chi2
    m_minuit->mnexcm("SET ERR", arglist, 1, ierflg);
    arglist[0] = 2;
    m_minuit->mnexcm("SET STR", arglist, 1, ierflg);
  
    double vstart[4];

    vstart[0] = m_refP.x();
    vstart[1] = m_refP.y();
    vstart[2] = theta;
    vstart[3] = phi;

    double step[4] = {.1,.1,0.002,0.002};
  
    m_minuit->mnparm(0, "Xentry", vstart[0], step[0],-99999.,99999.,ierflg);
    m_minuit->mnparm(1, "Yentry", vstart[1], step[1],-99999.,99999.,ierflg);
    m_minuit->mnparm(2, "Theta",  vstart[2], step[2], m_pi/2,  m_pi,ierflg);
    m_minuit->mnparm(3, "Phi",    vstart[3], step[3],  -m_pi,  m_pi,ierflg);
  
    // minimization
    arglist[0] = maxCall; // maxcalls
    arglist[1] = 1;   // tolerance
    m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
    double Xentry, XentryErr;
    double Yentry, YentryErr;
    double thetaErr,phiErr;

    m_minuit->GetParameter(0, Xentry, XentryErr);
    m_minuit->GetParameter(1, Yentry, YentryErr);
    m_minuit->GetParameter(2, theta , thetaErr );
    m_minuit->GetParameter(3, phi   , phiErr   );

    // Retrieve the chi-square and determine number of degrees of freedom
    int npar=m_minuit->GetNumFreePars();
    int ndof=m_ki2p-npar;

    log << MSG::DEBUG <<"m_ki2p ="<<m_ki2p<<endreq;
    log << MSG::DEBUG <<"npar ="<<npar<<endreq;
    log << MSG::DEBUG <<"ndof ="<<ndof<<endreq;

    double kidof=-5.;
    if (ndof>0)
      kidof=m_ki2/ndof;
  
    m_goodfit=false;

    log << MSG::DEBUG <<"kidof="<<kidof<<" - ki2Cut="<<ki2Cut<<endreq;

    if ((kidof>0 && kidof<ki2Cut) || (ndof==0 && m_ki2<ki2Cut))
    {
        m_goodfit = true;

        m_refP.setX(Xentry);
        m_refP.setY(Yentry);
        calMipTrack.setPoint(m_refP);

        m_dir.setX(sin(theta)*cos(phi));
        m_dir.setY(sin(theta)*sin(phi));
        m_dir.setZ(cos(theta));
        calMipTrack.setDir(m_dir);

        calMipTrack.setNdof(ndof);
        calMipTrack.setKi2(m_ki2);
    }

    arglist[0]=1;
    gMinuit->mnexcm("CLEAR",arglist,1,ierflg); 

    log << MSG::DEBUG <<"m_goodfit="<<m_goodfit<<endreq;
    log << MSG::DEBUG << " SG : fitTrack - End" << endreq;
    
    return;
}

//-----------------------------------------------------------------------------------------------------------------
double StdMipFindingTool::D2Edge(Event::CalMipXtal *calMipXtal)
{
    double x = calMipXtal->getXtal()->getPosition().x();
    double y = calMipXtal->getXtal()->getPosition().y();
    double z = calMipXtal->getXtal()->getPosition().z();

    double dd     = (x-m_xMinCalEdge)*(x-m_xMinCalEdge);
    double d2Edge = 99999;
  
    if (dd<d2Edge) d2Edge = dd;

    dd = (x-m_xMaxCalEdge)*(x-m_xMaxCalEdge);
    if (dd<d2Edge) d2Edge = dd;

    dd = (y-m_yMinCalEdge)*(y-m_yMinCalEdge);
    if (dd<d2Edge) d2Edge = dd;

    dd = (y-m_yMaxCalEdge)*(y-m_yMaxCalEdge);
    if (dd<d2Edge) d2Edge = dd;

    dd = (z-m_zMinCalEdge)*(z-m_zMinCalEdge);
    if (dd<d2Edge) d2Edge = dd;

    dd = (z-m_zMaxCalEdge)*(z-m_zMaxCalEdge);
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
        
    log << MSG::DEBUG << " SG : Entering findMIPCandidates" << endreq;
        
    //
    // Task 0: Retrieve CalXtalRecData objects
    //
    Event::CalXtalRecCol* calXtalRecCol = SmartDataPtr<Event::CalXtalRecCol>(m_dataSvc, EventModel::CalRecon::CalXtalRecCol); 

    //Get single cluster from TDS
    /// Pointer to the data service (valBase ????)
    Event::CalClusterCol* pCals = SmartDataPtr<Event::CalClusterCol>(m_dataSvc, EventModel::CalRecon::CalClusterCol);
    Event::CalCluster* calCluster = pCals->front();
        
    //
    // Task 1: Selection of Mip Hits in CalMipXtalVec
    //
    int numMipXtals = findMipXtals(calXtalRecCol, calCluster);
    //requires at least 4 mipXtals
    if (numMipXtals<4)
      {
	clearCalMipXtalVec(&m_calMipXtalVec);
	return sc;
      }
    log << MSG::DEBUG << " SG : numMipXtals=" << numMipXtals  << endreq;
        
    //
    // Task 2: Mip Tracks Finder
    //
    int numTracks = findMipTracks();
        
    //
    // Task 5: Store m_calMipTrackVec in TDS
    //
    sc=storeCalMipTracks(m_calMipTrackVec);
        
    //
    // Task 7: CalMipXtalVec Cleaner
    //
    clearCalMipXtalVec(&m_calMipXtalVec);
    //    readMipXtals(m_calMipXtalVec);
        
    //
    // Task 8: CalMipTrackVec Cleaner
    //
    clearCalMipTrackVec(&m_calMipTrackVec);
    //readMipTracks(m_calMipTrackVec);

    //
    // Task 9: check : read CalMipTrack & CalMIPsCol from TDS
    //
        readCalMipTrackCol();

    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::fcn (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    StatusCode sc = StatusCode::SUCCESS;
    //Calculate chisquare
    m_ki2p = 0;
    m_ki2  = 0.;

    double dx2,err2;

    m_refP.setX(par[0]);
    m_refP.setY(par[1]);
    m_refP.setZ(m_refPz);

    double ct,st,cp,sp;
    ct=cos(par[2]);
    st=sin(par[2]);
    cp=cos(par[3]);
    sp=sin(par[3]);
    m_dir.setX(st*cp);
    m_dir.setY(st*sp);
    m_dir.setZ(ct);

    Event::CalMipXtal      calMipXtal;
    Event::CalXtalRecData* xtalData;
    Point                  point;

    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

    for(int k=0;k<=m_ki2Type;k++)
    { 
        for (int ih=0; ih<calMipTrack.getNh(); ih++)
	    {
	        calMipXtal = calMipTrack.at(ih);
		if (k==0)
		  {
	            xtalData = calMipXtal.getXtal();
	            point    = xtalData->getPosition();
	            m_uu     = point-m_refP;
	            m_vv     = m_dir.cross(m_uu);
	            dx2      = m_vv*m_vv;
	            err2     = 289.;
		  }
	        else if (k==1)
	        {
	            double ec = Ecor(calMipXtal);
	            dx2  = pow(calMipXtal.getXtal()->getEnergy()/ec-11.3,2);
	            err2 = 5.;
	        }
	  
		m_ki2+=dx2/err2;
	        m_ki2p++;
	    }
    }

    f = m_ki2;

  return;
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
        for(Event::CalMipTrackCol::const_iterator calMipTrackIter=p_calMipTrackCol->begin(); calMipTrackIter != p_calMipTrackCol->end(); calMipTrackIter++)
        {
            log << "------------------------------------------------------------" << endreq;
            log << "counter col=" << counter << endreq;
            log << "----------------" << endreq;
            counter++;
            Event::CalMipTrack* p_calMipTrack    =  *calMipTrackIter;
            p_calMipTrack->writeOut(log);
            log << "------------------------------------------------------------" << endreq;
        }
    }
    log << MSG::DEBUG << " SG : readCalMipTrackCol End " << counter << endreq;
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findOldC1()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    log << MSG::DEBUG << " SG : findOldC1 in StdMipFindingTool" << endreq;

    m_hid  = -1;
    m_hid1 = -1;
    double dmin = 99999.;

    Vector u01,ddir;

    int itow0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
    int ilay0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();
    int icol0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getColumn();

    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;
        if (calMipXtal.getFree())
	    {
	        m_uu = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[m_hid0].getXtal()->getPosition();
		// impose downwards direction
		if (m_uu.z()>0)
		  m_uu=-1.*m_uu;
	        u01 = m_uu.unit();
	        m_dir = u01;

	        int itow1 = calMipXtal.getXtal()->getPackedId().getTower();
	        int ilay1 = calMipXtal.getXtal()->getPackedId().getLayer();
	        int icol1 = calMipXtal.getXtal()->getPackedId().getColumn();
	  
	        if (itow1==itow0 && ilay1!=ilay0)
	        {
	            double ec0   = Ecor(m_calMipXtalVec[m_hid0]);
		    double ener0 = m_calMipXtalVec[m_hid0].getXtal()->getEnergy()/ec0;

		    double ec1   = Ecor(calMipXtal);
		    double ener1 = calMipXtal.getXtal()->getEnergy()/ec1;

		    if (ener0>m_mipE1 && ener0<m_mipE2 && ec0>0.001 &&
		        ener1>m_mipE1 && ener1<m_mipE2 && ec1>0.001)
		      {
		            double d01 = sqrt(m_uu*m_uu);
			    if (d01<dmin && d01<m_distMaxBetweenHits1)
			      {
		                dmin   = d01;
		                m_hid1 = m_hid;
				ddir = u01;
			      }
		      }
	        }
	    }
    }

    if (m_hid1>=0)
    {
        Event::CalMipXtal& calMipXtalRef = m_calMipXtalVec[m_hid1];
        calMipXtalRef.setFree(false);

 	m_dir = ddir;

	log << MSG::DEBUG << " SG : findOldC1 End - True" << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " SG : findOldC1 End - False" << endreq;
        return false;
    }

}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findOldC2()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " SG : findOldC2 in StdMipFindingTool" << endreq;

    m_hid  = -1;
    m_hid2 = -1;
    double dmin = 99999;

    Vector u01,u02,u12;

    int itow1 = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getTower();
    int ilay1 = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getLayer();
    int icol1 = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getColumn();

    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;

        if (calMipXtal.getFree())
	    {
	        int itow2 = calMipXtal.getXtal()->getPackedId().getTower();
	        int ilay2 = calMipXtal.getXtal()->getPackedId().getLayer();
	        int icol2 = calMipXtal.getXtal()->getPackedId().getColumn();

		if (itow2==itow1 && ilay2!=ilay1)
	        {
	            double ec=Ecor(calMipXtal);
	            double ener=calMipXtal.getXtal()->getEnergy()/ec;

		    if (ener>m_mipE1 && ener<m_mipE2 && ec>0.001)
		        {
		            u01 = m_calMipXtalVec[m_hid1].getXtal()->getPosition()-m_calMipXtalVec[m_hid0].getXtal()->getPosition();
		            u02 = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[m_hid0].getXtal()->getPosition();
			    double d02=sqrt(u02*u02);
		            u12 = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[m_hid1].getXtal()->getPosition();
			    double d12=sqrt(u12*u12);
			    double d2TrackHit=(d02<d12) ? d02 : d12;
		            m_vv  = u01.cross(u02);
		            double d2D01=sqrt(m_vv*m_vv);
		            if (d2D01<dmin && d2TrackHit<m_distMaxBetweenHits1)
		            {
		                dmin   = d2D01;
		                m_hid2 = m_hid;
		            }
		        }
	        }
	    }
    }
    
    log << MSG::DEBUG <<"m_hid2="<<m_hid2<<endreq;

    if (m_hid2>=0)
    {
        Event::CalMipXtal& calMipXtalRef_m_hid2 = m_calMipXtalVec[m_hid2];
        calMipXtalRef_m_hid2.setFree(false);

	m_nbTracks++;
	m_calMipTrackVec.push_back(Event::CalMipTrack());
	Event::CalMipTrack& calMipTrackRef = m_calMipTrackVec.back();
	calMipTrackRef.setNdof(-1);
	calMipTrackRef.setKi2(-1);
    
        Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

        Event::CalMipXtal& calMipXtalRef0 = m_calMipXtalVec[m_hid0];
        calMipTrack.push_back(calMipXtalRef0);

        Event::CalMipXtal& calMipXtalRef1 = m_calMipXtalVec[m_hid1];
        calMipTrack.push_back(calMipXtalRef1);

        Event::CalMipXtal& calMipXtalRef2 = m_calMipXtalVec[m_hid2];
        calMipTrack.push_back(calMipXtalRef2);

        calMipTrack.setPoint(m_refP);
        calMipTrack.setDir(m_dir);

        log << MSG::DEBUG << " SG : findOldC2 End - True" << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " SG : findOldC2 End - False" << endreq;
        return false;
    }
}

//-----------------------------------------------------------------------------------------------------------------
bool StdMipFindingTool::findOldCn()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " SG : findOldCn in StdMipFindingTool" << endreq;

    int itowl   = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getTower();
    int ilayl   = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getLayer();
    int icoll   = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getColumn();

    m_hidn = -1;
    m_hid  = -1;
    double dmin = 99999;

    Event::CalMipTrack& calMipTrack = m_calMipTrackVec.back();

    m_refP = calMipTrack.getPoint();
    m_dir    = calMipTrack.getDir();

    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;
      
        if (calMipXtal.getFree())
	    {
	      double ec = Ecor(calMipXtal);
	      double ener=calMipXtal.getXtal()->getEnergy()/ec;
	      if (ener>m_mipE1 && ener<m_mipE2 && ec>0.001)
		{
		  m_uu = calMipXtal.getXtal()->getPosition()-m_refP;
		  m_vv = m_dir.cross(m_uu);
		  double d2dir = sqrt(m_vv*m_vv);
		  
		  double d2TrackHit=99999;
		  Event::CalMipXtal calMipXtal_ih;
		  for (int ih=0; ih<calMipTrack.getNh(); ih++)
		    {
		      calMipXtal_ih = calMipTrack.at(ih);
		      m_uu=calMipXtal.getXtal()->getPosition()-calMipXtal_ih.getXtal()->getPosition();
		      double dd=sqrt(m_uu*m_uu);
		      if (dd<d2TrackHit)
			d2TrackHit=dd;
		    }
		  if (d2dir<dmin && d2TrackHit<m_distMaxBetweenHits2)
		    {
		      dmin   = d2dir;
		      m_hidn = m_hid;
		    }
		}
	    }
    }
    
    if (m_hidn>=0)
      {
        Event::CalMipXtal& calMipXtalRef_m_hidn = m_calMipXtalVec[m_hidn];
        calMipXtalRef_m_hidn.setFree(false);
        calMipTrack.push_back(calMipXtalRef_m_hidn);
        log << MSG::DEBUG << " SG : findOldCn End - True" << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " SG : findOldCn End - False" << endreq;
        return false;
    }
}
