#include "IMipFindingTool.h"
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/DataSvc.h"
#include "GaudiKernel/ToolFactory.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/GaudiException.h" 
#include "GaudiKernel/SmartRefVector.h"

#include "Event/TopLevel/EventModel.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"
#include "Event/Recon/CalRecon/CalCluster.h"
#include "Event/Recon/CalRecon/CalClusterTab.h"
#include "Event/Recon/CalRecon/CalMipClasses.h"
//#include "Event/Recon/CalRecon/CalMIPs.h"

//#include "Event/Recon/CalRecom/CalMipClasses.h" //if #include "Event/Recon/CalRecon/CalMIPs.h" is disabled

#include "geometry/Vector.h"

#include "TMath.h"
#include "TMinuit.h"

#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Transform3D.h"

//#include "Event/Recon/CalRecon/CalMipClasses.h"

/**   
* @class StdMipFindingTool
*
* $Header$
*/

//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------
// define these for now
double ASYMER2;
double TRANSER2;
//-----------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------

//
// ******************
//
// Now ready to define the tool
//
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
    int                      findMipXtals(const Event::CalXtalRecCol* calXtalRecCol, const Event::CalCluster* calCluster);
    int                      findMipTracks();
    bool                     findC0();
    bool                     findC1();
    bool                     findC2();
    void                     fitTrack(double ki2Cut, int maxCall);
    bool                     findCn();
    double                   D2Edge(Event::CalMipXtal *calMipXtal);
    void                     readMipXtals (Event::CalMipXtalVec calMipXtalVec);
    void                     clearCalMipXtalVec(Event::CalMipXtalVec *calMipXtalVec);
    void                     readMipTracks(Event::CalMipTrackVec calMipTrackVec);
    void                     clearCalMipTrackVec(Event::CalMipTrackVec *calMipTrackVec);
    StatusCode               readGlastDet();

    static double            Ecor(Event::CalMipXtal calMipXtal);

    StatusCode               storeCalMipTracks(Event::CalMipTrackVec calMipTrackVec);

    void                     readCalMipTrackCol();

    /// These are used by MINUIT for fitting
    //static TMinuit* m_minuit;
    TMinuit*                 m_minuit;
    static void              fcn(int& npar, double* gin, double& f, double* par, int iflag);
    static double            func(double* par, double z, int icoord);
  
    /// Private data members
    /// Cut values for MIP energy
    double                   m_mipEneCutLow;
    double                   m_mipEneCutHigh;
    double                   m_fThetaMin;
    double                   m_fThetaMax;
    double                   m_fPhiMin;
    double                   m_fPhiMax;
    double                   m_radToDeg;

    double                   m_thetaMin;
    double                   m_thetaMax;
    double                   m_phiMin;
    double                   m_phiMax;

    double                   m_ki2Cut1;
    int                      m_maxCall1;
    double                   m_ki2Cut2;
    int                      m_maxCall2;

    double                   m_dx_max;
    double                   m_beta0;

    /// These will be accessed by minuit fit function
    static double            m_EntryPTracy[3];
    static double            m_EntryPTracyErr[3];
    static int               m_ki2Type;
    static double            m_asymer2;
    static double            m_transer2;
    
    /// Pointer to the Gaudi data provider service
    DataSvc*                 m_dataSvc;

    //  Mip Xtals candidate vector from recon hit list;
    static Event::CalMipXtalVec     m_calMipXtalVec;

    // Mip track candidate vector from Mip Xtals candidate vector
    static Event::CalMipTrackVec    m_calMipTrackVec;

    static Event::CalMipTrackCol*   m_calMipTrackCol;

    static Point             m_simpCluCent;
    // List of variables for findC0, findC1, findC2, finCn
    static int               m_fcnNb;
    static bool              m_first;
    static int               m_hid;
    static int               m_hid0;
    static int               m_hid1;
    static int               m_hid2;
    static int               m_hidn;
    static double            m_d2C;
    static double            m_d2Cmax;

    static bool              m_free;
    static double            m_dmin;
    static double            m_dE;
    static double            m_xE;
    static Vector            m_u01;
    static Vector            m_u02;
    static Vector            m_u12;
    static Vector            m_vv;
    static Vector            m_vect;
    static double            m_d;
    static double            m_d1;
    static double            m_d2;
    static double            m_beta;
  //static double            m_beta0;
    static Point             m_D01point;
    static Vector            m_D01vector;
    static int               m_nbTracks;

    static Point             m_entryP;
    static double            m_entryPz;

    static Point             m_E1;
    static Point             m_E2;
    static Point             m_C;
    static Point             m_centroid;

    /// the GlastDetSvc used for access to detector info
    static IGlastDetSvc*     m_detSvc;
    static double            m_coordXtalZ0;
    static double            m_CsILength;
    static double            m_CsIWidth;
    static double            m_CsIHeight;

    static double            m_XtalXlo[16][8][12];
    static double            m_XtalXhi[16][8][12];
    static double            m_XtalYlo[16][8][12];
    static double            m_XtalYhi[16][8][12];
    static double            m_XtalZlo[16][8][12];
    static double            m_XtalZhi[16][8][12];

    static double            m_xMinCalEdge;
    static double            m_xMaxCalEdge;
    static double            m_yMinCalEdge;
    static double            m_yMaxCalEdge;
    static double            m_zMinCalEdge;
    static double            m_zMaxCalEdge;

    static double            m_ki2cut;
    static Vector            m_dir;
    static double            m_theta;
    static double            m_thetaErr;
    static double            m_phi;
    static double            m_phiErr;
    static double            m_pi;
    static double            m_ki2p;
    static double            m_Ki2Tot;
    static double            m_kidof;
    static double            m_ki2;
    static double            m_ki2Tot;
    static double            m_thetaf;
    static double            m_thetafErr;
    static double            m_phif;
    static double            m_phifErr;
    static double            m_Xentry;
    static double            m_XentryErr;
    static double            m_Yentry;
    static double            m_YentryErr;
    static double            m_Xentryf;
    static double            m_XentryfErr;
    static double            m_Yentryf;
    static double            m_YentryfErr;
    static double            m_kidoff;

    static int               m_xNum;       ///< x tower number
    static int               m_eLATTowers; ///< the value of fLATObjects field, defining LAT towers 
    static int               m_eTowerCAL;  ///< value of fTowerObject field, defining cal. module 
    static int               m_eXtal;      ///< the value of fCellCmp field defining CsI crystal
    static int               m_nCsISeg;    ///< number of geometric segments per Xtal

    static int               m_comptor;
    static double            m_chiSquare;
    static int               m_nPar;
    static int               m_ndof;
    static bool              m_goodfit;


    static int               m_compteurTest;
} ;
//-----------------------------------------------------------------------------------------------------------------
static ToolFactory<StdMipFindingTool> s_factory;
const IToolFactory& StdMipFindingToolFactory = s_factory;
//-----------------------------------------------------------------------------------------------------------------
// Define the static variables
double                StdMipFindingTool::m_EntryPTracy[3];
double                StdMipFindingTool::m_EntryPTracyErr[3];
int                   StdMipFindingTool::m_ki2Type;
double                StdMipFindingTool::m_asymer2;
double                StdMipFindingTool::m_transer2;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_xNum;
int                   StdMipFindingTool::m_eLATTowers;
int                   StdMipFindingTool::m_eTowerCAL;
int                   StdMipFindingTool::m_eXtal;
int                   StdMipFindingTool::m_nCsISeg;
Event::CalMipXtalVec  StdMipFindingTool::m_calMipXtalVec;
//-----------------------------------------------------------------------------------------------------------------
// Mip track candidate vector from Mip Xtals candidate vector
Event::CalMipTrackVec StdMipFindingTool::m_calMipTrackVec;
//-----------------------------------------------------------------------------------------------------------------
//Event::CalMIPsCol*    StdMipFindingTool::m_calMIPsCol;
//-----------------------------------------------------------------------------------------------------------------
Event::CalMipTrackCol* StdMipFindingTool::m_calMipTrackCol;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_fcnNb;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_simpCluCent;
//-----------------------------------------------------------------------------------------------------------------
bool                  StdMipFindingTool::m_first;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_hid;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_hid0;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_hid1;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_hid2;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_hidn;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_d2C;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_d2Cmax;
//-----------------------------------------------------------------------------------------------------------------
bool                  StdMipFindingTool::m_free;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_dmin;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_dE;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_xE;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_u01;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_u02;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_u12;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_vv;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_vect;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_d;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_d1;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_d2;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_beta;
//-----------------------------------------------------------------------------------------------------------------
//double                StdMipFindingTool::m_beta0;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_D01point;
Vector                StdMipFindingTool::m_D01vector;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_nbTracks;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_entryP;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_entryPz;
//-----------------------------------------------------------------------------------------------------------------
IGlastDetSvc*         StdMipFindingTool::m_detSvc; 
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_coordXtalZ0;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_CsILength;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_CsIWidth;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_CsIHeight;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_ki2cut;
//-----------------------------------------------------------------------------------------------------------------
Vector                StdMipFindingTool::m_dir;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_theta;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_thetaErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_phi;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_phiErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_pi;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_ki2p;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_Ki2Tot;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_kidof;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_ki2;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_thetaf;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_thetafErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_phif;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_phifErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_Xentry;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XentryErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_Yentry;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_YentryErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_Xentryf;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XentryfErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_Yentryf;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_YentryfErr;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_kidoff;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_comptor;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_chiSquare;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_nPar;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_ndof;
//-----------------------------------------------------------------------------------------------------------------
bool                  StdMipFindingTool::m_goodfit;
//-----------------------------------------------------------------------------------------------------------------
int                   StdMipFindingTool::m_compteurTest;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_xMinCalEdge;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_xMaxCalEdge;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_yMinCalEdge;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_yMaxCalEdge;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_zMinCalEdge;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_zMaxCalEdge;
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalXlo[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalXhi[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalYlo[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalYhi[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalZlo[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
double                StdMipFindingTool::m_XtalZhi[16][8][12];
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_E1;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_E2;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_C;
//-----------------------------------------------------------------------------------------------------------------
Point                 StdMipFindingTool::m_centroid;
//-----------------------------------------------------------------------------------------------------------------
StdMipFindingTool::StdMipFindingTool(const std::string & type, 
                                     const std::string & name,
                                     const IInterface * parent ) : AlgTool( type, name, parent )
{ 
    declareInterface<IMipFindingTool>(this) ; 

    declareProperty("MIPEneCutLow",  m_mipEneCutLow  = 2.       );
    declareProperty("MIPEneCutHigh", m_mipEneCutHigh = 50.      );
    declareProperty("KI2TYPE",       m_ki2Type       = 0        );
    declareProperty("Asymer2",       m_asymer2       = ASYMER2  );
    declareProperty("Transer2",      m_transer2      = TRANSER2 );
    declareProperty("ThetaMin",      m_fThetaMin     = -10.     );
    declareProperty("ThetaMax",      m_fThetaMax     =  90.     );
    declareProperty("PhiMin",        m_fPhiMin       = -180.    );
    declareProperty("PhiMax",        m_fPhiMax       =  180.    );
    declareProperty("ThetaMin",      m_thetaMin      = -10.     );
    declareProperty("ThetaMax",      m_thetaMax      =  90.     );
    declareProperty("PhiMin",        m_phiMin        = -180.    );
    declareProperty("PhiMax",        m_phiMax        =  180.    );
    declareProperty("Ki2Cut1",       m_ki2Cut1       = 100.     );
    declareProperty("MaxCall1",      m_maxCall1      = 500      );
    declareProperty("Ki2Cut2",       m_ki2Cut2       = 99999.   );
    declareProperty("MaxCall2",      m_maxCall2      = 500      );
    declareProperty("Dx_max",        m_dx_max        = 4.       );
    declareProperty("Beta0",         m_beta0         = 0.71     );

    m_radToDeg = 180. / 3.14159;

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

    // Define Minuit
    m_minuit = new TMinuit(4);
    
    // Set the function
    m_minuit->SetFCN(fcn);

    m_compteurTest=0;

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

    m_detSvc=0;
    sc = service("GlastDetSvc", m_detSvc);
    
    //-----------------------------------
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
    
    m_detSvc=0;
    sc = service("GlastDetSvc", m_detSvc);
    
    if ( sc.isFailure() ) {
      log << MSG::ERROR << "StdMipFindingTool failed to get GlastDetSvc" << endreq;
      return sc;
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
    
    //  for (int itower=0;itower=0;itower<=15)
    //    for (int ilayer=0;ilayer=0;ilayer<=15)
    //      for (int icolumn=0;icolumn=0;icolumn<=15){
    int itower  = 0;
    int ilayer  = 0;
    int icolumn = 0;
    
    // create Volume Identifier for segment 0 of this crystal
    idents::VolumeIdentifier segm0Id;
    segm0Id.append(m_eLATTowers);
    segm0Id.append(itower/m_xNum);
    segm0Id.append(itower%m_xNum);
    segm0Id.append(m_eTowerCAL);
    segm0Id.append(ilayer);
    segm0Id.append(ilayer%2); 
    segm0Id.append(icolumn);
    segm0Id.append(m_eXtal);
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
    for(int ifield = 0; ifield<fSegment; ifield++)segm11Id.append(segm0Id[ifield]);
    segm11Id.append(m_nCsISeg-1); // set segment number for the last segment
    //get 3D transformation for the last segment of this crystal
    m_detSvc->getTransform3DByID(segm11Id,&transf);
    //get position of the center of the last segment
    Vector vect11 = transf.getTranslation();
    //      
    Point p0(0.,0.,0.);         
    // position of the crystal center
    Point pCenter = p0+(vect0+vect11)*0.5; 

    // populate crystal
    m_coordXtalZ0=pCenter.z();

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
	            segm0Id.append(m_eLATTowers);
	            segm0Id.append(itower/m_xNum);
	            segm0Id.append(itower%m_xNum);
	            segm0Id.append(m_eTowerCAL);
	            segm0Id.append(ilayer);
	            segm0Id.append(ilayer%2); 
	            segm0Id.append(icolumn);
	            segm0Id.append(m_eXtal);
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
	            //      
	            Point p0(0.,0.,0.);	      
	            // position of the crystal center
	            Point pCenter = p0+(vect0+vect11)*0.5; 

	            double delta_x=(ilayer%2==0) ? 0.5*m_CsILength : 0.5*m_CsIWidth;
	            double delta_y=(ilayer%2==1) ? 0.5*m_CsILength : 0.5*m_CsIWidth;

	            m_XtalXlo[itower][ilayer][icolumn]=pCenter.getX()-delta_x;
	            m_XtalXhi[itower][ilayer][icolumn]=pCenter.getX()+delta_x;
	            m_XtalYlo[itower][ilayer][icolumn]=pCenter.getY()-delta_y;
	            m_XtalYhi[itower][ilayer][icolumn]=pCenter.getY()+delta_y;
	            m_XtalZlo[itower][ilayer][icolumn]=pCenter.getZ()-0.5*m_CsIHeight;
	            m_XtalZhi[itower][ilayer][icolumn]=pCenter.getZ()+0.5*m_CsIHeight;
	
	            if (pCenter.getX()<m_xMinCalEdge) m_xMinCalEdge=pCenter.getX();
	            if (pCenter.getX()>m_xMaxCalEdge) m_xMaxCalEdge=pCenter.getX();
	            if (pCenter.getY()<m_yMinCalEdge) m_yMinCalEdge=pCenter.getY();
	            if (pCenter.getY()>m_yMaxCalEdge) m_yMaxCalEdge=pCenter.getY();
	            if (pCenter.getZ()<m_zMinCalEdge) m_zMinCalEdge=pCenter.getZ();
	            if (pCenter.getZ()>m_zMaxCalEdge) m_zMaxCalEdge=pCenter.getZ();
            }
        }
    }
  
    m_xMinCalEdge-=m_CsIWidth/2;
    m_yMinCalEdge-=m_CsIWidth/2;
    m_zMinCalEdge-=m_CsIHeight/2;
    m_xMaxCalEdge+=m_CsIWidth/2;
    m_yMaxCalEdge+=m_CsIWidth/2;
    m_zMaxCalEdge+=m_CsIHeight/2;

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

    if (m_dir.getZ()!=0)
    {
        // zlo plane
        lambda=(zlow-m_entryP.getZ())/m_dir.getZ();
        m_C=m_entryP+lambda*m_dir;

        if (m_C.getY()>ylow && m_C.getY()<yhig && m_C.getX()>xlow && m_C.getX()<xhig)
	    if (!l1)
	    {
	        m_E1=m_C;
	        l1=true;
	    }
	    else
	    {
	        m_E2=m_C;
	        l2=true;
	    }

        // zhi plane
        lambda=(zhig-m_entryP.getZ())/m_dir.getZ();
        m_C=m_entryP+lambda*m_dir;

        if (m_C.getY()>ylow && m_C.getY()<yhig && m_C.getX()>xlow && m_C.getX()<xhig)
	    if (!l1)
	    {
	        m_E1=m_C;
	        l1=true;
	    }
	    else
	    {
	        m_E2=m_C;
	        l2=true;
	    }
    }

    if (m_dir.getX()!=0 && !(l1 && l2))
    {
        //xlow plane
        lambda=(xlow-m_entryP.getX())/m_dir.getX();
        m_C=m_entryP+lambda*m_dir;
        if (m_C.getY()>ylow && m_C.getY()<yhig && m_C.getZ()>zlow && m_C.getZ()<zhig)
	    if (!l1)
	    {
	        m_E1=m_C;
	        l1=true;
	    }
	    else
	    {
	        m_E2=m_C;
	        l2=true;
	    }

        //xhig plane
        lambda=(xhig-m_entryP.getX())/m_dir.getX();
        m_C=m_entryP+lambda*m_dir;

        if (m_C.getY()>ylow && m_C.getY()<yhig && m_C.getZ()>zlow && m_C.getZ()<zhig)
	    if (!l1)
	    {
	        m_E1=m_C;
	        l1=true;
	    }
	    else
	    {
	        m_E2=m_C;
	        l2=true;
	    }
    }

    if (m_dir.getY()!=0 && !(l1 && l2))
    {
        //ylow plane
        lambda=(ylow-m_entryP.getY())/m_dir.getY();
        m_C=m_entryP+lambda*m_dir;

        if (m_C.getX()>xlow && m_C.getX()<xhig && m_C.getZ()>zlow && m_C.getZ()<zhig)
	    if (!l1)
	    {
	        m_E1=m_C;
	        l1=true;
	    }
	    else
	    {
	        m_E2=m_C;
	        l2=true;
	    }
    }

    // yhig plane
    lambda=(yhig-m_entryP.getY())/m_dir.getY();
    m_C=m_entryP+lambda*m_dir;

    if (m_C.getX()>xlow && m_C.getX()<xhig && m_C.getZ()>zlow && m_C.getZ()<zhig)
    if (!l1)
    {
	    m_E1=m_C;
	    l1=true;
    }
    else
    {
	    m_E2=m_C;
	    l2=true;
    }

    double ec=0.0001;

    if (l1 && l2)
    {
        m_vect=m_E2-m_E1;
        ec=sqrt(m_vect*m_vect)/m_CsIHeight;
    }

    return ec;
}

//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipXtals(const Event::CalXtalRecCol* calXtalRecCol, const Event::CalCluster* calCluster)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findMipXtals in StdMipFindingTool" << endreq;

    double emip1       = 6.0;
    double emip2       = 50.0;
    int    numMipXtals = 0;

    m_simpCluCent=calCluster->getPosition();
    //Loop over crystals in the collection

    idents::CalXtalId      xTalId;

    for(Event::CalXtalRecCol::const_iterator xTalIter=calXtalRecCol->begin(); xTalIter != calXtalRecCol->end(); xTalIter++)
    {
        Event::CalXtalRecData* xTalData = *xTalIter;
        double                 xTalE    =  xTalData->getEnergy();
        if ( xTalE>emip1 && xTalE<emip2)
        {
            m_calMipXtalVec.push_back(Event::CalMipXtal());
            numMipXtals++;

            xTalId = xTalData->getPackedId();

            Event::CalMipXtal& mipXtal = m_calMipXtalVec.back();
            //Update the values for this MipXtal
            mipXtal.setFree(true);
            mipXtal.setXtal(xTalData);
            //Update the distance between hit and single cluster centroid
            Point hitPos = xTalData->getPosition();
            Vector vHitClu = m_simpCluCent - hitPos;
            double d2C=sqrt(vHitClu*vHitClu);
            mipXtal.setD2C(d2C);

            //calXtal=0;
            //delete calXtal;
        }
    }
    log << MSG::DEBUG << " SG : findMipXtals - End" << endreq;
    return numMipXtals;
}
//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readMipXtals(Event::CalMipXtalVec calMipXtalVec)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readMipXtals in StdMipFindingTool" << endreq;
    int comptor=0;

    //Loop over crystals in the calMipXtalVec
    for(Event::CalMipXtalVec::iterator xTalIter=calMipXtalVec.begin(); xTalIter != calMipXtalVec.end(); xTalIter++)
    {
        Event::CalMipXtal calMipXtal=*xTalIter;

        log << MSG::DEBUG << "readMipXtals - MipXtal No" << comptor << "-------" << endreq;
        /*
        log << MSG::DEBUG << "readMipXtals - Free=" << calMipXtal.getFree()  << endreq;
        log << MSG::DEBUG << "readMipXtals - D2C =" << calMipXtal.getD2C()   << endreq;
        log << MSG::DEBUG << "readMipXtals - Ener=" << calMipXtal.getXtal().getEnergy()   << endreq;
        log << MSG::DEBUG << "readMipXtals - PosX=" << calMipXtal.getXtal().getPosition().x()   << endreq;
        log << MSG::DEBUG << "readMipXtals - PosY=" << calMipXtal.getXtal().getPosition().y()   << endreq;
        log << MSG::DEBUG << "readMipXtals - PosZ=" << calMipXtal.getXtal().getPosition().z()   << endreq;
        */
        //calMipXtal.print();
        comptor++;
    }
        
    log << MSG::DEBUG << " SG : readMipXtals End" << endreq;
}
//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::readMipTracks(Event::CalMipTrackVec calMipTrackVec)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : readMipTracks in StdMipFindingTool" << endreq;
    int comptor=0;
    //Loop over tracks in the calMipTrackVec
    for(Event::CalMipTrackVec::iterator xTalIter=calMipTrackVec.begin(); xTalIter != calMipTrackVec.end(); xTalIter++)
    {
        Event::CalMipTrack calMipTrack=*xTalIter;

        log << MSG::DEBUG << "readMipTrack - MipTrack No" << comptor << "-------" << endreq;
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

        comptor++;
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
    
    // Why isn't this all we need to do?
    // Since this a vector of objects (as is the inherited vector), the clear should clean everything up
    // properly... 
    calMipTrackVec->clear();

    log << MSG::DEBUG << " SG : clearCalMipXtalVec End" << endreq;
}
//-----------------------------------------------------------------------------------------------------------------
int StdMipFindingTool::findMipTracks()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG << " SG : findMipTracks in StdMipFindingTool" << endreq;

    m_hid0    = -1;
    m_d2Cmax  = -50;
    m_xE      = 0.5;
    m_ki2cut = 10000;
    m_pi      = TMath::Pi();

    sc=readGlastDet();
    ASYMER2  = pow(20.,2);
    TRANSER2 = pow(m_CsIWidth,2);

    m_nbTracks=-1;

    Event::CalMipXtal calMipXtal;

    while (findC0())
    {
        m_first=true;
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
      
        if (!findCn())
	    {
	        m_calMipXtalVec[m_hid0].setFree(true);
	        m_calMipXtalVec[m_hid1].setFree(true);
	        m_calMipXtalVec[m_hid2].setFree(true);
	        m_calMipTrackVec.pop_back();
	        m_nbTracks--;
	        continue;
	    }

        fitTrack(m_ki2Cut1, m_maxCall1);

        //sylvain change
        //m_goodfit=true;
        //

        if (!m_goodfit)
	    {
	        m_calMipXtalVec[m_hid0].setFree(true);
	        m_calMipXtalVec[m_hid1].setFree(true);
	        m_calMipXtalVec[m_hid2].setFree(true);
	        m_calMipXtalVec[m_hidn].setFree(true);
	        m_calMipTrackVec.pop_back();
	        m_nbTracks--;
	        continue;
	    }

        while (findCn()){}

        //sylvain commented
        //fitTrack(m_ki2Cut2, m_maxCall2);

        //sylvain change test
        if (m_calMipTrackVec.back().getNh()<4 && m_calMipTrackVec.back().getKi2()>=100)
	    //if (m_calMipTrackVec.back().getNh()<2 && m_calMipTrackVec.back().getKi2()>=1000000)
	    {
	        for (int ih=0; ih<m_calMipTrackVec.back().getNh(); ih++)
	        {
	            m_calMipXtalVec[ih].setFree(true);
	        }
	        
            m_calMipTrackVec.pop_back();
	        m_nbTracks--;
	    }

        if (m_nbTracks+1>0)
	    {
	        for (int itr=0; itr<=m_nbTracks; itr++)
	        {
	            m_kidof=-5;

	            Event::CalMipTrack calMipTrack=m_calMipTrackVec[itr];

	            //if (m_calMipTrackVec[itr].getNdof()>0)
	            if (calMipTrack.getNdof()>0)
		        {
		            //m_kidof=m_calMipTrackVec[itr].getKi2()/m_calMipTrackVec[itr].getNdof();
		            m_kidof=calMipTrack.getKi2()/calMipTrack.getNdof();
		        }
	            
	            m_dir  = calMipTrack.getDir();
	            m_C    = calMipTrack.getPoint();
	            m_vect = m_C-m_centroid;
	            m_vv   = m_vect.cross(m_dir);
	            m_d2C  = sqrt(m_vv*m_vv);
	            //m_calMipTrackVec[itr].setD2C(m_d2C);
	            calMipTrack.setD2C(m_d2C);

	            calMipTrack.setD2Edge(999999999);
	            calMipTrack.setLength(-999999999);

	            for (int ih=0; ih<m_calMipTrackVec[itr].getNh(); ih++)
		        {
		            Event::CalMipXtal* calMipXtal_ih=new Event::CalMipXtal();
		            *calMipXtal_ih = calMipTrack.at(ih);
		            double dd=D2Edge(calMipXtal_ih);

		            if (dd<calMipTrack.getD2Edge())
		            {
		                calMipTrack.setD2Edge(dd);
		            }
		  
                    for (int ihh=0; ihh<calMipTrack.getNh() &&ihh!=ih; ihh++)
		            {
		                Event::CalMipXtal* calMipXtal_ihh=new Event::CalMipXtal();
		                *calMipXtal_ihh=calMipTrack.at(ihh);
		                m_vv=calMipXtal_ih->getXtal()->getPosition()-calMipXtal_ihh->getXtal()->getPosition();

		                double d=sqrt(m_vv*m_vv);
		                if (d > m_calMipTrackVec.back().getLength())
			            {
			                m_calMipTrackVec.back().setLength(d);
			            }
		                calMipXtal_ihh=0;
		                delete calMipXtal_ihh;
		            }
		
		            calMipXtal_ih=0;
		            delete calMipXtal_ih;
		        }
	        }
	    }
    }
    
    log << MSG::DEBUG << " SG : findMipTracks - End" << endreq;
    return m_nbTracks = m_calMipTrackVec.size();
}
//-----------------------------------------------------------------------------------------------------------------
  bool StdMipFindingTool::findC0()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    log << MSG::DEBUG << " SG : findC0 in StdMipFindingTool" << endreq;

    m_hid           = -1;
    m_hid0          = -1;
    m_d2C           = -40;
    m_d2Cmax        = -50;
    m_dmin=500000;

    // I think this loop is looking for the xtal with the largest d2c... 
    for(Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;

        // This test will very nearly always be true because m_d2Cmax will always be -40?
        // So change this to set m_d2Cmax to the value for the selected crystal
        if (calMipXtal.getFree() && calMipXtal.getFreeC0() && calMipXtal.getD2C() > m_d2Cmax)
	    {
	        //m_d2Cmax = m_d2C;
	        m_d2Cmax = calMipXtal.getD2C();
	        m_hid0   = m_hid;
	    }
    }
    
    if (m_hid0 >= 0)
    {
        Event::CalMipXtal& calMipXtalRef = m_calMipXtalVec[m_hid0];
        calMipXtalRef.setFree(false);
        calMipXtalRef.setFreeC0(false);
        log << MSG::DEBUG << " SG : findC0 End - True" << endreq;
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
    m_dmin = 8*m_CsIWidth;

    int itow0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getTower();
    int ilay0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getLayer();
    int icol0 = m_calMipXtalVec[m_hid0].getXtal()->getPackedId().getColumn();


    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;
        m_dE=fabs(1-calMipXtal.getXtal()->getEnergy()/m_calMipXtalVec[m_hid0].getXtal()->getEnergy());
        if (calMipXtal.getFree())
	    {
	        m_vv     = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[m_hid0].getXtal()->getPosition();
	        m_u01    = m_vv.unit();
	        m_entryP = m_calMipXtalVec[m_hid0].getXtal()->getPosition();
	        m_dir    = m_u01;

	        int itow1  = m_calMipXtalVec[m_hid].getXtal()->getPackedId().getTower();
	        int ilay1  = m_calMipXtalVec[m_hid].getXtal()->getPackedId().getLayer();
	        int icol1  = m_calMipXtalVec[m_hid].getXtal()->getPackedId().getColumn();
	  
            //SG change
	        if (itow1==itow0)
	        {
	            double ec0 = Ecor(m_calMipXtalVec[m_hid0]);
	            double dx0 = fabs(m_calMipXtalVec[m_hid0].getXtal()->getEnergy()/ec0-13.2);
	      
                double ec1 = Ecor(m_calMipXtalVec[m_hid]);
	            double dx1 = fabs(m_calMipXtalVec[m_hid].getXtal()->getEnergy()/ec1-13.2);
	      
                //SG change
	            if (dx0<m_dx_max && ec0>0.001 && dx1<m_dx_max && ec1>0.001)
		        {
		            double d01 = sqrt(m_vv*m_vv);
		            m_vect = m_centroid-m_calMipXtalVec[m_hid0].getXtal()->getPosition();
		            double beta = m_vv*m_vect;
		            if (d01>0)
		            {
		                beta/=d01*sqrt(m_vect*m_vect);
		            }
		            else
		            {
		                beta=0;
		            }
		            
                    if (d01<m_dmin && beta>m_beta0)
		            {
		                m_dmin = d01;
		                m_hid1 = m_hid;
		            }
		        }
	        }
	    }
    }

    if (m_hid1>=0)
    {
        Event::CalMipXtal& calMipXtalRef = m_calMipXtalVec[m_hid1];
        calMipXtalRef.setFree(false);
        m_D01point=m_calMipXtalVec[m_hid0].getXtal()->getPosition();
        m_D01vector=m_u01;
        log << MSG::DEBUG << " SG : findC1 End - True" << endreq;
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
    log << MSG::DEBUG << " SG : findC2 in StdMipFindingTool" << endreq;

    m_hid  = -1;
    m_hid2 = -1;
    m_dmin = 9999999;

    m_entryP = m_D01point;
    m_dir    = m_D01vector;

    int itow1  = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getTower();
    int ilay1  = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getLayer();
    int icol1  = m_calMipXtalVec[m_hid1].getXtal()->getPackedId().getColumn();


    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;

        if (calMipXtal.getFree())
	    {
	        int itow2=calMipXtal.getXtal()->getPackedId().getTower();
	        int ilay2=calMipXtal.getXtal()->getPackedId().getLayer();
	        int icol2=calMipXtal.getXtal()->getPackedId().getColumn();
	        if (itow2==itow1 && abs(ilay2-ilay1)<3)
	        {
	            double ec=Ecor(calMipXtal);
	            double dx=fabs(calMipXtal.getXtal()->getEnergy()/ec-13.2)/6.;
	            
                if (dx<m_dx_max && ec>0.001)
		        {
		            m_u02 = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[m_hid0].getXtal()->getPosition();
		            m_u12 = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[m_hid1].getXtal()->getPosition();
		            m_vv  = m_u01.cross(m_u02);
		            double d2D01=sqrt(m_vv*m_vv);
		            double d12=sqrt(m_u12*m_u12);
		            //SG change
		            if (d2D01<m_dmin && d12<16*m_CsIWidth)
		            //		  if (d2D01<m_dmin && d12<8*m_CsIWidth)
		            {
		                m_dmin=d2D01;
		                m_hid2=m_hid;
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

        m_hidn=m_hid2;

        if (m_first)
	    {
	        m_nbTracks++;
	        m_calMipTrackVec.push_back(Event::CalMipTrack());
	        Event::CalMipTrack& calMipTrackRef =m_calMipTrackVec.back();
	        calMipTrackRef.setNdof(-1);
	        calMipTrackRef.setKi2(-1);
	    }
        Event::CalMipTrack& calMipTrack =m_calMipTrackVec.back();

        Event::CalMipXtal& calMipXtalRef0 = m_calMipXtalVec[m_hid0];
        calMipTrack.push_back(calMipXtalRef0);

        Event::CalMipXtal& calMipXtalRef1 = m_calMipXtalVec[m_hid1];
        calMipTrack.push_back(calMipXtalRef1);

        Event::CalMipXtal& calMipXtalRef2 = m_calMipXtalVec[m_hid2];
        calMipTrack.push_back(calMipXtalRef2);

        calMipTrack.setPoint(m_D01point);
        calMipTrack.setDir(m_D01vector);

        Vector dir = calMipTrack.getDir();
        Point  C   = calMipTrack.getPoint();

        if (dir.z()!=0)
	    {
	        m_entryPz=m_coordXtalZ0;
	        log << MSG::DEBUG <<"m_entryPz="<<m_entryPz<<" - doit etre -58,... val des readGlastDet"<<endreq;
	        m_d=(m_entryPz-C.z())/dir.z();
	        m_entryP.setX(C.x()+m_d*dir.x());
	        m_entryP.setY(C.y()+m_d*dir.y());
	        m_entryP.setZ(C.z()+m_d*dir.z());
	    }
        else
	    {
	        m_entryPz=C.z();
	        m_entryP=C;
	    }

        calMipTrack.setPoint(m_entryP);

        log << MSG::DEBUG << " SG : findC2 End - True" << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " SG : findC2 End - False" << endreq;
        return false;
    }
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::fitTrack(double ki2Cut, int maxCall)
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " SG : fitTrack in StdMipFindingTool" << endreq;

    m_fcnNb = 0;
    m_ki2   = 0;
  
    // Initial direction
    Event::CalMipTrack& calMipTrack=m_calMipTrackVec.back();
    m_dir      = calMipTrack.getDir();
    m_theta    = TMath::ACos(m_dir.z());
    m_thetaErr = 1.;

    if (m_dir.x() != 0) m_phi  = atan(m_dir.y()/m_dir.x());
    if (m_dir.x() < 0)  m_phi += m_pi;
    if (m_phi > m_pi)   m_phi -= 2 * m_pi;
  
    m_phiErr = 1.;
    m_entryP = calMipTrack.getPoint();
  
    // Beginning of m_minuit
    m_minuit->mninit(5,5,7);
  
    double arglist[10];
    int    ierflg = 0;
  
    arglist[0] = -1;
    m_minuit->mnexcm("SET PRI", arglist, 1, ierflg);
    m_minuit->mnexcm("SET NOW", arglist, 1, ierflg);
    arglist[0] = 1;    // up=1 car chi2
    m_minuit->mnexcm("SET ERR", arglist, 1, ierflg);
    arglist[0] = 2;
    m_minuit->mnexcm("SET STR", arglist, 1, ierflg);
  
    double vstart[4];

    vstart[0] = m_entryP.x();
    vstart[1] = m_entryP.y();
    vstart[2] = m_theta;
    vstart[3] = m_phi;

    double step[4] = {2,2,0.035,0.035};
  
    m_minuit->mnparm(0, "Xentry", vstart[0], step[0],-999999.,999999.,ierflg);
    m_minuit->mnparm(1, "Yentry", vstart[1], step[1],-999999.,999999.,ierflg);
    m_minuit->mnparm(2, "Theta",  vstart[2], step[2],m_thetaMin/m_radToDeg,m_thetaMax/m_radToDeg,ierflg);
    m_minuit->mnparm(3, "Phi",    vstart[3], step[3],m_phiMin/m_radToDeg,m_phiMax/m_radToDeg,ierflg);
  
    arglist[0] = 666;
  
    m_minuit->mnexcm("CALL FCN", arglist, 1, ierflg);
  
    // minimization
    arglist[0] = maxCall; // maxcalls
    arglist[1] = 1;   // tolerance
    m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
    m_minuit->GetParameter(0, m_Xentry, m_XentryErr);
    m_minuit->GetParameter(1, m_Yentry, m_YentryErr);
    m_minuit->GetParameter(2, m_theta , m_thetaErr );
    m_minuit->GetParameter(3, m_phi   , m_phiErr   );

    log << MSG::DEBUG <<"m_fcnNb="<<m_fcnNb<<endreq;
  
    // Retrieve the chi-square and determine number of degrees of freedom
    m_chiSquare = 0.;
    m_nPar      = 4;
    m_ndof      = -10;
    m_ndof      = int(m_ki2p) - m_minuit->GetNumFreePars();
    m_kidof     = -5.;

    log << MSG::DEBUG <<"m_ndof ="<<m_ndof<<endreq;
    if (m_ndof>0)
    {
        m_Ki2Tot++;
        m_kidof = m_ki2 / m_ndof;
    }
  
    m_goodfit=false;

    log << MSG::DEBUG <<"m_kidof="<<m_kidof<<" - m_ki2cut="<<m_ki2cut<<endreq;

    if (m_kidof>0 && m_kidof<ki2Cut)
    {
        m_goodfit = true;
        calMipTrack.setKi2(m_ki2);
        calMipTrack.setNdof(m_ndof);
        m_entryP.setX(m_Xentry);
        m_entryP.setY(m_Yentry);
        calMipTrack.setPoint(m_entryP);
        m_dir.setX(sin(m_theta)*cos(m_phi));
        m_dir.setY(sin(m_theta)*sin(m_phi));
        m_dir.setZ(cos(m_theta));
        calMipTrack.setDir(m_dir);
        m_thetaf     = m_theta*m_radToDeg;
        m_thetafErr  = m_thetaErr*m_radToDeg;
        m_phif       = m_phi*m_radToDeg;
        m_phifErr    = m_phiErr*m_radToDeg;
        m_thetaf     = 180-m_thetaf;
        m_phif      += 180;
        m_Xentryf    = m_Xentry;
        m_Yentryf    = m_Yentry;
        m_kidoff     = m_kidof;
    }
  
    m_theta     = m_theta*m_radToDeg;
    m_thetaErr *= m_radToDeg;
    m_phi      *= m_radToDeg;
    m_phiErr   *= m_radToDeg;
    m_theta     = 180-m_theta;
    m_phi      += 180;
    if (m_phi>180) m_phi-=360;
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
    double d2Edge = 99999999999;
  
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
bool StdMipFindingTool::findCn()
{
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
    log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
    log << MSG::DEBUG << " SG : findCn in StdMipFindingTool" << endreq;

    int hidlast = m_hidn;
    int itowl   = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getTower();
    int ilayl   = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getLayer();
    int icoll   = m_calMipXtalVec[m_hidn].getXtal()->getPackedId().getColumn();

    m_hidn = -1;
    m_hid  = -1;
    m_dmin = 9999999;

    Event::CalMipTrack& calMipTrack =m_calMipTrackVec.back();

    m_entryP = m_C=calMipTrack.getPoint();
    m_dir    = calMipTrack.getDir();

    for (Event::CalMipXtalVec::iterator xTalIter=m_calMipXtalVec.begin(); xTalIter != m_calMipXtalVec.end(); xTalIter++)
    {
        m_hid++;
        Event::CalMipXtal calMipXtal = *xTalIter;
      
        if (calMipXtal.getFree())
	    {
	        int itown = calMipXtal.getXtal()->getPackedId().getTower();
	        int ilayn = calMipXtal.getXtal()->getPackedId().getLayer();
	        int icoln = calMipXtal.getXtal()->getPackedId().getColumn();

	        //SG change
	        if (itown==itowl)
	        //	  if ((itown==itowl && abs(ilayn-ilayl)<3) || itown!=itowl)
	        {
	            double ec = Ecor(calMipXtal);
	            double dx = fabs(calMipXtal.getXtal()->getEnergy()/ec-13.2)/6.;
	            if (dx<3*m_dx_max && ec>0.001)
		        {
		            m_vect = calMipXtal.getXtal()->getPosition()-m_C;
		            m_vv   = m_dir.cross(m_vect);
		            double d2dir = sqrt(m_vect*m_vect);
		            m_vect = calMipXtal.getXtal()->getPosition()-m_calMipXtalVec[hidlast].getXtal()->getPosition();
		            double d2last = sqrt(m_vect*m_vect);
		            //SG change
		            if (d2dir<m_dmin && d2last<16*m_CsIWidth)
		            //		  if (d2dir<m_dmin && d2last<8*m_CsIWidth)
		            {
		                m_dmin = d2dir;
		                m_hidn = m_hid;
		            }
		        }
	        }
	    }
    }

    if (m_hidn>=0)
    {
        Event::CalMipXtal& calMipXtalRef_m_hidn = m_calMipXtalVec[m_hidn];
        calMipXtalRef_m_hidn.setFree(false);
        calMipTrack.push_back(calMipXtalRef_m_hidn);
        log << MSG::DEBUG << " SG : findCn End - True" << endreq;
        return true;
    }
    else
    {
        log << MSG::DEBUG << " SG : findCn End - False" << endreq;
        return false;
    }
}

//-----------------------------------------------------------------------------------------------------------------
// Main method
StatusCode StdMipFindingTool::findMIPCandidates()
{
    // Presume success unless something bad happens
    StatusCode sc = StatusCode::SUCCESS;
    MsgStream log(msgSvc(), name());
        
    log << MSG::DEBUG << " SG : Enter findMIPCandidates" << endreq;
        
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
    int numMipXtal = findMipXtals(calXtalRecCol, calCluster);
    if (numMipXtal<1) return  StatusCode::FAILURE;
    log << MSG::DEBUG << " SG : numMipXtals=" << numMipXtal  << endreq;
        
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
    readMipXtals(m_calMipXtalVec);
        
    //
    // Task 8: CalMipTrackVec Cleaner
    //
    clearCalMipTrackVec(&m_calMipTrackVec);
    readMipTracks(m_calMipTrackVec);

    //
    // Task 9: CalMipTrack &  CalMIPsCol Reader
    //
     readCalMipTrackCol();

    return sc;
}

//-----------------------------------------------------------------------------------------------------------------
void StdMipFindingTool::fcn (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    StatusCode sc = StatusCode::SUCCESS;
    //MsgStream log(msgSvc(), name());
    //if (m_fcnNb==0)
    //{
        //log << MSG::DEBUG <<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endreq;
        //log << MSG::DEBUG << "SG : fcn in StdMipFindingTool"<<endreq;
    //}
    m_fcnNb++;

    //Calculate chisquare
    m_ki2p = 0;
    m_ki2  = 0.;

    double dx2,err2;

    m_entryP.setX(par[0]);
    m_entryP.setY(par[1]);
    m_entryP.setZ(m_entryPz);

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
    double                 Emean,x,y,z;

    Event::CalMipTrack& calMipTrack =m_calMipTrackVec.back();

    for(int k=0;k<=m_ki2Type;k++)
    { 
        for (int ih=0; ih<calMipTrack.getNh(); ih++)
	    {
	        calMipXtal = calMipTrack.at(ih);
	  
	        int itow = calMipXtal.getXtal()->getPackedId().getTower();
	        int ilay = calMipXtal.getXtal()->getPackedId().getLayer();
	        int icol = calMipXtal.getXtal()->getPackedId().getColumn();

	        if (iflag==666)
	        {
	            xtalData = calMipXtal.getXtal();
	            Emean    = xtalData->getEnergy();
	            point    = xtalData->getPosition();
	            x        = point.x();
	            y        = point.y();
	            z        = point.z();
	        }

	        if (k==0)
	        {
	            xtalData = calMipXtal.getXtal();
	            point    = xtalData->getPosition();
	            m_vect   = point-m_entryP;
	            m_vv     = m_dir.cross(m_vect);
	            dx2      = m_vv*m_vv;
	            err2     = 100.;
	        }
	        else if (k==1)
	        {
	            double ec = Ecor(calMipXtal);
	            dx2  = pow(calMipXtal.getXtal()->getEnergy()/ec-13.2,2);
	            err2 = 36.*pow(ct,2);
	        }
	  
            m_ki2 += dx2/err2;
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

    int compteur=0;

    if (p_calMipTrackCol!=0)
    {
        for(Event::CalMipTrackCol::const_iterator calMipTrackIter=p_calMipTrackCol->begin(); calMipTrackIter != p_calMipTrackCol->end(); calMipTrackIter++)
        {
            log << "------------------------------------------------------------" << endreq;
            log << "compteur col=" << compteur << endreq;
            log << "----------------" << endreq;
            compteur++;
            Event::CalMipTrack* p_calMipTrack    =  *calMipTrackIter;
            p_calMipTrack->writeOut(log);
            log << "------------------------------------------------------------" << endreq;
        }

        compteur=0;
    }
    log << MSG::DEBUG << " SG : readCalMipTrackCol End" << endreq;
}

