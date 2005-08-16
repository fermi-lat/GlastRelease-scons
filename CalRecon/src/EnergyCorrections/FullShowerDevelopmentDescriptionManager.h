/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "GaudiKernel/MsgStream.h" 
#include "GaudiKernel/GaudiException.h" 
#include "GlastSvc/GlastDetSvc/IGlastDetSvc.h"

#define FSGM_XY_MAX 149800
#define FSGM_NPOINTS_MAX 100
#define FSGM_RPROF_R_MAX 1800
#define FSGM_RPROF_T_MAX 250

#define FSDD_NSTEPS_MAX 1300
#define FSDD_NMAX 20

/**   
* @class FullShowerDevelopmentDescription
* @author Philippe Bruel
*
* Tool that describes the shower developement in the calorimeter given
* the length in X0 seen in the tracker and the position of the shower maximum
*
* $Header$
*/

class FullShowerGeometryManager{
 private:
  /// Detector Service
  IGlastDetSvc *m_detSvc;
 private:
  double geom_xy_step;
  bool geom_x_mat[FSGM_XY_MAX][2]; // [0] = false : void; [0] = true : CsI; [1] = true : crack 
  bool geom_y_mat[FSGM_XY_MAX][2];
 public:
  // Geometry variables
  double FSGM_towerPitch;
  double FSGM_calZTop;
  double FSGM_calZBot;
  double FSGM_crackHalfWidth;
  double FSGM_cellVertPitch;
  double FSGM_cellHorPitch;
  double FSGM_CsIHeight;
  double FSGM_CsIWidth;
  double FSGM_CsILength;
 public:
  double FSGM_CalSafeBoundaries[3][2];
 public:
  int FSGM_NCircle;
  double FSGM_XCircle[FSGM_NPOINTS_MAX];
  double FSGM_YCircle[FSGM_NPOINTS_MAX];
  double FSGM_RCircle[FSGM_NPOINTS_MAX];
 private:
  // Parameters used for the radial shower development description
  double FSGM_CORE0;
  double FSGM_CORE1;
  double FSGM_TAIL0;
  double FSGM_TAIL1;
  double FSGM_TAIL2;
  double FSGM_TAIL3;
  double FSGM_PCT0;
  double FSGM_PCT1;
  double FSGM_PCT2;
  double FSGM_tmin;
  double FSGM_tmax;
 private:
  double RadProf_t_max;
  double RadProf_r_max;
  double RadProf_t_step;
  double RadProf_r_step;
  double RadProf[FSGM_RPROF_T_MAX][FSGM_RPROF_R_MAX];

 public:
  FullShowerGeometryManager(IGlastDetSvc *m_detSvc_input);
  virtual ~FullShowerGeometryManager();
 private:
  void Initialize();
  void FillXYPoints();
  double RCore(double x);
  double RTail(double x);
  double PCore(double x);
  void FillRadialProfile();
  void WhereInCalForGeom(double *xyz, int *whereincal);
  double RadialProfile(double t, double r);
 public:
  void WhereInCal(double *xyz, int *whereincal);
  double GetEffectiveRadius(double z, double radialcontainedfraction);
  double GetRadialProfile(double t, double r);
};

class FullShowerDevelopmentDescription{
 public:
  int Type;
  int NStep;
  double ZStep;
  double RadialContainedFraction;
  double dX0[FSDD_NSTEPS_MAX];
  double X0[FSDD_NSTEPS_MAX];
  double RM[FSDD_NSTEPS_MAX];
  double materialfraction[3][FSDD_NSTEPS_MAX];
  double layerfraction[8][FSDD_NSTEPS_MAX];
  double x0maxshower;
  double crackmaxfrac;
  double x0crackmaxfrac;
  double startx0;
  double lastx0;
  double totx0cal;
  double totx0crack;
  double totx0lay[8];
  double posx0lay[8];

 private:
  // geometry
  FullShowerGeometryManager *m_fsgm;

 private:
  // Densities, radiation lengths and Moliere radius
  double FSDD_MOLRAD;
  double FSDD_gCSI;
  double FSDD_gCRK;
  double FSDD_XCSI;
  double FSDD_XCRK;

 private:
  double wideningfactor;

 private:


 public:
  FullShowerDevelopmentDescription(FullShowerGeometryManager *fsgm_input, int type_input, double zstep_input, double radialcontainedfraction_input);
  virtual ~FullShowerDevelopmentDescription();

 private:
  void Initialize();
  void Reset();
  void RemoveEmptySteps();
  void GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend);

 public:
  bool Compute(double *pp, double *vv, double startx0_input, double x0maxshower_input);
  bool ConvertToFixedX0(double x0step, FullShowerDevelopmentDescription *shmm);
  void SetWideningFactor(double widfact);
};

class FullShowerDevelopmentDescriptionManager{
 public:
  /// Detector Service
  IGlastDetSvc *m_detSvc;
  // geometry
  FullShowerGeometryManager *m_fsgm;
 public:
  int NDevelopment;
  double DXMax;
  double XMax[FSDD_NMAX];
  double ZStep;
  double RadialContainedFraction;
  double X0Step;
  double mintotx0cal;
  double maxtotx0cal;
  double meantotx0lay[8];
  double meanposx0lay[8];
  FullShowerDevelopmentDescription *FSDDMM[FSDD_NMAX];
  FullShowerDevelopmentDescription *FSDDX0[FSDD_NMAX];
  FullShowerDevelopmentDescription *CurrentFSDD;
 public:
  FullShowerDevelopmentDescriptionManager(IGlastDetSvc *m_detSvc_input, int nxmax, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input, double x0step);
  virtual ~FullShowerDevelopmentDescriptionManager();
  bool Compute(double *pp, double *vv, double startx0_input);
  void FillCurrentFSDD(double showerxmax);
  void SetWideningFactor(double widfact);
};
