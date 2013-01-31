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
#define FSGM_NPOINTS_MAX 1100
#define FSGM_NPOINTS_NRING 6
#define FSGM_NPOINTS_IN_RING 10
#define FSGM_RPROF_R_MAX 5000
#define FSGM_RPROF_T_MAX 750

#define FSDD_NSTEPS_MAX 1300
#define FSDD_NMAX 20

#define FSDD_XTAL_NMAX 200

/**   
* @class NewFullShowerDevelopmentDescription
* @author Philippe Bruel
*
* Tool that describes the shower developement in the calorimeter given
* the length in X0 seen in the tracker and the position of the shower maximum
*
* $Header$
*/

class NewFullShowerGeometryManager{
 private:
  /// Detector Service
  IGlastDetSvc *m_detSvc;
 private:
  double geom_xy_step;
  bool geom_x_mat[FSGM_XY_MAX][2]; // [0] = false : void; [0] = true : CsI; [1] = true : crack 
  bool geom_y_mat[FSGM_XY_MAX][2];
  int geom_x_wic[FSGM_XY_MAX][4];
  int geom_y_wic[FSGM_XY_MAX][4];
 public:
  // Geometry variables
  int FSGM_numX;
  int FSGM_numY;
  bool FSGM_flight_geom;
  double FSGM_towerPitch;
  double FSGM_calZTop;
  double FSGM_calZBot;
  double FSGM_crackHalfWidth;
  double FSGM_crackPara;
  double FSGM_crackPerp;
  double FSGM_crackXtal;
  double FSGM_cellVertPitch;
  double FSGM_cellHorPitch;
  double FSGM_CsIHeight;
  double FSGM_CsIWidth;
  double FSGM_CsILength;
 public:
  double FSGM_CalSafeBoundaries[3][2];
  double FSGM_FakeCsIThickness;
 public:
  int FSGM_NCircle;
  int FSGM_InCircle[FSGM_NPOINTS_MAX];
  double FSGM_XCircle[FSGM_NPOINTS_MAX];
  double FSGM_YCircle[FSGM_NPOINTS_MAX];
  double FSGM_RCircle[FSGM_NPOINTS_MAX];
  double FSGM_WCircle[FSGM_NPOINTS_MAX];
  double FSGM_LambdaCircle[FSGM_NPOINTS_MAX];
  double FSGM_PCircle[FSGM_NPOINTS_MAX];
  double FSGM_CPCircle[FSGM_NPOINTS_MAX];
  double FSGM_SPCircle[FSGM_NPOINTS_MAX];
  double FSGM_TCircle[FSGM_NPOINTS_MAX];
  double FSGM_TCircle2[FSGM_NPOINTS_MAX];
  double FSGM_TPedCircle[FSGM_NPOINTS_MAX];
  double FSGM_RadCircle[FSGM_NPOINTS_MAX];
  double FSGM_ECircle[FSGM_NPOINTS_MAX];
  double FSGM_ScaleCircle[FSGM_NPOINTS_MAX][FSDD_NMAX];
  double FSGM_ZCrackCircle[FSGM_NPOINTS_MAX];
  int opthe;
 private:
  /* // Parameters used for the radial shower development description */
  /* double FSGM_CORE0; */
  /* double FSGM_CORE1; */
  /* double FSGM_TAIL0; */
  /* double FSGM_TAIL1; */
  /* double FSGM_TAIL2; */
  /* double FSGM_TAIL3; */
  /* double FSGM_PCT0; */
  /* double FSGM_PCT1; */
  /* double FSGM_PCT2; */

  double FSGM_tmin;
  double FSGM_tmax;
 private:
  double RadProf_t_max;
  double RadProf_r_max;
  double RadProf_t_step;
  double RadProf_r_step;
  double RadProf[FSGM_RPROF_T_MAX][FSGM_RPROF_R_MAX];
  double RadProf2[FSGM_RPROF_T_MAX][FSGM_RPROF_R_MAX];
  //
  double FSGM_RCORE_0;
  double FSGM_RCORE_1;
  double FSGM_RTAIL_0;
  double FSGM_RTAIL_1;
  double FSGM_RTAIL_2;
  double FSGM_RTAIL_3;
  double FSGM_RTAIL_4;
  double FSGM_RTAIL_5;
  double FSGM_PCORE_0;
  double FSGM_PCORE_1;
  double FSGM_PCORE_2;
  double FSGM_PCORE_3;
  
  double FSGM_RCORE_0_0;
  double FSGM_RCORE_0_1;
  double FSGM_RCORE_1_0;
  double FSGM_RCORE_1_1;
  double FSGM_RTAIL_0_0;
  double FSGM_RTAIL_0_1;
  //  double FSGM_RTAIL_0_2;
  double FSGM_RTAIL_1_0;
  double FSGM_RTAIL_1_1;
  double FSGM_RTAIL_2_0;
  double FSGM_RTAIL_3_0;
  double FSGM_RTAIL_3_1;
  double FSGM_RTAIL_4_0;
  double FSGM_RTAIL_5_0;
  double FSGM_RTAIL_5_1;
  double FSGM_PCORE_0_0;
  double FSGM_PCORE_0_1;
  double FSGM_PCORE_0_2;
  double FSGM_PCORE_0_3;
  double FSGM_PCORE_1_0;
  double FSGM_PCORE_1_1;
  double FSGM_PCORE_1_2;
  double FSGM_PCORE_2_0;
  double FSGM_PCORE_2_1;
  double FSGM_PCORE_2_2;
  double FSGM_PCORE_2_3;
  double FSGM_PCORE_3_0;
  double FSGM_PCORE_3_1;
  double FSGM_PCORE_3_2;
  double FSGM_PCORE_3_3;
 public:
  NewFullShowerGeometryManager(IGlastDetSvc *m_detSvc_input);
  virtual ~NewFullShowerGeometryManager();
  void FillRadialProfile(double loge, int numtab);
  void SetOptHE(int opthein);
 private:
  void Initialize();
  void FillXYPoints();
  void oldFillXYPoints();
  void PrepRadialProfile(double loge);
  double RCore(double x);
  double RTail(double x);
  double PCore(double x);
  void FillRadialProfile(int numtab);
  void WhereInCalForGeom(double *xyz, int *whereincal);
  double RadialProfile(double t, double r);
  double RadialProfileMod(double t, double r, double pp, double rtailmod);
 public:
  void WhereInCal(double *xyz, int *whereincal);
  void WhereInCalLT(double *xyz, int *whereincal);
  void WhereInCalLT2(double *xyz, int *whereincal);
  void WhereInCalForGeomCU(double *xyz, int *whereincal);
  double GetEffectiveRadius(double z, double radialcontainedfraction);
  double GetRadialProfile(double t, double r);
  double GetRadialProfile2(double t, double r);
  double GetRadialProfileMod2(double t, double r, double pp,double rtailmod);
  double GetZCrack(double *xyz, double *pp, double *vv);
};

class NewMultiFullShowerDevelopmentDescription{
 public:
  int NDevelopment;
  double DXMax;
  double XMax[FSDD_NMAX];
  //
  int Type;
  int NStep;
  double ZStep;
  double ZStepRef;
  double RadialContainedFraction;
  double dX0[FSDD_NMAX][FSDD_NSTEPS_MAX];
  double X0[FSDD_NMAX][FSDD_NSTEPS_MAX];
  double RM[FSDD_NMAX][FSDD_NSTEPS_MAX];
  double materialfraction[FSDD_NMAX][4][FSDD_NSTEPS_MAX];
  double layerfraction[FSDD_NMAX][8][FSDD_NSTEPS_MAX];
  double xtalfraction[FSDD_NMAX][FSDD_XTAL_NMAX][FSDD_NSTEPS_MAX];
  double x0maxshower[FSDD_NMAX];
  double crackmaxfrac[FSDD_NMAX];
  double x0crackmaxfrac[FSDD_NMAX];
  double startx0[FSDD_NMAX];
  double lastx0[FSDD_NMAX];
  double totx0cal[FSDD_NMAX];
  double totx0crack[FSDD_NMAX];
  double totx0lay[FSDD_NMAX][8];
  double posx0lay[FSDD_NMAX][8];
  double crackextinction;

 private:
  // geometry
  NewFullShowerGeometryManager *m_fsgm;

 private:
  // Densities, radiation lengths and Moliere radius
  double FSDD_MOLRAD;
  double FSDD_gCSI;
  double FSDD_gCRK;
  double FSDD_XCSI;
  double FSDD_XCRK;
  int currentwhereincal[4];

 private:
  double wideningfactor;

 public:
  int NXtal;
  int OffSatu[16][8][12]; // -1 = on; >=0 : index of xtalfraction vector

 public:
  NewMultiFullShowerDevelopmentDescription(NewFullShowerGeometryManager *fsgm_input, int nxmax, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input);
  virtual ~NewMultiFullShowerDevelopmentDescription();

 private:
  void Initialize();
  void Reset();
  void GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend);

 public:
  double GoThroughCrack(double *pp, double *vv, double zstep, int nstepmax);
  double GetCrackAngle(double *pp, double *vv, double lambdain, double *pptraj);
  bool Compute(double *pp, double *vv, double startx0_input, double zstep_input);
  void SetWideningFactor(double widfact);
};

class NewFullShowerDevelopmentDescription{
 public:
  int Type;
  int NStep;
  double ZStep;
  double ZStepRef;
  double RadialContainedFraction;
  double dX0[FSDD_NSTEPS_MAX];
  double X0[FSDD_NSTEPS_MAX];
  double RM[FSDD_NSTEPS_MAX];
  double materialfraction[4][FSDD_NSTEPS_MAX];
  double layerfraction[8][FSDD_NSTEPS_MAX];
  double xtalfraction[FSDD_XTAL_NMAX][FSDD_NSTEPS_MAX];
  double x0maxshower;
  double crackmaxfrac;
  double x0crackmaxfrac;
  double startx0;
  double lastx0;
  double totx0cal;
  double totx0crack;
  double totx0lay[8];
  double posx0lay[8];
  double crackextinction;

 private:
  // geometry
  NewFullShowerGeometryManager *m_fsgm;

 private:
  // Densities, radiation lengths and Moliere radius
  double FSDD_MOLRAD;
  double FSDD_gCSI;
  double FSDD_gCRK;
  double FSDD_XCSI;
  double FSDD_XCRK;

 private:
  double wideningfactor;

 public:
  int NXtal;
  int OffSatu[16][8][12]; // -1 = on; >=0 : index of xtalfraction vector

 public:
  NewFullShowerDevelopmentDescription(NewFullShowerGeometryManager *fsgm_input, int type_input, double zstep_input, double radialcontainedfraction_input);
  virtual ~NewFullShowerDevelopmentDescription();

 private:
  void Initialize();
  void Reset();
  void GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend);

 public:
  bool Compute(double *pp, double *vv, double startx0_input, double x0maxshower_input, double zstep_input);
  bool ConvertToFixedX0(double x0step, NewFullShowerDevelopmentDescription *shmm);
  void SetWideningFactor(double widfact);
  void RemoveEmptySteps();
  void FillFromMultiFullShowerDevelopmentDescription(NewMultiFullShowerDevelopmentDescription *mfsddmm, int index);
};

class NewFullShowerDevelopmentDescriptionManager{
 public:
  /// Detector Service
  IGlastDetSvc *m_detSvc;
  // geometry
  NewFullShowerGeometryManager *m_fsgm;
 public:
  int NDevelopment;
  double DXMax;
  double XMax[FSDD_NMAX];
  double ZStep;
  double ZStepRef;
  double RadialContainedFraction;
  double X0Step;
  double mintotx0cal;
  double maxtotx0cal;
  double meantotx0lay[8];
  double meanposx0lay[8];
  int NXtal;
  int OffSatu[16][8][12]; // 0 = on; 1= off or saturated
  //
  NewMultiFullShowerDevelopmentDescription *MFSDDMM;
  NewFullShowerDevelopmentDescription *FSDDMM[FSDD_NMAX];
  NewFullShowerDevelopmentDescription *FSDDX0[FSDD_NMAX];
  NewFullShowerDevelopmentDescription *CurrentFSDD;
 public:
  NewFullShowerDevelopmentDescriptionManager(IGlastDetSvc *m_detSvc_input, int nxmax, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input, double x0step);
  virtual ~NewFullShowerDevelopmentDescriptionManager();
  bool Compute(double *pp, double *vv, double startx0_input, double zstep_input);
  void FillCurrentFSDD(double showerxmax);
  void SetWideningFactor(double widfact);
  void SetCrackExtinctionFactor(double crackext);
};
