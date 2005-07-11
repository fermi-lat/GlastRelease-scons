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

#define FSDD_NSTEPS_MAX 1000
#define FSDD_NPOINTS_MAX 100
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

class FullShowerDevelopmentDescription{
 private:
  /// Detector Service
  IGlastDetSvc *m_detSvc;
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
  // Geometry
  double FSDD_towerPitch;
  double FSDD_calZTop;
  double FSDD_calZBot;
  double FSDD_crackHalfWidth;
  double FSDD_cellVertPitch;
  double FSDD_cellHorPitch;
  double FSDD_CsIHeight;
  double FSDD_CsIWidth;
  double FSDD_CsILength;
  // Radiation lengths
  double FSDD_MOLRAD;
  double FSDD_gCSI;
  double FSDD_gCRK;
  double FSDD_XCSI;
  double FSDD_XCRK;
  // Parameters used for the radial shower development description
  double FSDD_CORE0;
  double FSDD_CORE1;
  double FSDD_TAIL0;
  double FSDD_TAIL1;
  double FSDD_TAIL2;
  double FSDD_TAIL3;
  double FSDD_PCT0;
  double FSDD_PCT1;
  double FSDD_PCT2;

 private:
  double wideningfactor;

 private:
  double FSDD_CalSafeBoundaries[3][2];

 private:
  int FSDD_NCircle;
  double *FSDD_XCircle;
  double *FSDD_YCircle;
  double *FSDD_RCircle;

 public:
  FullShowerDevelopmentDescription(IGlastDetSvc *m_detSvc_input, int type_input, double zstep_input, double radialcontainedfraction_input);
  virtual ~FullShowerDevelopmentDescription();

 private:
  void Initialize();
  void Reset();
  void PrepXYPoints();
  double RCore(double x);
  double RTail(double x);
  double PCore(double x);
  double RadialProfile(double r, double x);
  double EffectiveRadius(double z, double radialcontainedfraction);
  void WhereInCal(double *xyz, int *whereincal);
  void RemoveEmptySteps();
  void GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend);

 public:
  bool Compute(double *pp, double *vv, double startx0_input, double x0maxshower_input);
  bool ConvertToFixedX0(double x0step, FullShowerDevelopmentDescription *shmm);
};

class FullShowerDevelopmentDescriptionManager{
 public:
  /// Detector Service
  IGlastDetSvc *m_detSvc;
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
};
