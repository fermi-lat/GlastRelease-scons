/**   
*
* THIS CODE IS STILL NOT FULLY DOCUMENTED !!!!!!!!!!!!!!!!!!!!
* DESCRIPTION COMMENTS WILL BE ADDED SOON.
* Philippe Bruel
*
*/

#include "NewFullShowerDevelopmentDescriptionManager.h"

#include "idents/TowerId.h"
#include "idents/VolumeIdentifier.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"

#include "TMath.h"

NewFullShowerGeometryManager::NewFullShowerGeometryManager(IGlastDetSvc *m_detSvc_input)
  :m_detSvc(m_detSvc_input)
{
  Initialize();
}

NewFullShowerGeometryManager::~NewFullShowerGeometryManager()
{
  
}

void NewFullShowerGeometryManager::Initialize()
{
  // Retrieve geometry from Detector Service

  FSGM_CsIHeight = 19.9;
  m_detSvc->getNumericConstByName(std::string("CsIHeight"),&FSGM_CsIHeight);
  
  FSGM_CsIWidth = 26.7;
  m_detSvc->getNumericConstByName(std::string("CsIWidth"),&FSGM_CsIWidth);

  FSGM_CsILength = 326.0;
  m_detSvc->getNumericConstByName(std::string("CsILength"),&FSGM_CsILength);

  FSGM_cellVertPitch = 21.35;
  m_detSvc->getNumericConstByName(std::string("cellVertPitch"),&FSGM_cellVertPitch);

  FSGM_cellHorPitch = 27.84;
  m_detSvc->getNumericConstByName(std::string("cellHorPitch"),&FSGM_cellHorPitch);

  FSGM_towerPitch = 374.5;
  m_detSvc->getNumericConstByName(std::string("towerPitch"),&FSGM_towerPitch);

  // FSGM_crackHalfWidth is the average of X and Y layers crackHalfWidths
  FSGM_crackHalfWidth = ((FSGM_towerPitch-FSGM_CsILength)/2 + (FSGM_towerPitch-12*FSGM_cellHorPitch)/2)/2;
  FSGM_crackPara = (FSGM_towerPitch-FSGM_CsILength)/2.0;
  FSGM_crackPerp = (FSGM_towerPitch-12*FSGM_cellHorPitch)/2.0+(FSGM_cellHorPitch-FSGM_CsIWidth)/2.0;
  FSGM_crackXtal = (FSGM_cellHorPitch-FSGM_CsIWidth)/2.0;

  // Retrieving the top position of the calorimeter as in TkrUtil/*/src/TkrGeometrySvc.cxx
  FSGM_calZTop = -47.395;

  int m_Tower = 0;
  idents::VolumeIdentifier topLayerId;
  topLayerId.init(0,0);
  topLayerId.append(0);               // in Tower
  idents::TowerId t(m_Tower);  
  topLayerId.append(t.iy());          // yTower
  topLayerId.append(t.ix());          // xTower
  topLayerId.append(0);  // CAL
  topLayerId.append(0);  // layer
  topLayerId.append(0);  // x view
  StatusCode sc;
  HepGeom::Transform3D transfTop;
  int count;
  for (count=0;count<3;++count) {
    topLayerId.append(0);
    if((sc = m_detSvc->getTransform3DByID(topLayerId,&transfTop)).isSuccess()) break;
  }
  if(sc.isSuccess()) 
    {
      FSGM_calZTop = (transfTop.getTranslation()).z()+ FSGM_cellVertPitch/2;
    }

  FSGM_calZBot = FSGM_calZTop-8*FSGM_cellVertPitch;

  // prepare safe volume in order to define the segment of the trajectory used by the propagator
  FSGM_FakeCsIThickness = 100.;
  double safedistance = FSGM_FakeCsIThickness;
  FSGM_CalSafeBoundaries[0][0] = -2*FSGM_towerPitch-safedistance;
  FSGM_CalSafeBoundaries[0][1] = +2*FSGM_towerPitch+safedistance;
  FSGM_CalSafeBoundaries[1][0] = -2*FSGM_towerPitch-safedistance;
  FSGM_CalSafeBoundaries[1][1] = +2*FSGM_towerPitch+safedistance;
  FSGM_CalSafeBoundaries[2][0] = FSGM_calZBot-safedistance;
  FSGM_CalSafeBoundaries[2][1] = FSGM_calZTop+safedistance;

  // get the number of towers
  m_detSvc->getNumericConstByName("xNum", &FSGM_numX);
  m_detSvc->getNumericConstByName("yNum", &FSGM_numY);
  FSGM_flight_geom = true;
  if(FSGM_numX==4 && FSGM_numY==1) FSGM_flight_geom = false;

  // Fill the look-up table for even layers
  geom_xy_step = (4*FSGM_towerPitch)/(double)(FSGM_XY_MAX);
  
  double xyz[3];
  int whereincal[4];

  int i,j,k;
  // reference = middle of first layer
  xyz[2] = FSGM_calZTop-FSGM_cellVertPitch/2;
  
  // reference = middle of a log
  xyz[1] = FSGM_towerPitch/2+FSGM_cellHorPitch/2;
  for(i=0;i<FSGM_XY_MAX;++i)
    {
      xyz[0] = -2*FSGM_towerPitch + geom_xy_step*(0.5+(double)i);
      WhereInCalForGeom(xyz,whereincal);
      switch(whereincal[0])
        {
        case 0:
          geom_x_mat[i][0] = false;
          geom_x_mat[i][1] = false;
          break;
        case 1:
          geom_x_mat[i][0] = true;
          geom_x_mat[i][1] = false;
          break;
        case 2:
          geom_x_mat[i][0] = false;
          geom_x_mat[i][1] = true;
          break;
        default:
          geom_x_mat[i][0] = false;
          geom_x_mat[i][1] = false;
          break;
        }
      for(k=0;k<4;++k) geom_x_wic[i][k] = whereincal[k];
      geom_x_wic[i][2] = (int)floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
    }
  
  // reference = middle of a log
  xyz[0] = FSGM_towerPitch/2+ FSGM_cellHorPitch/2;
  for(j=0;j<FSGM_XY_MAX;++j)
    {
      xyz[1] = -2*FSGM_towerPitch+geom_xy_step*(0.5+(double)j);
      WhereInCalForGeom(xyz,whereincal);
      switch(whereincal[0])
        {
        case 0:
          geom_y_mat[j][0] = false;
          geom_y_mat[j][1] = false;
          break;
        case 1:
          geom_y_mat[j][0] = true;
          geom_y_mat[j][1] = false;
          break;
        case 2:
          geom_y_mat[j][0] = false;
          geom_y_mat[j][1] = true;
          break;
        default:
          geom_y_mat[j][0] = false;
          geom_y_mat[j][1] = false;
          break;
        }
      for(k=0;k<4;++k) geom_y_wic[j][k] = whereincal[k];
      geom_y_wic[j][2] = (int)floor((xyz[1]+2*FSGM_towerPitch)/FSGM_towerPitch);
    }

  // Fill the test points used by the propagator
  FillXYPoints();

  // Limits for the relative position in the shower to shower maximum
  FSGM_tmin = 0.25;
  FSGM_tmax = 3.0;

  // PAPER version
  FSGM_RCORE_0_0 = 0.023745;
  FSGM_RCORE_0_1 = -0.004212;
  FSGM_RCORE_1_0 = 0.091398;
  FSGM_RCORE_1_1 = 0.006449;
  FSGM_RTAIL_0_0 = 0.781676;
  FSGM_RTAIL_0_1 = -0.025275;
  FSGM_RTAIL_1_0 = 1.052475;
  FSGM_RTAIL_1_1 = 0.860926;
  FSGM_RTAIL_2_0 = 1.473248;
  FSGM_RTAIL_3_0 = 0.151222;
  FSGM_RTAIL_3_1 = 0.081055;
  FSGM_RTAIL_4_0 = 1.629281;
  FSGM_RTAIL_5_0 = 0.573951;
  FSGM_RTAIL_5_1 = 0.008947;
  FSGM_PCORE_0_0 = 1.191034;
  FSGM_PCORE_0_1 = -0.174791;
  FSGM_PCORE_0_2 = 0.055140;
  FSGM_PCORE_0_3 = -0.005380;
  FSGM_PCORE_1_0 = -0.343568;
  FSGM_PCORE_1_1 = 0.026725;
  FSGM_PCORE_1_2 = -0.007592;
  FSGM_PCORE_2_0 = 1.855664;
  FSGM_PCORE_2_1 = -1.699304;
  FSGM_PCORE_2_2 = 0.616402;
  FSGM_PCORE_2_3 = -0.070756;
  FSGM_PCORE_3_0 = -0.505006;
  FSGM_PCORE_3_1 = 1.455674;
  FSGM_PCORE_3_2 = -1.616056;
  FSGM_PCORE_3_3 = -6.568445;

  // Fill radial profile look-up table
  opthe = 0;
  FillRadialProfile(2,0);
  FillRadialProfile(3,1);
}

void NewFullShowerGeometryManager::SetOptHE(int opthein)
{
  if(opthein==0)
    {
      opthe = 0;
      PrepRadialProfile(2);
    }
  else
    {
      opthe = 1;
      PrepRadialProfile(3);
    }
}

void NewFullShowerGeometryManager::FillXYPoints()
{
  int i,j,k;

  int nring = FSGM_NPOINTS_NRING;
  double rmin = 0;
  double rmax = 1;
  double rstep = (rmax-rmin)/(double)nring;
  double radius;
  double x,y,step,d,d1;
  double stepmin;

  int n = FSGM_NPOINTS_IN_RING;
  //  printf("BRUEL nring %d n = %d\n",nring,n);

  radius = rmin+rstep;
  step = radius/(double)n;
  double weightmin = step*step;

  FSGM_NCircle = 0;

  bool myindex[6][21][21];
  for(k=0;k<nring;++k)
    for(i=0;i<=2*n;++i)
      for(j=0;j<=2*n;++j)
        myindex[k][i][j] = false;
  myindex[0][1][6] = true;
  myindex[0][1][7] = true;
  myindex[0][1][8] = true;
  myindex[0][1][9] = true;
  myindex[0][1][10] = true;
  myindex[0][1][11] = true;
  myindex[0][1][12] = true;
  myindex[0][1][13] = true;
  myindex[0][1][14] = true;
  myindex[0][2][4] = true;
  myindex[0][2][5] = true;
  myindex[0][2][6] = true;
  myindex[0][2][7] = true;
  myindex[0][2][8] = true;
  myindex[0][2][9] = true;
  myindex[0][2][10] = true;
  myindex[0][2][11] = true;
  myindex[0][2][12] = true;
  myindex[0][2][13] = true;
  myindex[0][2][14] = true;
  myindex[0][2][15] = true;
  myindex[0][3][3] = true;
  myindex[0][3][4] = true;
  myindex[0][3][5] = true;
  myindex[0][3][6] = true;
  myindex[0][3][7] = true;
  myindex[0][3][8] = true;
  myindex[0][3][9] = true;
  myindex[0][3][10] = true;
  myindex[0][3][11] = true;
  myindex[0][3][12] = true;
  myindex[0][3][13] = true;
  myindex[0][3][14] = true;
  myindex[0][3][15] = true;
  myindex[0][3][16] = true;
  myindex[0][3][17] = true;
  myindex[0][4][2] = true;
  myindex[0][4][3] = true;
  myindex[0][4][4] = true;
  myindex[0][4][5] = true;
  myindex[0][4][6] = true;
  myindex[0][4][7] = true;
  myindex[0][4][8] = true;
  myindex[0][4][9] = true;
  myindex[0][4][10] = true;
  myindex[0][4][11] = true;
  myindex[0][4][12] = true;
  myindex[0][4][13] = true;
  myindex[0][4][14] = true;
  myindex[0][4][15] = true;
  myindex[0][4][16] = true;
  myindex[0][4][17] = true;
  myindex[0][5][2] = true;
  myindex[0][5][3] = true;
  myindex[0][5][4] = true;
  myindex[0][5][5] = true;
  myindex[0][5][6] = true;
  myindex[0][5][7] = true;
  myindex[0][5][8] = true;
  myindex[0][5][9] = true;
  myindex[0][5][10] = true;
  myindex[0][5][11] = true;
  myindex[0][5][12] = true;
  myindex[0][5][13] = true;
  myindex[0][5][14] = true;
  myindex[0][5][15] = true;
  myindex[0][5][16] = true;
  myindex[0][5][17] = true;
  myindex[0][5][18] = true;
  myindex[0][6][1] = true;
  myindex[0][6][2] = true;
  myindex[0][6][3] = true;
  myindex[0][6][4] = true;
  myindex[0][6][5] = true;
  myindex[0][6][6] = true;
  myindex[0][6][7] = true;
  myindex[0][6][8] = true;
  myindex[0][6][9] = true;
  myindex[0][6][10] = true;
  myindex[0][6][11] = true;
  myindex[0][6][12] = true;
  myindex[0][6][13] = true;
  myindex[0][6][14] = true;
  myindex[0][6][15] = true;
  myindex[0][6][16] = true;
  myindex[0][6][17] = true;
  myindex[0][6][18] = true;
  myindex[0][6][19] = true;
  myindex[0][7][1] = true;
  myindex[0][7][2] = true;
  myindex[0][7][3] = true;
  myindex[0][7][4] = true;
  myindex[0][7][5] = true;
  myindex[0][7][6] = true;
  myindex[0][7][7] = true;
  myindex[0][7][8] = true;
  myindex[0][7][9] = true;
  myindex[0][7][10] = true;
  myindex[0][7][11] = true;
  myindex[0][7][12] = true;
  myindex[0][7][13] = true;
  myindex[0][7][14] = true;
  myindex[0][7][15] = true;
  myindex[0][7][16] = true;
  myindex[0][7][17] = true;
  myindex[0][7][18] = true;
  myindex[0][7][19] = true;
  myindex[0][8][1] = true;
  myindex[0][8][2] = true;
  myindex[0][8][3] = true;
  myindex[0][8][4] = true;
  myindex[0][8][5] = true;
  myindex[0][8][6] = true;
  myindex[0][8][7] = true;
  myindex[0][8][8] = true;
  myindex[0][8][9] = true;
  myindex[0][8][10] = true;
  myindex[0][8][11] = true;
  myindex[0][8][12] = true;
  myindex[0][8][13] = true;
  myindex[0][8][14] = true;
  myindex[0][8][15] = true;
  myindex[0][8][16] = true;
  myindex[0][8][17] = true;
  myindex[0][8][18] = true;
  myindex[0][8][19] = true;
  myindex[0][9][1] = true;
  myindex[0][9][2] = true;
  myindex[0][9][3] = true;
  myindex[0][9][4] = true;
  myindex[0][9][5] = true;
  myindex[0][9][6] = true;
  myindex[0][9][7] = true;
  myindex[0][9][8] = true;
  myindex[0][9][9] = true;
  myindex[0][9][10] = true;
  myindex[0][9][11] = true;
  myindex[0][9][12] = true;
  myindex[0][9][13] = true;
  myindex[0][9][14] = true;
  myindex[0][9][15] = true;
  myindex[0][9][16] = true;
  myindex[0][9][17] = true;
  myindex[0][9][18] = true;
  myindex[0][9][19] = true;
  myindex[0][10][0] = true;
  myindex[0][10][1] = true;
  myindex[0][10][2] = true;
  myindex[0][10][3] = true;
  myindex[0][10][4] = true;
  myindex[0][10][5] = true;
  myindex[0][10][6] = true;
  myindex[0][10][7] = true;
  myindex[0][10][8] = true;
  myindex[0][10][9] = true;
  myindex[0][10][10] = true;
  myindex[0][10][11] = true;
  myindex[0][10][12] = true;
  myindex[0][10][13] = true;
  myindex[0][10][14] = true;
  myindex[0][10][15] = true;
  myindex[0][10][16] = true;
  myindex[0][10][17] = true;
  myindex[0][10][18] = true;
  myindex[0][10][19] = true;
  myindex[0][11][1] = true;
  myindex[0][11][2] = true;
  myindex[0][11][3] = true;
  myindex[0][11][4] = true;
  myindex[0][11][5] = true;
  myindex[0][11][6] = true;
  myindex[0][11][7] = true;
  myindex[0][11][8] = true;
  myindex[0][11][9] = true;
  myindex[0][11][10] = true;
  myindex[0][11][11] = true;
  myindex[0][11][12] = true;
  myindex[0][11][13] = true;
  myindex[0][11][14] = true;
  myindex[0][11][15] = true;
  myindex[0][11][16] = true;
  myindex[0][11][17] = true;
  myindex[0][11][18] = true;
  myindex[0][11][19] = true;
  myindex[0][12][1] = true;
  myindex[0][12][2] = true;
  myindex[0][12][3] = true;
  myindex[0][12][4] = true;
  myindex[0][12][5] = true;
  myindex[0][12][6] = true;
  myindex[0][12][7] = true;
  myindex[0][12][8] = true;
  myindex[0][12][9] = true;
  myindex[0][12][10] = true;
  myindex[0][12][11] = true;
  myindex[0][12][12] = true;
  myindex[0][12][13] = true;
  myindex[0][12][14] = true;
  myindex[0][12][15] = true;
  myindex[0][12][16] = true;
  myindex[0][12][17] = true;
  myindex[0][12][18] = true;
  myindex[0][12][19] = true;
  myindex[0][13][1] = true;
  myindex[0][13][2] = true;
  myindex[0][13][3] = true;
  myindex[0][13][4] = true;
  myindex[0][13][5] = true;
  myindex[0][13][6] = true;
  myindex[0][13][7] = true;
  myindex[0][13][8] = true;
  myindex[0][13][9] = true;
  myindex[0][13][10] = true;
  myindex[0][13][11] = true;
  myindex[0][13][12] = true;
  myindex[0][13][13] = true;
  myindex[0][13][14] = true;
  myindex[0][13][15] = true;
  myindex[0][13][16] = true;
  myindex[0][13][17] = true;
  myindex[0][13][18] = true;
  myindex[0][13][19] = true;
  myindex[0][14][1] = true;
  myindex[0][14][2] = true;
  myindex[0][14][3] = true;
  myindex[0][14][4] = true;
  myindex[0][14][5] = true;
  myindex[0][14][6] = true;
  myindex[0][14][7] = true;
  myindex[0][14][8] = true;
  myindex[0][14][9] = true;
  myindex[0][14][10] = true;
  myindex[0][14][11] = true;
  myindex[0][14][12] = true;
  myindex[0][14][13] = true;
  myindex[0][14][14] = true;
  myindex[0][14][15] = true;
  myindex[0][14][16] = true;
  myindex[0][14][17] = true;
  myindex[0][14][18] = true;
  myindex[0][14][19] = true;
  myindex[0][15][2] = true;
  myindex[0][15][3] = true;
  myindex[0][15][4] = true;
  myindex[0][15][5] = true;
  myindex[0][15][6] = true;
  myindex[0][15][7] = true;
  myindex[0][15][8] = true;
  myindex[0][15][9] = true;
  myindex[0][15][10] = true;
  myindex[0][15][11] = true;
  myindex[0][15][12] = true;
  myindex[0][15][13] = true;
  myindex[0][15][14] = true;
  myindex[0][15][15] = true;
  myindex[0][15][16] = true;
  myindex[0][15][17] = true;
  myindex[0][15][18] = true;
  myindex[0][16][3] = true;
  myindex[0][16][4] = true;
  myindex[0][16][5] = true;
  myindex[0][16][6] = true;
  myindex[0][16][7] = true;
  myindex[0][16][8] = true;
  myindex[0][16][9] = true;
  myindex[0][16][10] = true;
  myindex[0][16][11] = true;
  myindex[0][16][12] = true;
  myindex[0][16][13] = true;
  myindex[0][16][14] = true;
  myindex[0][16][15] = true;
  myindex[0][16][16] = true;
  myindex[0][16][17] = true;
  myindex[0][17][3] = true;
  myindex[0][17][4] = true;
  myindex[0][17][5] = true;
  myindex[0][17][6] = true;
  myindex[0][17][7] = true;
  myindex[0][17][8] = true;
  myindex[0][17][9] = true;
  myindex[0][17][10] = true;
  myindex[0][17][11] = true;
  myindex[0][17][12] = true;
  myindex[0][17][13] = true;
  myindex[0][17][14] = true;
  myindex[0][17][15] = true;
  myindex[0][17][16] = true;
  myindex[0][17][17] = true;
  myindex[0][18][5] = true;
  myindex[0][18][6] = true;
  myindex[0][18][7] = true;
  myindex[0][18][8] = true;
  myindex[0][18][9] = true;
  myindex[0][18][10] = true;
  myindex[0][18][11] = true;
  myindex[0][18][12] = true;
  myindex[0][18][13] = true;
  myindex[0][18][14] = true;
  myindex[0][18][15] = true;
  myindex[0][19][6] = true;
  myindex[0][19][7] = true;
  myindex[0][19][8] = true;
  myindex[0][19][9] = true;
  myindex[0][19][10] = true;
  myindex[0][19][11] = true;
  myindex[0][19][12] = true;
  myindex[0][19][13] = true;
  myindex[0][19][14] = true;
  myindex[1][1][6] = true;
  myindex[1][1][7] = true;
  myindex[1][1][8] = true;
  myindex[1][1][9] = true;
  myindex[1][1][10] = true;
  myindex[1][1][11] = true;
  myindex[1][1][12] = true;
  myindex[1][1][13] = true;
  myindex[1][1][14] = true;
  myindex[1][2][4] = true;
  myindex[1][2][5] = true;
  myindex[1][2][6] = true;
  myindex[1][2][7] = true;
  myindex[1][2][8] = true;
  myindex[1][2][9] = true;
  myindex[1][2][10] = true;
  myindex[1][2][11] = true;
  myindex[1][2][12] = true;
  myindex[1][2][13] = true;
  myindex[1][2][14] = true;
  myindex[1][2][15] = true;
  myindex[1][3][3] = true;
  myindex[1][3][4] = true;
  myindex[1][3][5] = true;
  myindex[1][3][6] = true;
  myindex[1][3][7] = true;
  myindex[1][3][8] = true;
  myindex[1][3][9] = true;
  myindex[1][3][10] = true;
  myindex[1][3][11] = true;
  myindex[1][3][12] = true;
  myindex[1][3][13] = true;
  myindex[1][3][14] = true;
  myindex[1][3][15] = true;
  myindex[1][3][16] = true;
  myindex[1][3][17] = true;
  myindex[1][4][2] = true;
  myindex[1][4][3] = true;
  myindex[1][4][4] = true;
  myindex[1][4][5] = true;
  myindex[1][4][6] = true;
  myindex[1][4][7] = true;
  myindex[1][4][8] = true;
  myindex[1][4][9] = true;
  myindex[1][4][10] = true;
  myindex[1][4][11] = true;
  myindex[1][4][12] = true;
  myindex[1][4][13] = true;
  myindex[1][4][14] = true;
  myindex[1][4][15] = true;
  myindex[1][4][16] = true;
  myindex[1][4][17] = true;
  myindex[1][5][2] = true;
  myindex[1][5][3] = true;
  myindex[1][5][4] = true;
  myindex[1][5][5] = true;
  myindex[1][5][6] = true;
  myindex[1][5][7] = true;
  myindex[1][5][8] = true;
  myindex[1][5][9] = true;
  myindex[1][5][11] = true;
  myindex[1][5][12] = true;
  myindex[1][5][13] = true;
  myindex[1][5][14] = true;
  myindex[1][5][15] = true;
  myindex[1][5][16] = true;
  myindex[1][5][17] = true;
  myindex[1][5][18] = true;
  myindex[1][6][1] = true;
  myindex[1][6][2] = true;
  myindex[1][6][3] = true;
  myindex[1][6][4] = true;
  myindex[1][6][5] = true;
  myindex[1][6][6] = true;
  myindex[1][6][13] = true;
  myindex[1][6][14] = true;
  myindex[1][6][15] = true;
  myindex[1][6][16] = true;
  myindex[1][6][17] = true;
  myindex[1][6][18] = true;
  myindex[1][6][19] = true;
  myindex[1][7][1] = true;
  myindex[1][7][2] = true;
  myindex[1][7][3] = true;
  myindex[1][7][4] = true;
  myindex[1][7][5] = true;
  myindex[1][7][14] = true;
  myindex[1][7][15] = true;
  myindex[1][7][16] = true;
  myindex[1][7][17] = true;
  myindex[1][7][18] = true;
  myindex[1][7][19] = true;
  myindex[1][8][1] = true;
  myindex[1][8][2] = true;
  myindex[1][8][3] = true;
  myindex[1][8][4] = true;
  myindex[1][8][5] = true;
  myindex[1][8][15] = true;
  myindex[1][8][16] = true;
  myindex[1][8][17] = true;
  myindex[1][8][18] = true;
  myindex[1][8][19] = true;
  myindex[1][9][1] = true;
  myindex[1][9][2] = true;
  myindex[1][9][3] = true;
  myindex[1][9][4] = true;
  myindex[1][9][5] = true;
  myindex[1][9][15] = true;
  myindex[1][9][16] = true;
  myindex[1][9][17] = true;
  myindex[1][9][18] = true;
  myindex[1][9][19] = true;
  myindex[1][10][0] = true;
  myindex[1][10][1] = true;
  myindex[1][10][2] = true;
  myindex[1][10][3] = true;
  myindex[1][10][4] = true;
  myindex[1][10][15] = true;
  myindex[1][10][16] = true;
  myindex[1][10][17] = true;
  myindex[1][10][18] = true;
  myindex[1][10][19] = true;
  myindex[1][11][1] = true;
  myindex[1][11][2] = true;
  myindex[1][11][3] = true;
  myindex[1][11][4] = true;
  myindex[1][11][5] = true;
  myindex[1][11][15] = true;
  myindex[1][11][16] = true;
  myindex[1][11][17] = true;
  myindex[1][11][18] = true;
  myindex[1][11][19] = true;
  myindex[1][12][1] = true;
  myindex[1][12][2] = true;
  myindex[1][12][3] = true;
  myindex[1][12][4] = true;
  myindex[1][12][5] = true;
  myindex[1][12][15] = true;
  myindex[1][12][16] = true;
  myindex[1][12][17] = true;
  myindex[1][12][18] = true;
  myindex[1][12][19] = true;
  myindex[1][13][1] = true;
  myindex[1][13][2] = true;
  myindex[1][13][3] = true;
  myindex[1][13][4] = true;
  myindex[1][13][5] = true;
  myindex[1][13][6] = true;
  myindex[1][13][14] = true;
  myindex[1][13][15] = true;
  myindex[1][13][16] = true;
  myindex[1][13][17] = true;
  myindex[1][13][18] = true;
  myindex[1][13][19] = true;
  myindex[1][14][1] = true;
  myindex[1][14][2] = true;
  myindex[1][14][3] = true;
  myindex[1][14][4] = true;
  myindex[1][14][5] = true;
  myindex[1][14][6] = true;
  myindex[1][14][7] = true;
  myindex[1][14][13] = true;
  myindex[1][14][14] = true;
  myindex[1][14][15] = true;
  myindex[1][14][16] = true;
  myindex[1][14][17] = true;
  myindex[1][14][18] = true;
  myindex[1][14][19] = true;
  myindex[1][15][2] = true;
  myindex[1][15][3] = true;
  myindex[1][15][4] = true;
  myindex[1][15][5] = true;
  myindex[1][15][6] = true;
  myindex[1][15][7] = true;
  myindex[1][15][8] = true;
  myindex[1][15][9] = true;
  myindex[1][15][10] = true;
  myindex[1][15][11] = true;
  myindex[1][15][12] = true;
  myindex[1][15][13] = true;
  myindex[1][15][14] = true;
  myindex[1][15][15] = true;
  myindex[1][15][16] = true;
  myindex[1][15][17] = true;
  myindex[1][15][18] = true;
  myindex[1][16][3] = true;
  myindex[1][16][4] = true;
  myindex[1][16][5] = true;
  myindex[1][16][6] = true;
  myindex[1][16][7] = true;
  myindex[1][16][8] = true;
  myindex[1][16][9] = true;
  myindex[1][16][10] = true;
  myindex[1][16][11] = true;
  myindex[1][16][12] = true;
  myindex[1][16][13] = true;
  myindex[1][16][14] = true;
  myindex[1][16][15] = true;
  myindex[1][16][16] = true;
  myindex[1][16][17] = true;
  myindex[1][17][3] = true;
  myindex[1][17][4] = true;
  myindex[1][17][5] = true;
  myindex[1][17][6] = true;
  myindex[1][17][7] = true;
  myindex[1][17][8] = true;
  myindex[1][17][9] = true;
  myindex[1][17][10] = true;
  myindex[1][17][11] = true;
  myindex[1][17][12] = true;
  myindex[1][17][13] = true;
  myindex[1][17][14] = true;
  myindex[1][17][15] = true;
  myindex[1][17][16] = true;
  myindex[1][17][17] = true;
  myindex[1][18][5] = true;
  myindex[1][18][6] = true;
  myindex[1][18][7] = true;
  myindex[1][18][8] = true;
  myindex[1][18][9] = true;
  myindex[1][18][10] = true;
  myindex[1][18][11] = true;
  myindex[1][18][12] = true;
  myindex[1][18][13] = true;
  myindex[1][18][14] = true;
  myindex[1][18][15] = true;
  myindex[1][19][6] = true;
  myindex[1][19][7] = true;
  myindex[1][19][8] = true;
  myindex[1][19][9] = true;
  myindex[1][19][10] = true;
  myindex[1][19][11] = true;
  myindex[1][19][12] = true;
  myindex[1][19][13] = true;
  myindex[1][19][14] = true;
  myindex[2][0][10] = true;
  myindex[2][1][6] = true;
  myindex[2][1][7] = true;
  myindex[2][1][8] = true;
  myindex[2][1][9] = true;
  myindex[2][1][10] = true;
  myindex[2][1][11] = true;
  myindex[2][1][12] = true;
  myindex[2][1][13] = true;
  myindex[2][1][14] = true;
  myindex[2][2][4] = true;
  myindex[2][2][5] = true;
  myindex[2][2][6] = true;
  myindex[2][2][7] = true;
  myindex[2][2][8] = true;
  myindex[2][2][9] = true;
  myindex[2][2][10] = true;
  myindex[2][2][11] = true;
  myindex[2][2][12] = true;
  myindex[2][2][13] = true;
  myindex[2][2][14] = true;
  myindex[2][2][15] = true;
  myindex[2][3][3] = true;
  myindex[2][3][4] = true;
  myindex[2][3][5] = true;
  myindex[2][3][6] = true;
  myindex[2][3][7] = true;
  myindex[2][3][8] = true;
  myindex[2][3][9] = true;
  myindex[2][3][10] = true;
  myindex[2][3][11] = true;
  myindex[2][3][12] = true;
  myindex[2][3][13] = true;
  myindex[2][3][14] = true;
  myindex[2][3][15] = true;
  myindex[2][3][16] = true;
  myindex[2][3][17] = true;
  myindex[2][4][2] = true;
  myindex[2][4][3] = true;
  myindex[2][4][4] = true;
  myindex[2][4][5] = true;
  myindex[2][4][6] = true;
  myindex[2][4][7] = true;
  myindex[2][4][13] = true;
  myindex[2][4][14] = true;
  myindex[2][4][15] = true;
  myindex[2][4][16] = true;
  myindex[2][4][17] = true;
  myindex[2][5][2] = true;
  myindex[2][5][3] = true;
  myindex[2][5][4] = true;
  myindex[2][5][5] = true;
  myindex[2][5][15] = true;
  myindex[2][5][16] = true;
  myindex[2][5][17] = true;
  myindex[2][5][18] = true;
  myindex[2][6][1] = true;
  myindex[2][6][2] = true;
  myindex[2][6][3] = true;
  myindex[2][6][4] = true;
  myindex[2][6][16] = true;
  myindex[2][6][17] = true;
  myindex[2][6][18] = true;
  myindex[2][6][19] = true;
  myindex[2][7][1] = true;
  myindex[2][7][2] = true;
  myindex[2][7][3] = true;
  myindex[2][7][4] = true;
  myindex[2][7][16] = true;
  myindex[2][7][17] = true;
  myindex[2][7][18] = true;
  myindex[2][7][19] = true;
  myindex[2][8][1] = true;
  myindex[2][8][2] = true;
  myindex[2][8][3] = true;
  myindex[2][8][17] = true;
  myindex[2][8][18] = true;
  myindex[2][8][19] = true;
  myindex[2][9][1] = true;
  myindex[2][9][2] = true;
  myindex[2][9][3] = true;
  myindex[2][9][17] = true;
  myindex[2][9][18] = true;
  myindex[2][9][19] = true;
  myindex[2][10][0] = true;
  myindex[2][10][1] = true;
  myindex[2][10][2] = true;
  myindex[2][10][3] = true;
  myindex[2][10][17] = true;
  myindex[2][10][18] = true;
  myindex[2][10][19] = true;
  myindex[2][11][1] = true;
  myindex[2][11][2] = true;
  myindex[2][11][3] = true;
  myindex[2][11][17] = true;
  myindex[2][11][18] = true;
  myindex[2][11][19] = true;
  myindex[2][12][1] = true;
  myindex[2][12][2] = true;
  myindex[2][12][3] = true;
  myindex[2][12][17] = true;
  myindex[2][12][18] = true;
  myindex[2][12][19] = true;
  myindex[2][13][1] = true;
  myindex[2][13][2] = true;
  myindex[2][13][3] = true;
  myindex[2][13][4] = true;
  myindex[2][13][16] = true;
  myindex[2][13][17] = true;
  myindex[2][13][18] = true;
  myindex[2][13][19] = true;
  myindex[2][14][1] = true;
  myindex[2][14][2] = true;
  myindex[2][14][3] = true;
  myindex[2][14][4] = true;
  myindex[2][14][16] = true;
  myindex[2][14][17] = true;
  myindex[2][14][18] = true;
  myindex[2][14][19] = true;
  myindex[2][15][2] = true;
  myindex[2][15][3] = true;
  myindex[2][15][4] = true;
  myindex[2][15][5] = true;
  myindex[2][15][15] = true;
  myindex[2][15][16] = true;
  myindex[2][15][17] = true;
  myindex[2][15][18] = true;
  myindex[2][16][3] = true;
  myindex[2][16][4] = true;
  myindex[2][16][5] = true;
  myindex[2][16][6] = true;
  myindex[2][16][7] = true;
  myindex[2][16][13] = true;
  myindex[2][16][14] = true;
  myindex[2][16][15] = true;
  myindex[2][16][16] = true;
  myindex[2][16][17] = true;
  myindex[2][17][3] = true;
  myindex[2][17][4] = true;
  myindex[2][17][5] = true;
  myindex[2][17][6] = true;
  myindex[2][17][7] = true;
  myindex[2][17][8] = true;
  myindex[2][17][9] = true;
  myindex[2][17][10] = true;
  myindex[2][17][11] = true;
  myindex[2][17][12] = true;
  myindex[2][17][13] = true;
  myindex[2][17][14] = true;
  myindex[2][17][15] = true;
  myindex[2][17][16] = true;
  myindex[2][17][17] = true;
  myindex[2][18][5] = true;
  myindex[2][18][6] = true;
  myindex[2][18][7] = true;
  myindex[2][18][8] = true;
  myindex[2][18][9] = true;
  myindex[2][18][10] = true;
  myindex[2][18][11] = true;
  myindex[2][18][12] = true;
  myindex[2][18][13] = true;
  myindex[2][18][14] = true;
  myindex[2][18][15] = true;
  myindex[2][19][6] = true;
  myindex[2][19][7] = true;
  myindex[2][19][8] = true;
  myindex[2][19][9] = true;
  myindex[2][19][10] = true;
  myindex[2][19][11] = true;
  myindex[2][19][12] = true;
  myindex[2][19][13] = true;
  myindex[2][19][14] = true;
  myindex[3][1][6] = true;
  myindex[3][1][7] = true;
  myindex[3][1][8] = true;
  myindex[3][1][9] = true;
  myindex[3][1][10] = true;
  myindex[3][1][11] = true;
  myindex[3][1][12] = true;
  myindex[3][1][13] = true;
  myindex[3][1][14] = true;
  myindex[3][2][4] = true;
  myindex[3][2][5] = true;
  myindex[3][2][6] = true;
  myindex[3][2][7] = true;
  myindex[3][2][8] = true;
  myindex[3][2][9] = true;
  myindex[3][2][10] = true;
  myindex[3][2][11] = true;
  myindex[3][2][12] = true;
  myindex[3][2][13] = true;
  myindex[3][2][14] = true;
  myindex[3][2][15] = true;
  myindex[3][3][3] = true;
  myindex[3][3][4] = true;
  myindex[3][3][5] = true;
  myindex[3][3][6] = true;
  myindex[3][3][7] = true;
  myindex[3][3][13] = true;
  myindex[3][3][14] = true;
  myindex[3][3][15] = true;
  myindex[3][3][16] = true;
  myindex[3][3][17] = true;
  myindex[3][4][2] = true;
  myindex[3][4][3] = true;
  myindex[3][4][4] = true;
  myindex[3][4][5] = true;
  myindex[3][4][15] = true;
  myindex[3][4][16] = true;
  myindex[3][4][17] = true;
  myindex[3][5][2] = true;
  myindex[3][5][3] = true;
  myindex[3][5][4] = true;
  myindex[3][5][16] = true;
  myindex[3][5][17] = true;
  myindex[3][5][18] = true;
  myindex[3][6][1] = true;
  myindex[3][6][2] = true;
  myindex[3][6][3] = true;
  myindex[3][6][17] = true;
  myindex[3][6][18] = true;
  myindex[3][6][19] = true;
  myindex[3][7][1] = true;
  myindex[3][7][2] = true;
  myindex[3][7][3] = true;
  myindex[3][7][17] = true;
  myindex[3][7][18] = true;
  myindex[3][7][19] = true;
  myindex[3][8][1] = true;
  myindex[3][8][2] = true;
  myindex[3][8][18] = true;
  myindex[3][8][19] = true;
  myindex[3][9][1] = true;
  myindex[3][9][2] = true;
  myindex[3][9][18] = true;
  myindex[3][9][19] = true;
  myindex[3][10][0] = true;
  myindex[3][10][1] = true;
  myindex[3][10][2] = true;
  myindex[3][10][18] = true;
  myindex[3][10][19] = true;
  myindex[3][11][1] = true;
  myindex[3][11][2] = true;
  myindex[3][11][18] = true;
  myindex[3][11][19] = true;
  myindex[3][12][1] = true;
  myindex[3][12][2] = true;
  myindex[3][12][18] = true;
  myindex[3][12][19] = true;
  myindex[3][13][1] = true;
  myindex[3][13][2] = true;
  myindex[3][13][3] = true;
  myindex[3][13][17] = true;
  myindex[3][13][18] = true;
  myindex[3][13][19] = true;
  myindex[3][14][1] = true;
  myindex[3][14][2] = true;
  myindex[3][14][3] = true;
  myindex[3][14][17] = true;
  myindex[3][14][18] = true;
  myindex[3][14][19] = true;
  myindex[3][15][2] = true;
  myindex[3][15][3] = true;
  myindex[3][15][4] = true;
  myindex[3][15][16] = true;
  myindex[3][15][17] = true;
  myindex[3][15][18] = true;
  myindex[3][16][3] = true;
  myindex[3][16][4] = true;
  myindex[3][16][5] = true;
  myindex[3][16][15] = true;
  myindex[3][16][16] = true;
  myindex[3][16][17] = true;
  myindex[3][17][3] = true;
  myindex[3][17][4] = true;
  myindex[3][17][5] = true;
  myindex[3][17][6] = true;
  myindex[3][17][7] = true;
  myindex[3][17][13] = true;
  myindex[3][17][14] = true;
  myindex[3][17][15] = true;
  myindex[3][17][16] = true;
  myindex[3][17][17] = true;
  myindex[3][18][5] = true;
  myindex[3][18][6] = true;
  myindex[3][18][7] = true;
  myindex[3][18][8] = true;
  myindex[3][18][9] = true;
  myindex[3][18][10] = true;
  myindex[3][18][11] = true;
  myindex[3][18][12] = true;
  myindex[3][18][13] = true;
  myindex[3][18][14] = true;
  myindex[3][18][15] = true;
  myindex[3][19][6] = true;
  myindex[3][19][7] = true;
  myindex[3][19][8] = true;
  myindex[3][19][9] = true;
  myindex[3][19][10] = true;
  myindex[3][19][11] = true;
  myindex[3][19][12] = true;
  myindex[3][19][13] = true;
  myindex[3][19][14] = true;
  myindex[4][0][10] = true;
  myindex[4][1][6] = true;
  myindex[4][1][7] = true;
  myindex[4][1][8] = true;
  myindex[4][1][9] = true;
  myindex[4][1][10] = true;
  myindex[4][1][11] = true;
  myindex[4][1][12] = true;
  myindex[4][1][13] = true;
  myindex[4][1][14] = true;
  myindex[4][2][4] = true;
  myindex[4][2][5] = true;
  myindex[4][2][6] = true;
  myindex[4][2][7] = true;
  myindex[4][2][8] = true;
  myindex[4][2][9] = true;
  myindex[4][2][10] = true;
  myindex[4][2][11] = true;
  myindex[4][2][12] = true;
  myindex[4][2][13] = true;
  myindex[4][2][14] = true;
  myindex[4][2][15] = true;
  myindex[4][3][3] = true;
  myindex[4][3][4] = true;
  myindex[4][3][5] = true;
  myindex[4][3][6] = true;
  myindex[4][3][14] = true;
  myindex[4][3][15] = true;
  myindex[4][3][16] = true;
  myindex[4][3][17] = true;
  myindex[4][4][2] = true;
  myindex[4][4][3] = true;
  myindex[4][4][4] = true;
  myindex[4][4][16] = true;
  myindex[4][4][17] = true;
  myindex[4][5][2] = true;
  myindex[4][5][3] = true;
  myindex[4][5][17] = true;
  myindex[4][5][18] = true;
  myindex[4][6][1] = true;
  myindex[4][6][2] = true;
  myindex[4][6][3] = true;
  myindex[4][6][17] = true;
  myindex[4][6][18] = true;
  myindex[4][6][19] = true;
  myindex[4][7][1] = true;
  myindex[4][7][2] = true;
  myindex[4][7][18] = true;
  myindex[4][7][19] = true;
  myindex[4][8][1] = true;
  myindex[4][8][2] = true;
  myindex[4][8][18] = true;
  myindex[4][8][19] = true;
  myindex[4][9][1] = true;
  myindex[4][9][2] = true;
  myindex[4][9][18] = true;
  myindex[4][9][19] = true;
  myindex[4][10][0] = true;
  myindex[4][10][1] = true;
  myindex[4][10][2] = true;
  myindex[4][10][18] = true;
  myindex[4][10][19] = true;
  myindex[4][11][1] = true;
  myindex[4][11][2] = true;
  myindex[4][11][18] = true;
  myindex[4][11][19] = true;
  myindex[4][12][1] = true;
  myindex[4][12][2] = true;
  myindex[4][12][18] = true;
  myindex[4][12][19] = true;
  myindex[4][13][1] = true;
  myindex[4][13][2] = true;
  myindex[4][13][18] = true;
  myindex[4][13][19] = true;
  myindex[4][14][1] = true;
  myindex[4][14][2] = true;
  myindex[4][14][3] = true;
  myindex[4][14][17] = true;
  myindex[4][14][18] = true;
  myindex[4][14][19] = true;
  myindex[4][15][2] = true;
  myindex[4][15][3] = true;
  myindex[4][15][17] = true;
  myindex[4][15][18] = true;
  myindex[4][16][3] = true;
  myindex[4][16][4] = true;
  myindex[4][16][16] = true;
  myindex[4][16][17] = true;
  myindex[4][17][3] = true;
  myindex[4][17][4] = true;
  myindex[4][17][5] = true;
  myindex[4][17][6] = true;
  myindex[4][17][14] = true;
  myindex[4][17][15] = true;
  myindex[4][17][16] = true;
  myindex[4][17][17] = true;
  myindex[4][18][5] = true;
  myindex[4][18][6] = true;
  myindex[4][18][7] = true;
  myindex[4][18][8] = true;
  myindex[4][18][9] = true;
  myindex[4][18][10] = true;
  myindex[4][18][11] = true;
  myindex[4][18][12] = true;
  myindex[4][18][13] = true;
  myindex[4][18][14] = true;
  myindex[4][18][15] = true;
  myindex[4][19][6] = true;
  myindex[4][19][7] = true;
  myindex[4][19][8] = true;
  myindex[4][19][9] = true;
  myindex[4][19][10] = true;
  myindex[4][19][11] = true;
  myindex[4][19][12] = true;
  myindex[4][19][13] = true;
  myindex[4][19][14] = true;
  myindex[5][0][10] = true;
  myindex[5][1][6] = true;
  myindex[5][1][7] = true;
  myindex[5][1][8] = true;
  myindex[5][1][9] = true;
  myindex[5][1][10] = true;
  myindex[5][1][11] = true;
  myindex[5][1][12] = true;
  myindex[5][1][13] = true;
  myindex[5][1][14] = true;
  myindex[5][2][4] = true;
  myindex[5][2][5] = true;
  myindex[5][2][6] = true;
  myindex[5][2][7] = true;
  myindex[5][2][13] = true;
  myindex[5][2][14] = true;
  myindex[5][2][15] = true;
  myindex[5][3][3] = true;
  myindex[5][3][4] = true;
  myindex[5][3][5] = true;
  myindex[5][3][15] = true;
  myindex[5][3][16] = true;
  myindex[5][3][17] = true;
  myindex[5][4][2] = true;
  myindex[5][4][3] = true;
  myindex[5][4][4] = true;
  myindex[5][4][16] = true;
  myindex[5][4][17] = true;
  myindex[5][5][2] = true;
  myindex[5][5][3] = true;
  myindex[5][5][17] = true;
  myindex[5][5][18] = true;
  myindex[5][6][1] = true;
  myindex[5][6][2] = true;
  myindex[5][6][18] = true;
  myindex[5][6][19] = true;
  myindex[5][7][1] = true;
  myindex[5][7][2] = true;
  myindex[5][7][18] = true;
  myindex[5][7][19] = true;
  myindex[5][8][1] = true;
  myindex[5][8][19] = true;
  myindex[5][9][1] = true;
  myindex[5][9][19] = true;
  myindex[5][10][0] = true;
  myindex[5][10][1] = true;
  myindex[5][10][19] = true;
  myindex[5][11][1] = true;
  myindex[5][11][19] = true;
  myindex[5][12][1] = true;
  myindex[5][12][19] = true;
  myindex[5][13][1] = true;
  myindex[5][13][2] = true;
  myindex[5][13][18] = true;
  myindex[5][13][19] = true;
  myindex[5][14][1] = true;
  myindex[5][14][2] = true;
  myindex[5][14][18] = true;
  myindex[5][14][19] = true;
  myindex[5][15][2] = true;
  myindex[5][15][3] = true;
  myindex[5][15][17] = true;
  myindex[5][15][18] = true;
  myindex[5][16][3] = true;
  myindex[5][16][4] = true;
  myindex[5][16][16] = true;
  myindex[5][16][17] = true;
  myindex[5][17][3] = true;
  myindex[5][17][4] = true;
  myindex[5][17][5] = true;
  myindex[5][17][15] = true;
  myindex[5][17][16] = true;
  myindex[5][17][17] = true;
  myindex[5][18][5] = true;
  myindex[5][18][6] = true;
  myindex[5][18][7] = true;
  myindex[5][18][13] = true;
  myindex[5][18][14] = true;
  myindex[5][18][15] = true;
  myindex[5][19][6] = true;
  myindex[5][19][7] = true;
  myindex[5][19][8] = true;
  myindex[5][19][9] = true;
  myindex[5][19][10] = true;
  myindex[5][19][11] = true;
  myindex[5][19][12] = true;
  myindex[5][19][13] = true;
  myindex[5][19][14] = true;

  for(k=0;k<nring;++k)
    {
      radius = rmin+rstep*(1.+(double)k);
      step = radius/(double)n;
      //      printf("BRUEL %d radius %f (%f RM %f mm) step = %f (%f RM %f mm) -> %f\n",k,radius,radius*3.,radius*3.*35.,step,step*3.,step*3.*35.,step*step/weightmin);
      //
      for(i=0;i<=2*n;++i)
        for(j=0;j<=2*n;++j)
          {
            x = -radius+step*(double)i;
            y = -radius+step*(double)j;
            d = sqrt(x*x+y*y);
            d1 = d;
            //
            // PhB 20141202: point definition is now hardcoded in order to ensure that we get the exactly the same points in rhel6 as in rhel5
            //            if(d>radius || (k>0 && d<radius-rstep)) continue;
            if(myindex[k][i][j]==false) continue;
            //
            FSGM_XCircle[FSGM_NCircle] = x;
            FSGM_YCircle[FSGM_NCircle] = y;
            FSGM_RCircle[FSGM_NCircle] = d;
            FSGM_PCircle[FSGM_NCircle] = atan2(y,x);
            FSGM_CPCircle[FSGM_NCircle] = cos(FSGM_PCircle[FSGM_NCircle]);
            FSGM_SPCircle[FSGM_NCircle] = sin(FSGM_PCircle[FSGM_NCircle]);
            FSGM_WCircle[FSGM_NCircle] = step*step/weightmin;
            //   printf("BRUEL x = %f %f %f %f\n",x,y,d,step*step/weightmin);
            //   printf("BRUEL FSGM_PCircle[%d] = %f %f %f %f %f %f\n",FSGM_NCircle,x,y,d,FSGM_PCircle[FSGM_NCircle],FSGM_CPCircle[FSGM_NCircle],FSGM_SPCircle[FSGM_NCircle]);
            ++FSGM_NCircle;
          }
      //      printf("BRUEL FSGM_NCircle = %d\n",FSGM_NCircle);
    }
}

void NewFullShowerGeometryManager::oldFillXYPoints()
{
  int n = (int)sqrt(1./3.2*(double)FSGM_NPOINTS_MAX);
  int i,j;
  double step = 1./(double)(n);

  double x,y;
  FSGM_NCircle = 0;
  for(i=0;i<=2*n;++i)
    {
      x = -1+step*(double)i;
      for(j=0;j<=2*n;++j)
        {
          y = -1+step*(double)j;
          if(sqrt(x*x+y*y)>1) continue;
          FSGM_XCircle[FSGM_NCircle] = x;
          FSGM_YCircle[FSGM_NCircle] = y;
          FSGM_RCircle[FSGM_NCircle] = sqrt(x*x+y*y);
          FSGM_WCircle[FSGM_NCircle] = 1;
          ++FSGM_NCircle;
        }
    }
}

void NewFullShowerGeometryManager::WhereInCalForGeom(double *xyz, int *whereincal)
{
  whereincal[0] = 3; // material
  whereincal[1] = -1; // layer
  whereincal[2] = 0; // tower
  whereincal[3] = 0; // column
  //

  if(fabs(xyz[0])>=2*FSGM_towerPitch) return;
  if(fabs(xyz[1])>=2*FSGM_towerPitch) return;
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;
  whereincal[0] = 0;
  //
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  whereincal[1] = ilayer;
  //  if(fabs(xyz[2]- (FSGM_calZTop-FSGM_cellVertPitch*(0.5+(double)ilayer)))>FSGM_CsIHeight*0.5) return;
  //
  whereincal[2] = (int)floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  whereincal[3] = (int)floor((xyz[1]+2*FSGM_towerPitch)/FSGM_towerPitch);
  whereincal[2] = 4*whereincal[3]+whereincal[2];
  //
  double xmod = xyz[0]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double ymod = xyz[1]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[1]+2*FSGM_towerPitch)/FSGM_towerPitch);
  if(ilayer%2==0)
    whereincal[3] = (int)floor((ymod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  else
    whereincal[3] = (int)floor((xmod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  //

  if(FSGM_towerPitch-xmod<xmod) xmod = FSGM_towerPitch-xmod;
  if(FSGM_towerPitch-ymod<ymod) ymod = FSGM_towerPitch-ymod;

  double xpara = xmod;
  double xperp = ymod;
  if(ilayer%2==1)
    {
      xpara = ymod;
      xperp = xmod;
    }

  if(xpara<FSGM_crackPara || xperp<FSGM_crackPerp)
    {
      whereincal[0] = 2;
      return;
    }

  xperp -= (FSGM_towerPitch/2-6*FSGM_cellHorPitch);
  double xymod = xperp-FSGM_cellHorPitch*floor(xperp/FSGM_cellHorPitch);
  if(FSGM_cellHorPitch-xymod<xymod) xymod = FSGM_cellHorPitch-xymod;

  if(xymod<FSGM_crackXtal) return;

  //   double zlayer = (FSGM_calZTop-xyz[2])-FSGM_cellVertPitch*(double)ilayer;
  //   if(FSGM_cellVertPitch-zlayer>zlayer) zlayer = FSGM_cellVertPitch-zlayer;
  //   if(zlayer<(FSGM_cellVertPitch-FSGM_CsIHeight)/2) return;
  whereincal[0] = 1;
  whereincal[1] = ilayer;
}

void NewFullShowerGeometryManager::WhereInCalForGeomCU(double *xyz, int *whereincal)
{
  whereincal[0] = 0;
  whereincal[1] = 0;
  whereincal[2] = 0;
  whereincal[3] = 0;
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  if(xyz[0]>2*FSGM_towerPitch || xyz[0]<-FSGM_towerPitch) return;
  if(fabs(xyz[1])>FSGM_towerPitch/2) return;
  whereincal[2] = (int)floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double xmod = xyz[0]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double ymod = xyz[1]+FSGM_towerPitch/2;
  if(ilayer%2==0)
    whereincal[3] = (int)floor((ymod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  else
    whereincal[3] = (int)floor((xmod-(FSGM_towerPitch/2-6*FSGM_cellHorPitch))/FSGM_cellHorPitch);
  if(FSGM_towerPitch-xmod<xmod) xmod = FSGM_towerPitch-xmod;
  if(FSGM_towerPitch-ymod<ymod) ymod = FSGM_towerPitch-ymod;
  if(xmod<FSGM_crackHalfWidth || ymod<FSGM_crackHalfWidth)
    {
      whereincal[0] = 2;
      return;
    }
  double xy = ymod;
  if(ilayer%2==1) xy = xmod;
  xy -= (FSGM_towerPitch/2-6*FSGM_cellHorPitch);
  double xymod = xy-FSGM_cellHorPitch*floor(xy/FSGM_cellHorPitch);
  if(FSGM_cellHorPitch-xymod<xymod) xymod = FSGM_cellHorPitch-xymod;
  if(xymod<(FSGM_cellHorPitch-FSGM_CsIWidth)/2) return;
/*   double zlayer = (FSGM_calZTop-xyz[2])-FSGM_cellVertPitch*(double)ilayer; */
/*   if(FSGM_cellVertPitch-zlayer>zlayer) zlayer = FSGM_cellVertPitch-zlayer; */
/*   if(zlayer<(FSGM_cellVertPitch-FSGM_CsIHeight)/2) return; */
  whereincal[0] = 1;
  whereincal[1] = ilayer;
}

void NewFullShowerGeometryManager::WhereInCalLT(double *xyz, int *whereincal)
{
  whereincal[0] = 0;
  whereincal[1] = 0;
  // tower and column currently not used so not filled for the LT case
  whereincal[2] = 0; // tower
  whereincal[3] = 0; // column

  // test if outside z_dimension
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;

  // test if outside xy_dimensions
  if(fabs(xyz[0])>=2*FSGM_towerPitch) return;
  if(fabs(xyz[1])>=2*FSGM_towerPitch) return;

  // get layer number
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  whereincal[1] = ilayer;

  // retrieve position in look-up table
  int i = (int)((xyz[0]+2*FSGM_towerPitch)/geom_xy_step);
  int j = (int)((xyz[1]+2*FSGM_towerPitch)/geom_xy_step);
  
  if(ilayer%2==0)  // look-up table set for even layers, invert x and y otherwise
    {
      if(geom_x_mat[i][1] || geom_y_mat[j][1])
        whereincal[0] = 2;
      else
        whereincal[0] = (int)(geom_x_mat[i][0] && geom_y_mat[j][0]);
    }
  else
    {
      if(geom_y_mat[i][1] || geom_x_mat[j][1])
        whereincal[0] = 2;
      else
        whereincal[0] = (int)(geom_y_mat[i][0] && geom_x_mat[j][0]);
    }
}

void NewFullShowerGeometryManager::WhereInCalLT2(double *xyz, int *whereincal)
{
  whereincal[0] = 3;
  whereincal[1] = -1;
  whereincal[2] = 0; // tower
  whereincal[3] = 0; // column

  if(fabs(xyz[0])>=2*FSGM_towerPitch) return;
  if(fabs(xyz[1])>=2*FSGM_towerPitch) return;
  if(xyz[2]>=FSGM_calZTop) return;
  if(xyz[2]<FSGM_calZBot) return;

  whereincal[0] = 0;

  // get layer number
  int ilayer = (int)floor( (FSGM_calZTop-xyz[2])/FSGM_cellVertPitch );
  whereincal[1] = ilayer;

  // retrieve position in look-up table
  int i = (int)((xyz[0]+2*FSGM_towerPitch)/geom_xy_step);
  int j = (int)((xyz[1]+2*FSGM_towerPitch)/geom_xy_step);

  // get tower number before x and y switch
  whereincal[2] = 4*geom_y_wic[j][2]+geom_x_wic[i][2];

  int k;
  if(ilayer%2!=0) // switch i and j
    {
      k = i;
      i = j;
      j = k;
    }

  if(geom_x_wic[i][0]==2||geom_y_wic[j][0]==2)
    whereincal[0] = 2;
  else if(geom_x_wic[i][0]==1&&geom_y_wic[j][0]==1)
    whereincal[0] = 1;

  whereincal[3] = geom_y_wic[j][3];
}

void NewFullShowerGeometryManager::WhereInCal(double *xyz, int *whereincal)
{
  if(FSGM_flight_geom)
    {
      //      WhereInCalForGeom(xyz,whereincal); // not needed for saturation handling
      WhereInCalLT2(xyz,whereincal);
    }
  else
    WhereInCalForGeomCU(xyz,whereincal);
}

void NewFullShowerGeometryManager::PrepRadialProfile(double loge)
{
  FSGM_RCORE_0 = FSGM_RCORE_0_0+FSGM_RCORE_0_1*loge;
  FSGM_RCORE_1 = FSGM_RCORE_1_0+FSGM_RCORE_1_1*loge;
  //
  FSGM_RTAIL_0 = FSGM_RTAIL_0_0+FSGM_RTAIL_0_1*loge;
  FSGM_RTAIL_1 = FSGM_RTAIL_1_0+FSGM_RTAIL_1_1*loge;
  FSGM_RTAIL_2 = FSGM_RTAIL_2_0;
  FSGM_RTAIL_3 = FSGM_RTAIL_3_0+FSGM_RTAIL_3_1*loge;
  FSGM_RTAIL_4 = FSGM_RTAIL_4_0;
  FSGM_RTAIL_5 = FSGM_RTAIL_5_0+FSGM_RTAIL_5_1*loge;
  //
  FSGM_PCORE_0 = FSGM_PCORE_0_0+FSGM_PCORE_0_1*loge+FSGM_PCORE_0_2*loge*loge+FSGM_PCORE_0_3*loge*loge*loge;
  FSGM_PCORE_1 = FSGM_PCORE_1_0+FSGM_PCORE_1_1*loge+FSGM_PCORE_1_2*loge*loge;
  FSGM_PCORE_2 = FSGM_PCORE_2_0+FSGM_PCORE_2_1*loge+FSGM_PCORE_2_2*loge*loge+FSGM_PCORE_2_3*loge*loge*loge;
  FSGM_PCORE_3 = FSGM_PCORE_3_0+FSGM_PCORE_3_1*loge+FSGM_PCORE_3_2*loge*loge;
  if(FSGM_PCORE_3<FSGM_PCORE_3_3) FSGM_PCORE_3 = FSGM_PCORE_3_3;
}

double NewFullShowerGeometryManager::RCore(double t)
{
  double tt = t;
  if(tt<FSGM_tmin) tt = FSGM_tmin;
  if(tt>FSGM_tmax) tt = FSGM_tmax;
  return (FSGM_RCORE_0+FSGM_RCORE_1*tt);
}

double NewFullShowerGeometryManager::RTail(double t)
{
  double tt = t;
  if(tt<FSGM_tmin) tt = FSGM_tmin;
  if(tt>FSGM_tmax) tt = FSGM_tmax;
  double myval = FSGM_RTAIL_5;
  if(tt<FSGM_RTAIL_0)
    myval += FSGM_RTAIL_1*pow(fabs(tt-FSGM_RTAIL_0),FSGM_RTAIL_2);
  else
    myval += FSGM_RTAIL_3*pow(tt-FSGM_RTAIL_0,FSGM_RTAIL_4);
  return myval;
}

double NewFullShowerGeometryManager::PCore(double t)
{
  double tt = t;
  if(tt<FSGM_tmin) tt = FSGM_tmin;
  if(tt>FSGM_tmax) tt = FSGM_tmax;

  double arg;
  if(tt>FSGM_PCORE_2)
    {
      arg = FSGM_PCORE_0+FSGM_PCORE_1*tt;
      if(arg<0) arg = 0;
      if(arg>1) arg = 1;
      return arg;
    }

  double x0 = FSGM_PCORE_2;
  double g1 = FSGM_PCORE_3;
  double b1 = FSGM_PCORE_1-2*g1*x0;
  double a1 = FSGM_PCORE_0+FSGM_PCORE_1*x0-b1*x0-g1*x0*x0;
  arg = a1+b1*tt+g1*tt*tt;
  if(arg<0) arg = 0;
  if(arg>1) arg = 1;
  return arg;
}

double NewFullShowerGeometryManager::RadialProfile(double t, double r)
{
  double rcore2 = RCore(t);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(t);
  rtail2 = rtail2*rtail2;
  double pcore = PCore(t);
  return pcore*(2*rcore2/(r*r+rcore2)/(r*r+rcore2)) + (1-pcore)*(2*rtail2/(r*r+rtail2)/(r*r+rtail2));
}

double NewFullShowerGeometryManager::RadialProfileMod(double t, double r, double pp, double rtailmod)
{
  double rcore2 = RCore(t);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(t)*rtailmod;
  rtail2 = rtail2*rtail2;
  double pcore = PCore(t)-pp;
  if(pcore<0) pcore = 0;
  if(pcore>1) pcore = 1;
  return pcore*(2*rcore2/(r*r+rcore2)/(r*r+rcore2)) + (1-pcore)*(2*rtail2/(r*r+rtail2)/(r*r+rtail2));
}

double NewFullShowerGeometryManager::GetEffectiveRadius(double t, double radialcontainedfraction)
{
  double rcore2 = RCore(t);
  rcore2 = rcore2*rcore2;
  double rtail2 = RTail(t);
  rtail2 = rtail2*rtail2;
  double pcore = PCore(t);
  
  double a = (1-radialcontainedfraction)/rcore2/rtail2;
  double b = (1-radialcontainedfraction)*(1/rcore2+1/rtail2) - pcore/rtail2 - (1-pcore)/rcore2;
  double c = -radialcontainedfraction;
  double delta = b*b-4*a*c;
  double r1 = (-b+sqrt(delta))/2/a;
  r1 = sqrt(r1);

  if(r1>RadProf_r_max) r1 = RadProf_r_max;

  return r1;
}

void NewFullShowerGeometryManager::FillRadialProfile(int numtab)
{
  // Defining limits and steps
  RadProf_r_max = 5.;
  RadProf_t_max = 3.;  
  RadProf_r_step = RadProf_r_max/(double)FSGM_RPROF_R_MAX;
  RadProf_t_step = RadProf_t_max/(double)FSGM_RPROF_T_MAX;

  int i,j;
  double r,t;
  for(i=0;i<FSGM_RPROF_T_MAX;++i)
    {
      t = RadProf_t_step*(0.5+(double)i);
      for(j=0;j<FSGM_RPROF_R_MAX;++j)
        {
          r = RadProf_r_step*(0.5+(double)j);
          if(numtab==0)
            RadProf[i][j] = RadialProfile(t,r);
          else
            RadProf2[i][j] = RadialProfile(t,r);
        }
    }
}

void NewFullShowerGeometryManager::FillRadialProfile(double loge, int numtab)
{
  PrepRadialProfile(loge);
  FillRadialProfile(numtab);
}

double NewFullShowerGeometryManager::GetRadialProfile(double t, double r)
{
  int i_t = (int)(t/RadProf_t_step);
  if(i_t>=FSGM_RPROF_T_MAX)
    i_t = FSGM_RPROF_T_MAX-1;

  int i_r = (int)(r/RadProf_r_step);
  if(i_r>=FSGM_RPROF_R_MAX)
    i_r = FSGM_RPROF_R_MAX-1;

  if(opthe==0) return RadProf[i_t][i_r];
  return RadProf2[i_t][i_r];
}

double NewFullShowerGeometryManager::GetRadialProfile2(double t_in, double r_in)
{
  double t = t_in;
  if(t>RadProf_t_max) t = RadProf_t_max;
  double r = r_in;
  if(r>RadProf_r_max) r = RadProf_r_max;
  return RadialProfile(t,r);
}

double NewFullShowerGeometryManager::GetRadialProfileMod2(double t_in, double r_in, double pp, double rtailmod)
{
  double t = t_in;
  if(t>RadProf_t_max) t = RadProf_t_max;
  double r = r_in;
  if(r>RadProf_r_max) r = RadProf_r_max;
  return RadialProfileMod(t,r,pp,rtailmod);
}

double NewFullShowerGeometryManager::GetZCrack(double *xyz, double *pp, double *vv)
{
  double zcrack = 0;
  if(fabs(xyz[0])>=2*FSGM_towerPitch) return zcrack;
  if(fabs(xyz[1])>=2*FSGM_towerPitch) return zcrack;
  if(xyz[2]>=FSGM_calZTop) return zcrack;
  if(xyz[2]<FSGM_calZBot) return zcrack;

  double xmod = xyz[0]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[0]+2*FSGM_towerPitch)/FSGM_towerPitch);
  double ymod = xyz[1]+2*FSGM_towerPitch-FSGM_towerPitch*floor((xyz[1]+2*FSGM_towerPitch)/FSGM_towerPitch);

  xmod -= FSGM_towerPitch/2.;
  ymod -= FSGM_towerPitch/2.;

  double vvp[3];
  vvp[0] = 1;
  vvp[1] = 0;
  vvp[2] = 0;

  if(fabs(ymod)>fabs(xmod))
    {
      vvp[0] = 0;
      vvp[1] = 1;
      vvp[2] = 0;
    }

  double crossprod = vvp[0]*vv[0]+vvp[1]*vv[1]+vvp[2]*vv[2];
  if(crossprod==0) return zcrack;

  double lambda = (vvp[0]*(xyz[0]-pp[0])+vvp[1]*(xyz[1]-pp[1])+vvp[2]*(xyz[2]-pp[2]))/crossprod;

  double ppp[3];

  int i;
  for(i=0;i<3;++i) ppp[i] = pp[i]+lambda*vv[i];

  zcrack = ppp[2];
  if(zcrack>0) zcrack = 0;
  if(zcrack<FSGM_calZBot) zcrack = FSGM_calZBot;

  return zcrack;
}


/**   
* @class NewFullShowerDevelopmentDescription
* @author Philippe Bruel
*
* Tool that describes the shower developement in the calorimeter given
* the length in X0 seen in the tracker and the position of the shower maximum
*
* $Header$
*/

NewFullShowerDevelopmentDescription::NewFullShowerDevelopmentDescription(NewFullShowerGeometryManager *fsgm_input, int type_input, double zstep_input, double radialcontainedfraction_input)
  :m_fsgm(fsgm_input)
{
  Initialize();
  Reset();
  Type = type_input;
  ZStep = zstep_input;
  ZStepRef = zstep_input;
  RadialContainedFraction = radialcontainedfraction_input;
  wideningfactor = 1.;
}

NewFullShowerDevelopmentDescription::~NewFullShowerDevelopmentDescription()
{

}

void NewFullShowerDevelopmentDescription::Initialize()
{

  // Radiation lengths
  FSDD_MOLRAD = 35.;
  FSDD_gCSI = 4.51;
  FSDD_gCRK = 0.607287449392712508; // 10.*2.7/(2.*22.23) // 10mm
  FSDD_XCSI = 18.5*(m_fsgm->FSGM_cellVertPitch/m_fsgm->FSGM_CsIHeight);
  FSDD_XCRK = 395.364666666666722; // 24.01/FSDD_gCRK*10. // *10 to get it in mm

  Reset();
}

void NewFullShowerDevelopmentDescription::Reset()
{

  int i,j,k;
  NStep = 0;
  ZStep = 0;
  ZStepRef = 0;
  x0maxshower = 0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  startx0 = 0;
  lastx0 = 0;
  totx0cal = 0;
  totx0crack = 0;
  for(i=0;i<8;++i) totx0lay[i] = 0;
  for(i=0;i<8;++i) posx0lay[i] = 0;
  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      dX0[i] = 0;
      X0[i] = 0;
      RM[i] = 0;
      for(j=0;j<4;++j) materialfraction[j][i] = 0;
      for(j=0;j<8;++j) layerfraction[j][i] = 0;
      for(j=0;j<FSDD_XTAL_NMAX;++j) xtalfraction[j][i] = 0;
    }
  NXtal = 0;
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        OffSatu[i][j][k] = -1;
}

void NewFullShowerDevelopmentDescription::GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend)
{
  int i,j;
  double lambda;
  double pp0[3][3];
  double pp1[3][3];
  double segmentlength[3];
  // looking for the best axis
  int idir = -1;
  double segmentlengthmin = 99999999;

  for(i=0;i<3;++i)
    {
      if(vv[i]==0) continue;
      lambda = (m_fsgm->FSGM_CalSafeBoundaries[i][0]-pp[i])/vv[i];
      for(j=0;j<3;++j)
        {
          if(j!=i) pp0[i][j] = pp[j]+lambda*vv[j];
          else pp0[i][j] = m_fsgm->FSGM_CalSafeBoundaries[i][0];
        }
      lambda = (m_fsgm->FSGM_CalSafeBoundaries[i][1]-pp[i])/vv[i];
      for(j=0;j<3;++j)
        {
          if(j!=i) pp1[i][j] = pp[j]+lambda*vv[j];
          else pp1[i][j] = m_fsgm->FSGM_CalSafeBoundaries[i][1];
        }
      segmentlength[i] = 0;
      for(j=0;j<3;++j) segmentlength[i] += (pp0[i][j]-pp1[i][j])*(pp0[i][j]-pp1[i][j]);
      segmentlength[i] = sqrt(segmentlength[i]);
      if(segmentlength[i]<segmentlengthmin)
        {
          segmentlengthmin = segmentlength[i];
          idir = i;
        }
    }

  double checkdirection = (pp1[idir][0]-pp0[idir][0])*vv[0]+(pp1[idir][1]-pp0[idir][1])*vv[1]+(pp1[idir][2]-pp0[idir][2])*vv[2];
  
  if(checkdirection>0)
    {
      for(i=0;i<3;++i)
        {
          ppstart[i] = pp0[idir][i];
          ppend[i] = pp1[idir][i];
        }
    }
  else
    {
      for(i=0;i<3;++i)
        {
          ppstart[i] = pp1[idir][i];
          ppend[i] = pp0[idir][i];
        }
    }
}

bool NewFullShowerDevelopmentDescription::Compute(double *pp, double *vv, double startx0_input, double x0maxshower_input, double zstep_input, double totrlnmax)
{
  ZStepRef = zstep_input;

  if(Type!=0)
    {
      //  NewFullShowerDevelopmentDescription WRONG TYPE
      return false;
    }

  int i,j;
  int whereincal[4];
  // whereincal[0] = material
  // whereincal[1] = layer
  // whereincal[2] = tower
  // whereincal[3] = column

  NStep = 0;
  x0maxshower = x0maxshower_input;
  startx0 = startx0_input;
  lastx0 = 0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  if(x0maxshower<2)
    {
      // NewFullShowerDevelopmentDescription x0maxshower is too small
      x0maxshower = 2;
    }

  double ppc[3];
  double ppcc[3];
  double pp0[3];
  double pp1[3];
  double vv0[3];
  double vv1[3];
  double vv2[3];
  double vnorm = sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
  if(vnorm==0) return false;

  vv2[0] = vv[0]/vnorm;
  vv2[1] = vv[1]/vnorm;
  vv2[2] = vv[2]/vnorm;

  GetTrajectorySegment(pp,vv2,pp0,pp1);

  vv0[0] = vv2[1];
  vv0[1] = -vv2[0];
  vv0[2] = 0;
  vnorm = sqrt(vv0[0]*vv0[0]+vv0[1]*vv0[1]+vv0[2]*vv0[2]);
  if(vnorm>0)
    {
      vv0[0] = vv0[0]/vnorm;
      vv0[1] = vv0[1]/vnorm;
      vv0[2] = vv0[2]/vnorm;
    }
  else
    {
      vv0[0] = 1;
      vv0[1] = 0;
      vv0[2] = 0;
    }
  vv1[0] = vv2[1]*vv0[2]-vv2[2]*vv0[1];
  vv1[1] = vv2[2]*vv0[0]-vv2[0]*vv0[2];
  vv1[2] = vv2[0]*vv0[1]-vv2[1]*vv0[0];

  double lambda = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );

  int nhist = (int)floor(lambda/ZStepRef);
  if(nhist>=FSDD_NSTEPS_MAX)
    {
      // NewFullShowerDevelopmentDescription::Compute : nhist >= FSDD_NSTEPS_MAX
      return false;
    }
  NStep = nhist;
  ZStep = ZStepRef;

  double ppref[3];
  ppref[2] = 0;
  double lambdaref = -100;
  double zreftmp;
  for(i=0;i<=8;++i)
    {
      zreftmp = m_fsgm->FSGM_calZTop-m_fsgm->FSGM_cellVertPitch*(double)i;
      lambdaref = -100;
      if(vv2[2]!=0)
        lambdaref = (zreftmp-pp0[2])/vv2[2];
      if(lambdaref>0)
        {
          for(j=0;j<3;++j) ppref[j] = pp0[j]+lambdaref*vv2[j];
          break;
        }
    }
  double myzstep;
  int mynstep;
  double lengthref;

  if(ppref[2]<0 && vv2[2]<-0.001)
    {
      lengthref = fabs(m_fsgm->FSGM_cellVertPitch/vv2[2]);
      mynstep = (int)floor(lengthref/ZStepRef)+1;
      ZStep = lengthref/(double)mynstep;
      // change pp0
      lambdaref = sqrt( (ppref[0]-pp0[0])*(ppref[0]-pp0[0])+(ppref[1]-pp0[1])*(ppref[1]-pp0[1])+(ppref[2]-pp0[2])*(ppref[2]-pp0[2]) );
      mynstep = (int)floor(lambdaref/ZStep)+1;
      lambdaref = ZStep*(double)mynstep;
      for(j=0;j<3;++j) pp0[j] = ppref[j]-lambdaref*vv2[j];
      // change pp1
      lambdaref = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );
      mynstep = (int)floor(lambdaref/ZStep)+1;
      lambdaref = ZStep*(double)mynstep;
      for(j=0;j<3;++j) pp1[j] = pp0[j]+lambdaref*vv2[j];

      lambda = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );
      
      nhist = (int)floor(lambda/ZStep);
      if(nhist>=FSDD_NSTEPS_MAX)
        {
          // NewMultiFullShowerDevelopmentDescription::Compute : nhist >= FSDD_NSTEPS_MAX
          return false;
        }
      NStep = nhist;
    }

  double etotdep;

  double z2X0mat;

  X0[0] = startx0_input;
  lastx0 = startx0_input;

  double radius, relradius;
  double effradius, releffradius;
  double radialprofile;

  double x0position = X0[0];
  double relx0position;

  totx0cal = 0;
  totx0crack = 0;
  for(j=0;j<8;++j) totx0lay[j] = 0;
  for(j=0;j<8;++j) posx0lay[j] = 0;

  int isatxtal;

  bool InCal = false;
  bool OutofCal = false;

  effradius = FSDD_MOLRAD;

  etotdep = 0;
  for(i=0;i<m_fsgm->FSGM_NCircle;++i)
    {
      m_fsgm->FSGM_InCircle[i] = 1;
      //
      m_fsgm->FSGM_ECircle[i] = m_fsgm->GetRadialProfile(1,effradius*m_fsgm->FSGM_RCircle[i]/FSDD_MOLRAD/wideningfactor)*m_fsgm->FSGM_WCircle[i];
      etotdep += m_fsgm->FSGM_ECircle[i];
      //
      m_fsgm->FSGM_TCircle[i] = X0[0];
      m_fsgm->FSGM_TPedCircle[i] = X0[0];
      m_fsgm->FSGM_RadCircle[i] = effradius*m_fsgm->FSGM_RCircle[i]*exp((m_fsgm->FSGM_TCircle[i]-m_fsgm->FSGM_TPedCircle[i])/x0maxshower-1);
      //
      // proj onto cal top
      //
      ppcc[0] = pp0[0]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[0]+m_fsgm->FSGM_SPCircle[i]*vv1[0]);
      ppcc[1] = pp0[1]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[1]+m_fsgm->FSGM_SPCircle[i]*vv1[1]);
      ppcc[2] = pp0[2]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[2]+m_fsgm->FSGM_SPCircle[i]*vv1[2]);
      m_fsgm->FSGM_LambdaCircle[i] = (m_fsgm->FSGM_calZTop-ppcc[2])/vv2[2]-ZStep/2;
      ppcc[0] = pp0[0]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[0]+m_fsgm->FSGM_SPCircle[i]*vv1[0])+m_fsgm->FSGM_LambdaCircle[i]*vv2[0];
      ppcc[1] = pp0[1]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[1]+m_fsgm->FSGM_SPCircle[i]*vv1[1])+m_fsgm->FSGM_LambdaCircle[i]*vv2[1];
      ppcc[2] = pp0[2]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[2]+m_fsgm->FSGM_SPCircle[i]*vv1[2])+m_fsgm->FSGM_LambdaCircle[i]*vv2[2];
    }

  // normalize ECircle to 1
  if(etotdep>0)
    {
      for(i=0;i<m_fsgm->FSGM_NCircle;++i)
        {
          m_fsgm->FSGM_ECircle[i] /= etotdep;
        }
    }

  int nstepincrackmax = 1000;
  nstepincrackmax = (int)floor(FSDD_XCRK/FSDD_XCSI);

  double lengthincrackmean;
  double lengthincrackmeanw;
  int nincrack = 0;
  for(j=0;j<NStep;++j)
    {
      if(j>0)  X0[j] = dX0[j-1] + X0[j-1];
      //
      for(i=0;i<4;++i) materialfraction[i][j] = 0;
      for(i=0;i<8;++i) layerfraction[i][j] = 0;
      for(i=0;i<NXtal;++i) xtalfraction[i][j] = 0;
      //
      lengthincrackmean = 0.;
      lengthincrackmeanw = 0.;
      //
      for(i=0;i<m_fsgm->FSGM_NCircle;++i)
        {
          if(m_fsgm->FSGM_InCircle[i])
            {
              m_fsgm->FSGM_LambdaCircle[i] += ZStep;
              m_fsgm->FSGM_TCircle[i] += ZStep/FSDD_XCSI;
              m_fsgm->FSGM_RadCircle[i] = effradius*m_fsgm->FSGM_RCircle[i]*exp((m_fsgm->FSGM_TCircle[i]-m_fsgm->FSGM_TPedCircle[i])/x0maxshower-1);
              ppcc[0] = pp0[0]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[0]+m_fsgm->FSGM_SPCircle[i]*vv1[0])+m_fsgm->FSGM_LambdaCircle[i]*vv2[0];
              ppcc[1] = pp0[1]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[1]+m_fsgm->FSGM_SPCircle[i]*vv1[1])+m_fsgm->FSGM_LambdaCircle[i]*vv2[1];
              ppcc[2] = pp0[2]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[2]+m_fsgm->FSGM_SPCircle[i]*vv1[2])+m_fsgm->FSGM_LambdaCircle[i]*vv2[2];
              m_fsgm->WhereInCal(ppcc,whereincal);
              //
              //
              if(whereincal[0]==3)
                {
                  m_fsgm->FSGM_InCircle[i] = 0;
                }
              else if(whereincal[0]==2)
                {
                  nincrack = 0;
                  while(nincrack<nstepincrackmax)
                    {
                      m_fsgm->FSGM_LambdaCircle[i] += ZStep;
                      //                      m_fsgm->FSGM_TCircle[i] += ZStep/FSDD_XCSI;
                      m_fsgm->FSGM_TCircle[i] += ZStep/FSDD_XCSI/2.;
                      m_fsgm->FSGM_RadCircle[i] = effradius*m_fsgm->FSGM_RCircle[i]*exp((m_fsgm->FSGM_TCircle[i]-m_fsgm->FSGM_TPedCircle[i])/x0maxshower-1);
                      ppcc[0] = pp0[0]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[0]+m_fsgm->FSGM_SPCircle[i]*vv1[0])+m_fsgm->FSGM_LambdaCircle[i]*vv2[0];
                      ppcc[1] = pp0[1]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[1]+m_fsgm->FSGM_SPCircle[i]*vv1[1])+m_fsgm->FSGM_LambdaCircle[i]*vv2[1];
                      ppcc[2] = pp0[2]+m_fsgm->FSGM_RadCircle[i]*(m_fsgm->FSGM_CPCircle[i]*vv0[2]+m_fsgm->FSGM_SPCircle[i]*vv1[2])+m_fsgm->FSGM_LambdaCircle[i]*vv2[2];
                      m_fsgm->WhereInCal(ppcc,whereincal);
                      if(whereincal[0]!=2)
                        {
                          m_fsgm->FSGM_LambdaCircle[i] -= ZStep;
                          m_fsgm->FSGM_TCircle[i] -= ZStep/FSDD_XCSI/2.;
                          whereincal[0] = 2;
                          break;
                        }
                      ++nincrack;
                    }
                  lengthincrackmean += m_fsgm->FSGM_ECircle[i]*(double)nincrack;
                  lengthincrackmeanw += m_fsgm->FSGM_ECircle[i];
                }
            }
          else
            whereincal[0] = 3;

          radialprofile = m_fsgm->FSGM_ECircle[i];
          etotdep += radialprofile;
          materialfraction[whereincal[0]][j] += radialprofile;
          if(whereincal[0]==1)
            {
              layerfraction[whereincal[1]][j] += radialprofile;
              //
              isatxtal = OffSatu[whereincal[2]][whereincal[1]][whereincal[3]];
              if(isatxtal>=0) xtalfraction[isatxtal][j] += radialprofile;
              //              if(isatxtal>=0) xtalpoint[isatxtal] += 1;
            }
        }

      // if(etotdep>0)
      //        {
      //          for(i=0;i<4;++i)
      //            materialfraction[i][j] /= etotdep;
      //          for(i=0;i<8;++i)
      //            layerfraction[i][j] /= etotdep;
      //          for(i=0;i<NXtal;++i)
      //            xtalfraction[i][j] /= etotdep;
      //        }

      //      if(lengthincrackmeanw>0) lengthincrackmean /= lengthincrackmeanw;

      // Ignore the FakeCsI volume contribution when computing lastx0
      z2X0mat = materialfraction[1][j]/FSDD_XCSI + lengthincrackmean/FSDD_XCRK;
      //      z2X0mat = materialfraction[1][j]/FSDD_XCSI + materialfraction[2][j]*lengthincrackmean/FSDD_XCRK;
      lastx0 += ZStep*z2X0mat;
      // But take into account for the shower development
      if(OutofCal) z2X0mat += materialfraction[3][j]/FSDD_XCSI;
        
      dX0[j] = ZStep*z2X0mat;
      x0position = X0[j];
      totx0cal += materialfraction[1][j]*dX0[j];
      totx0crack += materialfraction[2][j]*dX0[j];
      for(i=0;i<8;++i) totx0lay[i] += layerfraction[i][j]*dX0[j];

      if(materialfraction[2][j]>crackmaxfrac)
        {
          x0crackmaxfrac = X0[j];
          crackmaxfrac = materialfraction[2][j];
        }
    }

  dX0[NStep] = 0;
  X0[NStep] = dX0[NStep-1] + X0[NStep-1];
  materialfraction[0][NStep] = 1;
  for(i=1;i<4;++i) materialfraction[i][NStep] = 0;
  for(i=0;i<8;++i) layerfraction[i][NStep] = 0;
  for(i=0;i<NXtal;++i) xtalfraction[i][NStep] = 0;

  RemoveEmptySteps();

  return true;
}

void NewFullShowerDevelopmentDescription::RemoveEmptySteps()
{
  int i,j;
  int iwo0 = 0;
  for(i=0;i<NStep;++i)
    {
      if(dX0[i]==0) continue;
      if(iwo0==i) continue;
      X0[iwo0] = X0[i];
      dX0[iwo0] = dX0[i];
      RM[iwo0] = RM[i];
      for(j=0;j<4;++j) materialfraction[j][iwo0] = materialfraction[j][i];
      for(j=0;j<8;++j) layerfraction[j][iwo0] = layerfraction[j][i];
      for(j=0;j<NXtal;++j) xtalfraction[j][iwo0] = xtalfraction[j][i];
      ++iwo0;
    }
  dX0[iwo0] = 0;
  X0[iwo0] = X0[NStep];
  for(j=0;j<4;++j) materialfraction[j][iwo0] = materialfraction[j][NStep];
  for(j=0;j<8;++j) layerfraction[j][iwo0] = layerfraction[j][NStep];
  for(j=0;j<NXtal;++j) xtalfraction[j][iwo0] = xtalfraction[j][NStep];
  NStep = iwo0;
}

bool NewFullShowerDevelopmentDescription::ConvertToFixedX0(double x0step, NewFullShowerDevelopmentDescription *shmm)
{
  if(shmm==NULL) return false;
  if(shmm->NStep==0) return false;

  int i,j,ix0;
  double check[FSDD_NSTEPS_MAX];
  int icheck[FSDD_NSTEPS_MAX];

  NStep = 0;
  x0maxshower = shmm->x0maxshower;
  RadialContainedFraction = shmm->RadialContainedFraction ;
  startx0 = shmm->startx0;
  lastx0 = shmm->lastx0;
  crackmaxfrac = 0;
  x0crackmaxfrac = 0;
  totx0cal = shmm->totx0cal;
  totx0crack = 0;
  for(i=0;i<8;++i) totx0lay[i] = shmm->totx0lay[i];
  for(i=0;i<8;++i) posx0lay[i] = -1;

  ZStep = -1;

  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      check[i] = 0;
      icheck[i] = 0;
      X0[i] = x0step*(double)i;
      dX0[i] = x0step;
      for(j=0;j<4;++j) materialfraction[j][i] = 0;
      for(j=0;j<8;++j) layerfraction[j][i] = 0;  
      for(j=0;j<NXtal;++j) xtalfraction[j][i] = 0;  
      RM[i] = m_fsgm->GetEffectiveRadius(X0[i]/x0maxshower,RadialContainedFraction);
    }

  NStep = (int)floor(shmm->X0[shmm->NStep]/x0step)+1;
  if(NStep>=FSDD_NSTEPS_MAX)
    {
      //  NewFullShowerDevelopmentDescription::Convert : NStep >= FSDD_NSTEPS_MAX
      return false;
    }

  int start_ix0 = (int)floor(shmm->X0[0]/x0step);
  for(i=0;i<start_ix0;++i)
    materialfraction[0][i] = 1;

  int dimm = 0;
  int lastdimm = 0;

  double last_materialfraction[4];
  double last_layerfraction[8];
  double last_xtalfraction[FSDD_XTAL_NMAX];

  last_materialfraction[0] = 1;
  for(j=1;j<4;++j) last_materialfraction[j] = 0;
  for(j=0;j<8;++j) last_layerfraction[j] = 0;
  for(j=0;j<NXtal;++j) last_xtalfraction[j] = 0;

  int imm = 0;
  for(ix0=start_ix0;ix0<NStep;++ix0)
    {
      lastdimm = dimm;
      dimm = 0;
      while(shmm->X0[imm]<X0[ix0+1] && imm<shmm->NStep)
        {
          if(dimm==0)
            {
              check[ix0] += shmm->X0[imm]-X0[ix0];
              for(j=0;j<4;++j) materialfraction[j][ix0] += (shmm->X0[imm]-X0[ix0])*last_materialfraction[j];
              for(j=0;j<8;++j) layerfraction[j][ix0] += (shmm->X0[imm]-X0[ix0])*last_layerfraction[j];
              for(j=0;j<NXtal;++j) xtalfraction[j][ix0] += (shmm->X0[imm]-X0[ix0])*last_xtalfraction[j];
            }
          icheck[ix0] += 1;
          if(shmm->X0[imm+1]<X0[ix0+1])
            {
              check[ix0] += shmm->dX0[imm];
              for(j=0;j<4;++j) materialfraction[j][ix0] += shmm->dX0[imm]*shmm->materialfraction[j][imm];
              for(j=0;j<8;++j) layerfraction[j][ix0] += shmm->dX0[imm]*shmm->layerfraction[j][imm];
              for(j=0;j<NXtal;++j) xtalfraction[j][ix0] += shmm->dX0[imm]*shmm->xtalfraction[j][imm];
            }
          else
            {
              check[ix0] += (X0[ix0+1]-shmm->X0[imm]);
              for(j=0;j<4;++j) materialfraction[j][ix0] += (X0[ix0+1]-shmm->X0[imm])*shmm->materialfraction[j][imm];
              for(j=0;j<8;++j) layerfraction[j][ix0] += (X0[ix0+1]-shmm->X0[imm])*shmm->layerfraction[j][imm];
              for(j=0;j<NXtal;++j) xtalfraction[j][ix0] += (X0[ix0+1]-shmm->X0[imm])*shmm->xtalfraction[j][imm];
            }
          for(j=0;j<4;++j) last_materialfraction[j] = shmm->materialfraction[j][imm];
          for(j=0;j<8;++j) last_layerfraction[j] = shmm->layerfraction[j][imm];
          for(j=0;j<NXtal;++j) last_xtalfraction[j] = shmm->xtalfraction[j][imm];
          ++imm;
          ++dimm;
        }
      if(dimm==0 && ix0<NStep-1)
        {
          for(j=0;j<4;++j) materialfraction[j][ix0] = last_materialfraction[j];
          for(j=0;j<8;++j) layerfraction[j][ix0] = last_layerfraction[j];
          for(j=0;j<NXtal;++j) xtalfraction[j][ix0] = last_xtalfraction[j];
        }
//       if(dimm==0)
//         printf("BRUEL ix0 %d X0 %f dimm %d check %f frac %f %f %f %f %f %f %f %f\n",ix0,X0[ix0],dimm,check[ix0],
//                layerfraction[0][ix0],layerfraction[1][ix0],layerfraction[2][ix0],layerfraction[3][ix0],layerfraction[4][ix0],layerfraction[5][ix0],layerfraction[6][ix0],layerfraction[7][ix0]);
//       else
//         printf("BRUEL ix0 %d X0 %f dimm %d check %f frac %f %f %f %f %f %f %f %f\n",ix0,X0[ix0],dimm,check[ix0],
//                layerfraction[0][ix0]/x0step,layerfraction[1][ix0]/x0step,layerfraction[2][ix0]/x0step,layerfraction[3][ix0]/x0step,layerfraction[4][ix0]/x0step,layerfraction[5][ix0]/x0step,layerfraction[6][ix0]/x0step,layerfraction[7][ix0]/x0step);
    }

  ix0 = NStep-1;
  if(icheck[ix0]==0)
    {
      check[ix0] = shmm->X0[shmm->NStep]-X0[ix0];
      for(j=0;j<4;++j) materialfraction[j][ix0] = (shmm->X0[shmm->NStep]-X0[ix0])*last_materialfraction[j];
      for(j=0;j<8;++j) layerfraction[j][ix0] = (shmm->X0[shmm->NStep]-X0[ix0])*last_layerfraction[j];
      for(j=0;j<NXtal;++j) xtalfraction[j][ix0] = (shmm->X0[shmm->NStep]-X0[ix0])*last_xtalfraction[j];
    }
  check[ix0] += (X0[NStep]-shmm->X0[shmm->NStep]);
  for(j=0;j<4;++j) materialfraction[j][ix0] += (X0[NStep]-shmm->X0[shmm->NStep])*shmm->materialfraction[j][shmm->NStep];
  for(j=0;j<8;++j) layerfraction[j][ix0] += (X0[NStep]-shmm->X0[shmm->NStep])*shmm->layerfraction[j][shmm->NStep];
  for(j=0;j<NXtal;++j) xtalfraction[j][ix0] += (X0[NStep]-shmm->X0[shmm->NStep])*shmm->xtalfraction[j][shmm->NStep];
  
  for(i=NStep;i<FSDD_NSTEPS_MAX;++i)
    materialfraction[0][i] = 1;

  for(ix0=start_ix0;ix0<NStep;++ix0)
    {
      if(check[ix0]==0) continue;
      if(fabs(check[ix0]-x0step)>0.0001)
        {
          // NewFullShowerDevelopmentDescription::Convert PROBLEM DURING CONVERSION
          return false;
        }
      for(j=0;j<4;++j) materialfraction[j][ix0] /= x0step;
      for(j=0;j<8;++j) layerfraction[j][ix0] /= x0step;
      for(j=0;j<NXtal;++j) xtalfraction[j][ix0] /= x0step;
    }

  double meanpos,meanweight;
  for(j=0;j<8;++j)
    {
      meanpos = 0;
      meanweight = 0;
      for(ix0=0;ix0<NStep;++ix0)
        {
          if(layerfraction[j][ix0]>0)
            {
              meanpos += layerfraction[j][ix0]*X0[ix0];
              meanweight += layerfraction[j][ix0];
            }
        }
      if(meanweight==0)
        posx0lay[j] = -1;
      else
        posx0lay[j] = meanpos/meanweight;
    }

  return true;
}

void NewFullShowerDevelopmentDescription::SetWideningFactor(double widfact)
{
  wideningfactor = widfact;
}

void NewFullShowerDevelopmentDescription::FillFromMultiFullShowerDevelopmentDescription(NewMultiFullShowerDevelopmentDescription *mfsddmm, int index)
{
  int i,j,k;
  Type = 0;
  NStep = mfsddmm->NStep;
  ZStep = mfsddmm->ZStep;
  ZStepRef = mfsddmm->ZStepRef;
  x0maxshower = mfsddmm->x0maxshower[index];
  crackmaxfrac = mfsddmm->crackmaxfrac[index];
  x0crackmaxfrac = mfsddmm->x0crackmaxfrac[index];
  startx0 = mfsddmm->startx0[index];
  lastx0 = mfsddmm->lastx0[index];
  totx0cal = mfsddmm->totx0cal[index];
  totx0crack = mfsddmm->totx0crack[index];
  for(i=0;i<8;++i) totx0lay[i] = mfsddmm->totx0lay[index][i];
  for(i=0;i<8;++i) posx0lay[i] = mfsddmm->posx0lay[index][i];
  for(i=0;i<FSDD_NSTEPS_MAX;++i)
    {
      dX0[i] = mfsddmm->dX0[index][i];
      X0[i] = mfsddmm->X0[index][i];
      RM[i] = mfsddmm->RM[index][i];
      for(j=0;j<4;++j) materialfraction[j][i] = mfsddmm->materialfraction[index][j][i];
      for(j=0;j<8;++j) layerfraction[j][i] = mfsddmm->layerfraction[index][j][i];
      for(j=0;j<NXtal;++j) xtalfraction[j][i] = mfsddmm->xtalfraction[index][j][i];
    }
  NXtal = mfsddmm->NXtal;
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        OffSatu[i][j][k] = mfsddmm->OffSatu[i][j][k];
}

////////////////////////////////////////////////////////////////////////////////////////////////////

NewMultiFullShowerDevelopmentDescription::NewMultiFullShowerDevelopmentDescription(NewFullShowerGeometryManager *fsgm_input, int nxmax, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input)
  :m_fsgm(fsgm_input)
{
  Initialize();
  Reset();

  NDevelopment = nxmax;
  DXMax = dxmax;
  ZStep = zstep_input;
  ZStepRef = zstep_input;
  
  RadialContainedFraction = radialcontainedfraction_input;
  wideningfactor = 1.;
  crackextinction = 2.0;

  int i;
  for(i=0;i<=NDevelopment;++i)
    XMax[i] = xmax0+DXMax*(double)i;
}

NewMultiFullShowerDevelopmentDescription::~NewMultiFullShowerDevelopmentDescription()
{

}

void NewMultiFullShowerDevelopmentDescription::Initialize()
{

  // Radiation lengths
  FSDD_MOLRAD = 35.;
  FSDD_gCSI = 4.51;
  FSDD_gCRK = 0.607287449392712508; // 10.*2.7/(2.*22.23) // 10mm
  FSDD_XCSI = 18.5*(m_fsgm->FSGM_cellVertPitch/m_fsgm->FSGM_CsIHeight);
  FSDD_XCRK = 395.364666666666722; // 24.01/FSDD_gCRK*10. // *10 to get it in mm

  Reset();
}

void NewMultiFullShowerDevelopmentDescription::Reset()
{

  int i,j,k,ii;
  NStep = 0;
  ZStep = 0;
  ZStepRef = 0;
  for(ii=0;ii<FSDD_NMAX;++ii)
    {
      x0maxshower[ii] = 0;
      crackmaxfrac[ii] = 0;
      x0crackmaxfrac[ii] = 0;
      startx0[ii] = 0;
      lastx0[ii] = 0;
      totx0cal[ii] = 0;
      totx0crack[ii] = 0;
    }

  for(ii=0;ii<FSDD_NMAX;++ii)
    for(i=0;i<8;++i)
      totx0lay[ii][i] = 0;

  for(ii=0;ii<FSDD_NMAX;++ii)
    for(i=0;i<8;++i)
      posx0lay[ii][i] = 0;

  for(ii=0;ii<FSDD_NMAX;++ii)
    for(i=0;i<FSDD_NSTEPS_MAX;++i)
      {
        dX0[ii][i] = 0;
        X0[ii][i] = 0;
        RM[ii][i] = 0;
        for(j=0;j<4;++j) materialfraction[ii][j][i] = 0;
        for(j=0;j<8;++j) layerfraction[ii][j][i] = 0;
        for(j=0;j<FSDD_XTAL_NMAX;++j) xtalfraction[ii][j][i] = 0;
      }

  NXtal = 0;
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        OffSatu[i][j][k] = -1;
}

void NewMultiFullShowerDevelopmentDescription::GetTrajectorySegment(double *pp, double *vv, double *ppstart, double *ppend)
{
  int i,j;
  double lambda;
  double pp0[3][3];
  double pp1[3][3];
  double segmentlength[3];
  // looking for the best axis
  int idir = -1;
  double segmentlengthmin = 99999999;

  for(i=0;i<3;++i)
    {
      if(vv[i]==0) continue;
      lambda = (m_fsgm->FSGM_CalSafeBoundaries[i][0]-pp[i])/vv[i];
      for(j=0;j<3;++j)
        {
          if(j!=i) pp0[i][j] = pp[j]+lambda*vv[j];
          else pp0[i][j] = m_fsgm->FSGM_CalSafeBoundaries[i][0];
        }
      lambda = (m_fsgm->FSGM_CalSafeBoundaries[i][1]-pp[i])/vv[i];
      for(j=0;j<3;++j)
        {
          if(j!=i) pp1[i][j] = pp[j]+lambda*vv[j];
          else pp1[i][j] = m_fsgm->FSGM_CalSafeBoundaries[i][1];
        }
      segmentlength[i] = 0;
      for(j=0;j<3;++j) segmentlength[i] += (pp0[i][j]-pp1[i][j])*(pp0[i][j]-pp1[i][j]);
      segmentlength[i] = sqrt(segmentlength[i]);
      if(segmentlength[i]<segmentlengthmin)
        {
          segmentlengthmin = segmentlength[i];
          idir = i;
        }
    }

  double checkdirection = (pp1[idir][0]-pp0[idir][0])*vv[0]+(pp1[idir][1]-pp0[idir][1])*vv[1]+(pp1[idir][2]-pp0[idir][2])*vv[2];
  
  if(checkdirection>0)
    {
      for(i=0;i<3;++i)
        {
          ppstart[i] = pp0[idir][i];
          ppend[i] = pp1[idir][i];
        }
    }
  else
    {
      for(i=0;i<3;++i)
        {
          ppstart[i] = pp1[idir][i];
          ppend[i] = pp0[idir][i];
        }
    }
}

double NewMultiFullShowerDevelopmentDescription::GoThroughCrack(double *pp, double *vv, double zstep, int nstepmax)
{
  double ppp[3];
  double lambda = 0;

  int i,j;
  for(i=0;i<nstepmax;++i)
    {
      lambda += zstep;
      for(j=0;j<3;++j) ppp[j] = pp[j]+lambda*vv[j];
      m_fsgm->WhereInCal(ppp,currentwhereincal);
      if(currentwhereincal[0]!=2)
        {
          lambda -= zstep;
          break;
        }
    }
  return lambda;
}

double NewMultiFullShowerDevelopmentDescription::GetCrackAngle(double *pp, double *vv, double lambdain, double *pptraj)
{
  double pp0[3];
  double pp1[3];
  double ppb[3];
  double vvcrack[3];
  double vvtraj[3];
  int i,j;
  for(j=0;j<3;++j) pp0[j] = pp[j];
  for(j=0;j<3;++j) pp1[j] = pp[j]+lambdain*vv[j];

  for(j=0;j<3;++j) ppb[j] = pp1[j];
  ppb[2] = m_fsgm->FSGM_calZBot;

  double mynorm = 0;

  for(j=0;j<3;++j) vvcrack[j] = ppb[j]-pp0[j];
  mynorm = 0;
  for(j=0;j<3;++j) mynorm += vvcrack[j]*vvcrack[j];
  mynorm = sqrt(mynorm);
  if(mynorm>0)
    for(j=0;j<3;++j) vvcrack[j] /= mynorm;

  for(j=0;j<3;++j) vvtraj[j] = pp1[j]-pptraj[j];
  mynorm = 0;
  for(j=0;j<3;++j) mynorm += vvtraj[j]*vvtraj[j];
  mynorm = sqrt(mynorm);
  if(mynorm>0)
    for(j=0;j<3;++j) vvtraj[j] /= mynorm;

  double myprod = 0;
  for(j=0;j<3;++j) myprod += vvcrack[j]*vvtraj[j];

  double myang = acos(myprod);

  return myang;
}


bool NewMultiFullShowerDevelopmentDescription::Compute(double *pp, double *vv, double startx0_input, double zstep_input, double totrlnmax)
{
  ZStepRef = zstep_input;

  int i,j,k,ii;
  int whereincal[4];
  int whereincal2[4];
  // whereincal[0] = material
  // whereincal[1] = layer
  // whereincal[2] = tower
  // whereincal[3] = column

  NStep = 0;

  for(ii=0;ii<=NDevelopment;++ii)
    {
      x0maxshower[ii] = XMax[ii];
      if(x0maxshower[ii]<2)
        {
          // NewMultiFullShowerDevelopmentDescription x0maxshower is too small
          x0maxshower[ii] = 2;
        }
      startx0[ii] = startx0_input;
      lastx0[ii] = 0;
      crackmaxfrac[ii] = 0;
      x0crackmaxfrac[ii] = 0;
    }

  double ppc[3];
  double ppcc[3];
  double ppcc2[3];
  double pp0[3];
  double pp1[3];
  double vv0[3];
  double vv1[3];
  double vv2[3];
  double vnorm = sqrt(vv[0]*vv[0]+vv[1]*vv[1]+vv[2]*vv[2]);
  if(vnorm==0) return false;

  vv2[0] = vv[0]/vnorm;
  vv2[1] = vv[1]/vnorm;
  vv2[2] = vv[2]/vnorm;

  GetTrajectorySegment(pp,vv2,pp0,pp1);

  vv0[0] = vv2[1];
  vv0[1] = -vv2[0];
  vv0[2] = 0;
  vnorm = sqrt(vv0[0]*vv0[0]+vv0[1]*vv0[1]+vv0[2]*vv0[2]);
  if(vnorm>0)
    {
      vv0[0] = vv0[0]/vnorm;
      vv0[1] = vv0[1]/vnorm;
      vv0[2] = vv0[2]/vnorm;
    }
  else
    {
      vv0[0] = 1;
      vv0[1] = 0;
      vv0[2] = 0;
    }
  vv1[0] = vv2[1]*vv0[2]-vv2[2]*vv0[1];
  vv1[1] = vv2[2]*vv0[0]-vv2[0]*vv0[2];
  vv1[2] = vv2[0]*vv0[1]-vv2[1]*vv0[0];

  if(vv1[2]>0)
    {
      for(i=0;i<3;++i)
        {
          vv0[i] = -vv0[i];
          vv1[i] = -vv1[i];
        }
    }

  double lambda = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );

  int nhist = (int)floor(lambda/ZStepRef);
  if(nhist>=FSDD_NSTEPS_MAX)
    {
      // NewMultiFullShowerDevelopmentDescription::Compute : nhist >= FSDD_NSTEPS_MAX
      return false;
    }
  NStep = nhist;
  ZStep = ZStepRef;

  double ppref[3];
  ppref[2] = 0;
  double lambdaref = -100;
  double zreftmp;
  for(i=0;i<=8;++i)
    {
      zreftmp = m_fsgm->FSGM_calZTop-m_fsgm->FSGM_cellVertPitch*(double)i;
      lambdaref = -100;
      if(vv2[2]!=0)
        lambdaref = (zreftmp-pp0[2])/vv2[2];
      if(lambdaref>0)
        {
          for(j=0;j<3;++j) ppref[j] = pp0[j]+lambdaref*vv2[j];
          break;
        }
    }
  double myzstep;
  int mynstep;
  double lengthref;

  if(ppref[2]<0 && vv2[2]<-0.001)
    {
      lengthref = fabs(m_fsgm->FSGM_cellVertPitch/vv2[2]);
      mynstep = (int)floor(lengthref/ZStepRef)+1;
      ZStep = lengthref/(double)mynstep;
      // change pp0
      lambdaref = sqrt( (ppref[0]-pp0[0])*(ppref[0]-pp0[0])+(ppref[1]-pp0[1])*(ppref[1]-pp0[1])+(ppref[2]-pp0[2])*(ppref[2]-pp0[2]) );
      mynstep = (int)floor(lambdaref/ZStep)+1;
      lambdaref = ZStep*(double)mynstep;
      for(j=0;j<3;++j) pp0[j] = ppref[j]-lambdaref*vv2[j];
      // change pp1
      lambdaref = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );
      mynstep = (int)floor(lambdaref/ZStep)+1;
      lambdaref = ZStep*(double)mynstep;
      for(j=0;j<3;++j) pp1[j] = pp0[j]+lambdaref*vv2[j];

      lambda = sqrt( (pp1[0]-pp0[0])*(pp1[0]-pp0[0])+(pp1[1]-pp0[1])*(pp1[1]-pp0[1])+(pp1[2]-pp0[2])*(pp1[2]-pp0[2]) );      
      nhist = (int)floor(lambda/ZStep);
      if(nhist>=FSDD_NSTEPS_MAX)
        {
          // NewMultiFullShowerDevelopmentDescription::Compute : nhist >= FSDD_NSTEPS_MAX
          return false;
        }
      NStep = nhist;
    }

  int nstepincrackmax = (int)floor(FSDD_XCRK/FSDD_XCSI);

  double z2X0mat;

  for(ii=0;ii<=NDevelopment;++ii)
    {
      X0[ii][0] = startx0_input;
      lastx0[ii] = startx0_input;
    }

  double radius, relradius;
  double effradius, releffradius;
  double radialprofile;

  double etotdep[FSDD_NMAX];
  double x0position[FSDD_NMAX];
  for(ii=0;ii<=NDevelopment;++ii)
    x0position[ii] = X0[ii][0];
  double relx0position[FSDD_NMAX];
  double meancracklength[FSDD_NMAX];
  double meancracklengthw[FSDD_NMAX];

  for(ii=0;ii<=NDevelopment;++ii)
    {
      totx0cal[ii] = 0;
      totx0crack[ii] = 0;
    }
  for(ii=0;ii<=NDevelopment;++ii)
    for(j=0;j<8;++j)
      totx0lay[ii][j] = 0;
  for(ii=0;ii<=NDevelopment;++ii)
    for(j=0;j<8;++j)
      posx0lay[ii][j] = 0;

  int isatxtal;

  bool InCal = false;

  for(i=0;i<m_fsgm->FSGM_NCircle;++i)
    {
      m_fsgm->FSGM_LambdaCircle[i] = -0.5*ZStep;
      m_fsgm->FSGM_TCircle[i] = 1000.0;
      m_fsgm->FSGM_TCircle2[i] = -1000.0;
      m_fsgm->FSGM_ZCrackCircle[i] = 0.0;
      for(j=0;j<FSDD_NMAX;++j)
        m_fsgm->FSGM_ScaleCircle[i][j] = 1.0;
    }
  double lambdacrack,lambdacrack2,myscale,myskin,myevap,myevapmax,myevap2,myrtailright;
  double wfcrack,myang,angshower;
  double pptraj[3];

  double explosionangle,explosionradius,explosionradius2;
  double mynorm;
  double ppbot[3];
  double vvcrack[3];
  double vvperp[3];
  double vvshower[3];

  double gapzmean = 0;
  double gapzmeanw = 0;
  for(j=0;j<NStep;++j)
    {
      lambda = ZStep*((double)j+0.5);
      ppc[0] = pp0[0]+vv2[0]*lambda;
      ppc[1] = pp0[1]+vv2[1]*lambda;
      ppc[2] = pp0[2]+vv2[2]*lambda;
      m_fsgm->WhereInCal(ppc,whereincal);
      if(whereincal[0]==2)
        {
          gapzmean += ppc[2];
          gapzmeanw += 1.;
        }
    }
  if(gapzmeanw>0) gapzmean /= gapzmeanw;

  myevapmax= 0;

  double evapmaxpar0 = -0.242958+0.640567*fabs(vv[2]);
  if(evapmaxpar0<0.21) evapmaxpar0 = 0.21;
  double evapmaxpar1 = 130.;

  double tailrightpar0 = 0.5;
  double tailrightpar1 = 0.6;
  double tailrightpar2 = 22.0;

  double tailrightratiomax = tailrightpar0;
  if(fabs(vv[2])>tailrightpar1) tailrightratiomax = tailrightpar0+tailrightpar2*(fabs(vv[2])-tailrightpar1)*(fabs(vv[2])-tailrightpar1);

  for(j=0;j<NStep;++j)
    {
      if(j>0)
        {
          for(ii=0;ii<=NDevelopment;++ii)
            X0[ii][j] = dX0[ii][j-1] + X0[ii][j-1];
        }

      lambda = ZStep*((double)j+0.5);
      ppc[0] = pp0[0]+vv2[0]*lambda;
      ppc[1] = pp0[1]+vv2[1]*lambda;
      ppc[2] = pp0[2]+vv2[2]*lambda;
      // Find out if the trajectory gets out of the cal
      m_fsgm->WhereInCal(ppc,whereincal);
      if(!InCal && whereincal[0]==1 || whereincal[0]==2) InCal = true;
      //
      for(ii=0;ii<=NDevelopment;++ii)
        for(i=0;i<4;++i)
          materialfraction[ii][i][j] = 0;
      for(ii=0;ii<=NDevelopment;++ii)
        for(i=0;i<8;++i)
          layerfraction[ii][i][j] = 0;
      for(ii=0;ii<=NDevelopment;++ii)
        for(i=0;i<NXtal;++i)
          xtalfraction[ii][i][j] = 0;
      //
      for(ii=0;ii<=NDevelopment;++ii)
        {
          etotdep[ii] = 0;
          relx0position[ii] = x0position[ii]/XMax[ii];
          releffradius = wideningfactor*m_fsgm->GetEffectiveRadius(relx0position[ii],RadialContainedFraction);
          RM[ii][j] = releffradius;
          meancracklength[ii] = 0;
          meancracklengthw[ii] = 0;
        }
      effradius = 3.0*FSDD_MOLRAD;
      //
      for(i=0;i<m_fsgm->FSGM_NCircle;++i)
        {
          lambdacrack = 0.0;

          m_fsgm->FSGM_LambdaCircle[i] += ZStep;

          ppcc[0] = pp0[0]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[0]+m_fsgm->FSGM_YCircle[i]*vv1[0])+m_fsgm->FSGM_LambdaCircle[i]*vv2[0];
          ppcc[1] = pp0[1]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[1]+m_fsgm->FSGM_YCircle[i]*vv1[1])+m_fsgm->FSGM_LambdaCircle[i]*vv2[1];
          ppcc[2] = pp0[2]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[2]+m_fsgm->FSGM_YCircle[i]*vv1[2])+m_fsgm->FSGM_LambdaCircle[i]*vv2[2];

          m_fsgm->WhereInCal(ppcc,whereincal);

          // evap before crack : not used yet
//           if(whereincal[0]==1 && m_fsgm->FSGM_TCircle2[i]<0)
//             {
//               ppcc2[0] = ppcc[0]+FSDD_XCSI*vv2[0];
//               ppcc2[1] = ppcc[1]+FSDD_XCSI*vv2[1];
//               ppcc2[2] = ppcc[2]+FSDD_XCSI*vv2[2];

//               m_fsgm->WhereInCal(ppcc2,whereincal2);

//            if(whereincal2[0]==2)
//              {
//                m_fsgm->FSGM_TCircle2[i] = 1;
//                //
//                lambdacrack2 = GoThroughCrack(ppcc2,vv2,ZStep,1000);
//                ppcc2[0] = ppcc[0]+(FSDD_XCSI+lambdacrack2)*vv2[0];
//                ppcc2[1] = ppcc[1]+(FSDD_XCSI+lambdacrack2)*vv2[1];
//                ppcc2[2] = ppcc[2]+(FSDD_XCSI+lambdacrack2)*vv2[2];
//                m_fsgm->FSGM_ZCrackCircle[i] = m_fsgm->GetZCrack(ppcc2,pp0,vv2);
//              }
//             }
//           else
//             m_fsgm->FSGM_TCircle2[i] -= ZStep/FSDD_XCSI;

          if(whereincal[0]==2)
            {
              lambdacrack = GoThroughCrack(ppcc,vv2,ZStep,nstepincrackmax);
              m_fsgm->FSGM_TCircle[i] = -ZStep/FSDD_XCSI;
              m_fsgm->FSGM_TCircle2[i] = -1000;
              //
              m_fsgm->FSGM_LambdaCircle[i] += lambdacrack;
              ppcc[0] = pp0[0]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[0]+m_fsgm->FSGM_YCircle[i]*vv1[0])+m_fsgm->FSGM_LambdaCircle[i]*vv2[0];
              ppcc[1] = pp0[1]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[1]+m_fsgm->FSGM_YCircle[i]*vv1[1])+m_fsgm->FSGM_LambdaCircle[i]*vv2[1];
              ppcc[2] = pp0[2]+effradius*(m_fsgm->FSGM_XCircle[i]*vv0[2]+m_fsgm->FSGM_YCircle[i]*vv1[2])+m_fsgm->FSGM_LambdaCircle[i]*vv2[2];
              m_fsgm->WhereInCal(ppcc,whereincal);

              m_fsgm->FSGM_ZCrackCircle[i] = m_fsgm->GetZCrack(ppcc,pp0,vv2);         
            }
          else
            m_fsgm->FSGM_TCircle[i] += ZStep/FSDD_XCSI;

          myscale = 1.0;
          myrtailright = 1.0;
          if(m_fsgm->FSGM_TCircle[i]<=4)
            {
              myscale = 0.5*exp(-m_fsgm->FSGM_TCircle[i]);
              if(m_fsgm->FSGM_YCircle[i]<0)
                myrtailright = 1+tailrightratiomax*exp(-m_fsgm->FSGM_TCircle[i]/2.);
            }
          myevap = 1;
          if(m_fsgm->FSGM_TCircle[i]<=4 && ppcc[2]>m_fsgm->FSGM_calZBot)
            {
              myevapmax= evapmaxpar0*(1-(m_fsgm->FSGM_ZCrackCircle[i]-m_fsgm->FSGM_calZBot)/evapmaxpar1);
              if(myevapmax<0) myevapmax= 0;
              myevap = 1-myevapmax*exp(-m_fsgm->FSGM_TCircle[i]);
              if(myevap<0) myevap = 0;
              if(myevap>1) myevap = 1;
            }
          myevap2 = 1;

          // evap before crack : not used yet
//           if(m_fsgm->FSGM_TCircle2[i]>=0 && ppcc[2]>m_fsgm->FSGM_calZBot)
//             {
//            myevapmax= evapmaxpar0*(1-(m_fsgm->FSGM_ZCrackCircle[i]-m_fsgm->FSGM_calZBot)/evapmaxpar1);
//               myevap2 = 1-myevapmax*exp(-m_fsgm->FSGM_TCircle2[i]*2.);
//               if(myevap2<0) myevap2 = 0;
//               if(myevap2>1) myevap2 = 1;
//          }
//        myevap *= myevap2;

          
          radius = m_fsgm->FSGM_RCircle[i]*effradius;
          relradius = radius/FSDD_MOLRAD;
          //
          for(ii=0;ii<=NDevelopment;++ii)
            {
              if(m_fsgm->FSGM_TCircle[i]>4)
                radialprofile = m_fsgm->GetRadialProfile(relx0position[ii],relradius/wideningfactor)*m_fsgm->FSGM_WCircle[i];
              else
                radialprofile = m_fsgm->GetRadialProfileMod2(relx0position[ii],relradius/wideningfactor,myscale,myrtailright)*m_fsgm->FSGM_WCircle[i];
              etotdep[ii] += radialprofile;
              //
              materialfraction[ii][whereincal[0]][j] += radialprofile;
              if(whereincal[0]==1)
                {
                  layerfraction[ii][whereincal[1]][j] += radialprofile*myevap;
                  //
                  isatxtal = OffSatu[whereincal[2]][whereincal[1]][whereincal[3]];
                  if(isatxtal>=0) xtalfraction[ii][isatxtal][j] += radialprofile*myevap;
                }
              else if(whereincal[0]==2)
                meancracklength[ii] += (lambdacrack/ZStep+1.0)*radialprofile; // in ZStep units
            }
        }
      //
      for(ii=0;ii<=NDevelopment;++ii)
        {
          if(etotdep[ii]>0)
            {
              for(i=0;i<4;++i)
                materialfraction[ii][i][j] /= etotdep[ii];
              for(i=0;i<8;++i)
                layerfraction[ii][i][j] /= etotdep[ii];
              for(i=0;i<NXtal;++i)
                xtalfraction[ii][i][j] /= etotdep[ii];
              //
              meancracklength[ii] /= etotdep[ii];
            }
          
          // Ignore the FakeCsI volume contribution when computing lastx0
          z2X0mat = materialfraction[ii][1][j]/FSDD_XCSI + meancracklength[ii]/FSDD_XCRK;

          lastx0[ii] += ZStep*z2X0mat;
          // But take into account for the shower development
          if(InCal) z2X0mat += materialfraction[ii][3][j]/FSDD_XCSI;
          //
          dX0[ii][j] = ZStep*z2X0mat;
          x0position[ii] = X0[ii][j];
          totx0cal[ii] += materialfraction[ii][1][j]*dX0[ii][j];
          totx0crack[ii] += materialfraction[ii][2][j]*dX0[ii][j];
          for(i=0;i<8;++i) totx0lay[ii][i] += layerfraction[ii][i][j]*dX0[ii][j];
          
          if(materialfraction[ii][2][j]>crackmaxfrac[ii])
            {
              x0crackmaxfrac[ii] = X0[ii][j];
              crackmaxfrac[ii] = materialfraction[ii][2][j];
            }
        }
      if(X0[NDevelopment][j]>totrlnmax) 
        {
          NStep = j+1;
          break;
        }
    }

  for(ii=0;ii<=NDevelopment;++ii)
    {
      dX0[ii][NStep] = 0;
      X0[ii][NStep] = dX0[ii][NStep-1] + X0[ii][NStep-1];
      materialfraction[ii][0][NStep] = 1;
      for(i=1;i<4;++i) materialfraction[ii][i][NStep] = 0;
      for(i=0;i<8;++i) layerfraction[ii][i][NStep] = 0;
      for(i=0;i<NXtal;++i) xtalfraction[ii][i][NStep] = 0;
    }

  return true;
}

void NewMultiFullShowerDevelopmentDescription::SetWideningFactor(double widfact)
{
  wideningfactor = widfact;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

NewFullShowerDevelopmentDescriptionManager::NewFullShowerDevelopmentDescriptionManager(IGlastDetSvc *m_detSvc_input, int nxmax_input, double xmax0, double dxmax, double zstep_input, double radialcontainedfraction_input, double x0step)
  :m_detSvc(m_detSvc_input)
{

  // Initialize geometry
  m_fsgm = new NewFullShowerGeometryManager(m_detSvc);

  NDevelopment = 0;
  mintotx0cal = 0;
  maxtotx0cal = 0;
  int i;

  for(i=0;i<8;++i)
    {
      meantotx0lay[i] = 0;;
      meanposx0lay[i] = -1;
    }
  

  for(i=0;i<FSDD_NMAX;++i)
    {
      FSDDMM[i] = NULL;
      FSDDX0[i] = NULL;
    }
  CurrentFSDD = NULL;

  int nxmax = nxmax_input;
  if(nxmax>FSDD_NMAX-1)
    {
      // NewFullShowerDevelopmentDescriptionManager nxmax>FSDD_NMAX-1
      nxmax = FSDD_NMAX-1;
    }

  NDevelopment = nxmax;
  DXMax = dxmax;
  ZStep = zstep_input;
  ZStepRef = zstep_input;
  RadialContainedFraction = radialcontainedfraction_input;
  X0Step = x0step;

  MFSDDMM = new NewMultiFullShowerDevelopmentDescription(m_fsgm,NDevelopment,xmax0,DXMax,ZStepRef,RadialContainedFraction);
  
  for(i=0;i<=NDevelopment;++i)
    {
      XMax[i] = xmax0+DXMax*(double)i;
      FSDDMM[i] = new NewFullShowerDevelopmentDescription(m_fsgm,0,ZStepRef,RadialContainedFraction);
      FSDDX0[i] = new NewFullShowerDevelopmentDescription(m_fsgm,1,X0Step,RadialContainedFraction);
    }
  CurrentFSDD = new NewFullShowerDevelopmentDescription(m_fsgm,1,X0Step,RadialContainedFraction);
}

NewFullShowerDevelopmentDescriptionManager::~NewFullShowerDevelopmentDescriptionManager()
{
  int i;
  if(MFSDDMM) delete MFSDDMM;
  for(i=0;i<FSDD_NMAX;++i)
    {
      if(FSDDMM[i]!=NULL) delete FSDDMM[i];
      if(FSDDX0[i]!=NULL) delete FSDDX0[i];
    }
  if(CurrentFSDD!=NULL) delete CurrentFSDD;
}

bool NewFullShowerDevelopmentDescriptionManager::Compute(double *pp, double *vv, double startx0_input, double zstep_input, double totrlnmax)
{
  ZStepRef = zstep_input;
  
  int i,j,k,l;
  mintotx0cal = 99999999;
  maxtotx0cal = -99999999;

  MFSDDMM->NXtal = NXtal;
  for(i=0;i<16;++i)
    for(j=0;j<8;++j)
      for(k=0;k<12;++k)
        MFSDDMM->OffSatu[i][j][k] = OffSatu[i][j][k];

  for(l=0;l<=NDevelopment;++l)
    {
      FSDDMM[l]->NXtal = NXtal;
      FSDDX0[l]->NXtal = NXtal;
      //
      for(i=0;i<16;++i)
        for(j=0;j<8;++j)
          for(k=0;k<12;++k)
            {
              FSDDMM[l]->OffSatu[i][j][k] = OffSatu[i][j][k];
              FSDDX0[l]->OffSatu[i][j][k] = OffSatu[i][j][k];
            }
    }

  int optmulti = 1;

  if(optmulti) MFSDDMM->Compute(pp,vv,startx0_input,ZStepRef,totrlnmax);
  
  for(i=0;i<=NDevelopment;++i)
    {
      if(!optmulti)
        {
          if(!FSDDMM[i]->Compute(pp,vv,startx0_input,XMax[i],ZStepRef,totrlnmax)) return false;
        }
      else
        {
          FSDDMM[i]->FillFromMultiFullShowerDevelopmentDescription(MFSDDMM,i);
          FSDDMM[i]->RemoveEmptySteps();
        }
      //
      if(!FSDDX0[i]->ConvertToFixedX0(X0Step,FSDDMM[i])) return false;
      if(FSDDX0[i]->totx0cal<mintotx0cal) mintotx0cal = FSDDX0[i]->totx0cal;
      if(FSDDX0[i]->totx0cal>maxtotx0cal) maxtotx0cal = FSDDX0[i]->totx0cal;
    }
  
  double meanw = 0;
  for(j=0;j<8;++j)
    {
      meantotx0lay[j] = 0;;
      meanposx0lay[j] = 0;
      meanw = 0;
      for(i=0;i<=NDevelopment;++i)
        {
          meantotx0lay[j] += FSDDX0[i]->totx0lay[j];
          if(FSDDX0[i]->posx0lay[j]>-1)
            {
              meanw += 1.;
              meanposx0lay[j] += FSDDX0[i]->posx0lay[j];
            }
        }
      meantotx0lay[j] /= (double)(NDevelopment+1);
      if(meanw>0)
        meanposx0lay[j] /= meanw;
      else
        meanposx0lay[j] = -1;
    }
  
  return true;
}

void NewFullShowerDevelopmentDescriptionManager::FillCurrentFSDD(double showerxmax)
{  
  int ish;
  double interpol;

  if(showerxmax<=XMax[0])
    {
      ish = 0;
      interpol = 0;
    }
  else if(showerxmax>=XMax[NDevelopment])
    {
      ish = NDevelopment-1;
      interpol = 1;
    }
  else
    {
      ish = (int)floor((showerxmax-XMax[0])/DXMax);
      interpol = (showerxmax-XMax[ish])/DXMax;
    }

  int i,j;
  CurrentFSDD->NStep = FSDDX0[ish]->NStep;
  if(FSDDX0[ish+1]->NStep>CurrentFSDD->NStep)
    CurrentFSDD->NStep = FSDDX0[ish+1]->NStep;
  CurrentFSDD->x0maxshower = showerxmax;
  CurrentFSDD->crackmaxfrac = (1-interpol)*FSDDX0[ish]->crackmaxfrac+interpol*FSDDX0[ish+1]->crackmaxfrac;
  CurrentFSDD->x0crackmaxfrac = (1-interpol)*FSDDX0[ish]->x0crackmaxfrac+interpol*FSDDX0[ish+1]->x0crackmaxfrac;
  CurrentFSDD->startx0 = (1-interpol)*FSDDX0[ish]->startx0+interpol*FSDDX0[ish+1]->startx0;
  CurrentFSDD->lastx0 = (1-interpol)*FSDDX0[ish]->lastx0+interpol*FSDDX0[ish+1]->lastx0;
  CurrentFSDD->totx0cal = (1-interpol)*FSDDX0[ish]->totx0cal+interpol*FSDDX0[ish+1]->totx0cal;
  CurrentFSDD->totx0crack = (1-interpol)*FSDDX0[ish]->totx0crack+interpol*FSDDX0[ish+1]->totx0crack;
  for(i=0;i<8;++i) CurrentFSDD->totx0lay[i] = (1-interpol)*FSDDX0[ish]->totx0lay[i]+interpol*FSDDX0[ish+1]->totx0lay[i];
  for(i=0;i<8;++i) CurrentFSDD->posx0lay[i] = (1-interpol)*FSDDX0[ish]->posx0lay[i]+interpol*FSDDX0[ish+1]->posx0lay[i];
  
  for(i=0;i<CurrentFSDD->NStep;++i)
    {
      CurrentFSDD->dX0[i] = FSDDX0[ish]->dX0[i];
      CurrentFSDD->X0[i] = FSDDX0[ish]->X0[i];
      CurrentFSDD->RM[i] = FSDDX0[ish]->RM[i];
      for(j=0;j<4;++j)
        CurrentFSDD->materialfraction[j][i] = (1-interpol)*FSDDX0[ish]->materialfraction[j][i]
          +interpol*FSDDX0[ish+1]->materialfraction[j][i];
      for(j=0;j<8;++j)
        CurrentFSDD->layerfraction[j][i] = (1-interpol)*FSDDX0[ish]->layerfraction[j][i]
          +interpol*FSDDX0[ish+1]->layerfraction[j][i];
      for(j=0;j<NXtal;++j)
        CurrentFSDD->xtalfraction[j][i] = (1-interpol)*FSDDX0[ish]->xtalfraction[j][i]
          +interpol*FSDDX0[ish+1]->xtalfraction[j][i];
    }

  for(i=CurrentFSDD->NStep;i<FSDD_NSTEPS_MAX;++i)
    {
      CurrentFSDD->dX0[i] = FSDDX0[ish]->dX0[i];
      CurrentFSDD->X0[i] = FSDDX0[ish]->X0[i];
      CurrentFSDD->RM[i] = FSDDX0[ish]->RM[i];
      CurrentFSDD->materialfraction[0][i] = 1;
      for(j=1;j<3;++j) CurrentFSDD->materialfraction[j][i] = 0;
      for(j=0;j<8;++j) CurrentFSDD->layerfraction[j][i] = 0;
      break;
    }

  return;
}

void NewFullShowerDevelopmentDescriptionManager::SetWideningFactor(double widfact)
{  
  MFSDDMM->SetWideningFactor(widfact);
  int i;
  for(i=0;i<=NDevelopment;++i)
    {
      FSDDMM[i]->SetWideningFactor(widfact);
      FSDDX0[i]->SetWideningFactor(widfact);
    }
  CurrentFSDD->SetWideningFactor(widfact);
}

void NewFullShowerDevelopmentDescriptionManager::SetCrackExtinctionFactor(double crackext)
{  
  MFSDDMM->crackextinction = crackext;
  int i;
  for(i=0;i<=NDevelopment;++i)
    {
      FSDDMM[i]->crackextinction = crackext;
      FSDDX0[i]->crackextinction = crackext;
    }
  CurrentFSDD->crackextinction = crackext;
}
