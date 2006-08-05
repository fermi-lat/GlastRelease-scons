#ifndef ANCILLARYGEOMETRY_HH
#define ANCILLARYGEOMETRY_HH
#include <string>
#include "AdfEvent/TaggerParameters.h"


/*
Simple Geometry configuration file for ancyllary system:
0.48
ModId    X      Y       Z       V1      D1      V2      D2
0       0.0     0.0     1500.   Y       +       X       +
1       5.618   -1.46   1000.   Y       +       X       +
2       11.55   0.0     -1000.  X       +       Y       -
3       34.55   0.0     -1300.  X       +       Y       -
 */
// NOTICE THAT:
// Y is: strip along Z (1); Z is strip along Y (0)



namespace AncillaryData
{
  class AncillaryGeometry
    {
    public:
      AncillaryGeometry(std::string GeometryfileName);
      double getX(unsigned int Module) {return X[Module];}
      double getY(unsigned int Module) {return Y[Module];}
      double getZ(unsigned int Module) {return Z[Module];}
      unsigned int getView1(unsigned int Module){return V1[Module];}
      unsigned int getView2(unsigned int Module){return V2[Module];}
      unsigned int getDirection1(unsigned int Module){return D1[Module];}
      unsigned int getDirection2(unsigned int Module){return D2[Module];}
      double getBL(){return BL;}
      void print();
    private:
      double BL;
      char title[100];
      double X[N_MODULES],Y[N_MODULES],Z[N_MODULES];
      unsigned int V1[N_MODULES], V2[N_MODULES], D1[N_MODULES], D2[N_MODULES]; 
    };

}

#endif
