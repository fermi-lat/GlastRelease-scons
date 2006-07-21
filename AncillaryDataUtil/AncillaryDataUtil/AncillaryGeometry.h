#ifndef ANCILLARYGEOMETRY_HH
#define ANCILLARYGEOMETRY_HH
#include <string>
#include "AdfEvent/TaggerParameters.h"


/*
ModId 	X	Y	Z	Tx	Ty	Tz	View(Y||Z)
1	-1500	0	0	0	0	0	Y
2	-1000	0	0	0	0	0	Y
3	1000	230	0	0	0	0	Y
4	1300	240	0	0	0	0	Y
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
      double getTx(unsigned int Module){return Tx[Module];}
      double getTy(unsigned int Module){return Ty[Module];}
      double getTz(unsigned int Module){return Tz[Module];}
      double getBL(){return BL;}
      unsigned int getView(unsigned int Module){return View[Module];}
      void print();
    private:
      double BL;
      char title[100];
      double X[N_MODULES],Y[N_MODULES],Z[N_MODULES],Tx[N_MODULES],Ty[N_MODULES],Tz[N_MODULES];
      unsigned int View[N_MODULES]; 
    };

}

#endif
