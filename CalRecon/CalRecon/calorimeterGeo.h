
#ifndef CALORIMETERGEO_H
#define CALORIMETERGEO_H 1

#include "TkrRecon/detGeo.h"


//----------------------------------------------
//
//   calorimeterGeo
//
//     It stores the calorimeter geometry.
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------
//##########################################################
class calorimeterGeo 
//##########################################################
{
public:

	static int numViews()     {return m_nviews;}
	static int numLayers()    {return m_nlayers;}
	static int numLogs()      {return m_nLogs;}

	static double Z0()          {return m_Z0;}
	static double layerWidth()  {return m_layerWidth;}
	static double layerHeight() {return m_layerHeight;}

	static double logWidth()       {return m_logWidth;}
	static double logHeight()      {return m_logHeight;}
	static double logLength()      {return m_logLength;}
	static double logGap()         {return m_logGap;}

	static detGeo getLayer(int ilayer, detGeo::axis a);
	static detGeo getLog(int ilayer, detGeo::axis a, int ilog); 




private:

	static int m_nviews;
	static int m_nlayers;
	static int m_nLogs;

	static double m_Z0;

	static double m_layerWidth;
	static double m_layerHeight;

	static double m_logWidth;
	static double m_logLength;
	static double m_logHeight;
	static double m_logGap;

};

#endif
