
#include "CalRecon/calorimeterGeo.h"

//--------------- parameters -----------------------
int calorimeterGeo::m_nviews  = 2;
int calorimeterGeo::m_nlayers = 4;
int calorimeterGeo::m_nLogs   = 10;

double calorimeterGeo::m_Z0          = -25.9602; 
double calorimeterGeo::m_layerWidth  =   31.05;
double calorimeterGeo::m_layerHeight =    2.614;

double calorimeterGeo::m_logWidth    =  3.05;
double calorimeterGeo::m_logLength   = 31.05;
double calorimeterGeo::m_logHeight   =  2.35;
double calorimeterGeo::m_logGap      =  0.06;

//----------------------------------------------------
//##############################################
detGeo calorimeterGeo::getLayer(int ilayer, detGeo::axis a)
//##############################################
{
	double xpos = 0.;
	double ypos = 0.;
	double zfar = 0.;
	if (a == detGeo::Y) zfar = layerHeight();
	double zpos = Z0()+2.*ilayer*layerHeight()+zfar;

	double xsize = calorimeterGeo::layerWidth();
	double ysize = calorimeterGeo::layerWidth();
	double zsize = calorimeterGeo::layerHeight();

	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*zsize);

	detGeo layer(ilayer,a,ilayer,P,S);
	return layer;
}

//##############################################
detGeo calorimeterGeo::getLog(int ilayer, detGeo::axis a, int ilog)
//##############################################
{
	detGeo layer = getLayer(ilayer,a);
	double xpos = layer.position().x();
	double ypos = layer.position().y();
	double zpos = layer.position().z();
	
	int nlogs = calorimeterGeo::numLogs();
	double posRef = -0.5*(nlogs-1)*(calorimeterGeo::logWidth()+
		calorimeterGeo::logGap());
	double pos = posRef+ilog*(calorimeterGeo::logWidth()+
		calorimeterGeo::logGap());
	if (fabs(pos) < 1e-5) pos =0.; 
	if (a == detGeo::X) xpos = pos;
	else ypos = pos;

	double xsize = calorimeterGeo::logLength();
	double ysize = calorimeterGeo::logLength();
	if (a == detGeo::X) xsize = calorimeterGeo::logWidth();
	else ysize = calorimeterGeo::logWidth();
	double zsize = calorimeterGeo::logHeight();
	
	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*zsize);

	detGeo log(ilayer,a,ilog,P,S);

	return log;
}
