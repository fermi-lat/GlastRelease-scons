
#include "CalRecon/calorimeterGeo.h"

//--------------- parameters -----------------------
int calorimeterGeo::m_nviews  = 2;
int calorimeterGeo::m_nlayers = 4;
// int calorimeterGeo::m_nLogs   = 10;
int calorimeterGeo::m_nLogs   = 12;
int calorimeterGeo::m_nmodx   = 4;
int calorimeterGeo::m_nmody   = 4;

/*
// for tb geometry
double calorimeterGeo::m_Z0          = -25.9602; 
double calorimeterGeo::m_layerWidth  =   31.05;
double calorimeterGeo::m_layerHeight =    2.614;

double calorimeterGeo::m_logWidth    =  3.05;
double calorimeterGeo::m_logLength   = 31.05;
double calorimeterGeo::m_logHeight   =  2.35;
double calorimeterGeo::m_logGap      =  0.06;
*/

// for flight geometry

double calorimeterGeo::m_modWidth    = 39.37;
double calorimeterGeo::m_Z0          = -34.094;   //  corrected

// double calorimeterGeo::m_Z0          = -25.9602;   // not yet corrected
double calorimeterGeo::m_layerWidth  =   36.99;
double calorimeterGeo::m_layerHeight =    2.446;

double calorimeterGeo::m_logWidth    =  3.00;
double calorimeterGeo::m_logLength   = 36.99;
double calorimeterGeo::m_logHeight   =  2.10;
double calorimeterGeo::m_logGap      =  0.084;


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
detGeo calorimeterGeo::getLog(int ilayer, detGeo::axis a, int ilog, idents::ModuleId mod)
//##############################################
{

	double xmod = (mod.ix()-(m_nmodx-1)*0.5)*m_modWidth;
	double ymod = (mod.iy()-(m_nmody-1)*0.5)*m_modWidth;
	Vector modcenter(xmod,ymod,0);
	detGeo log = getLog(ilayer, a, ilog);
	log.setPosition(log.position()+modcenter);
	return log;
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
