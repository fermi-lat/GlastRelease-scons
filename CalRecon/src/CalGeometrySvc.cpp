
#include "Gaudi/MessageSvc/MsgStream.h"
#include "Gaudi/Kernel/SvcFactory.h"

#include "CalRecon/CalGeometrySvc.h"

static const SvcFactory<CalGeometrySvc> s_factory;
const ISvcFactory& CalGeometrySvcFactory = s_factory;


CalGeometrySvc::CalGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator)

:Service(name, pSvcLocator) 
{
	declareProperty("geometryType", m_geoType = 0);
	
	
		
// testbeam geometry (default)


		m_nviews  = 2;
		m_nlayers = 4;
		m_nLogs   = 10;
		m_nmodx   = 1;
		m_nmody   = 1;

		m_Z0          = -25.9602; 
		m_layerWidth  =   31.05;
		m_layerHeight =    2.614;

		m_logWidth    =  3.05;
		m_logLength   = 31.05;
		m_logHeight   =  2.35;
		m_logGap      =  0.06;



	



}


//##############################################
StatusCode CalGeometrySvc::initialize()
//##############################################
{
    StatusCode sc = StatusCode::SUCCESS;
	Service::initialize();
	setProperties();

	if (m_geoType == 1)
	{
	
		// Flight geometry
		m_nviews  = 2;
		m_nlayers = 4;
		m_nLogs   = 12;
		m_nmodx   = 4;
		m_nmody   = 4;

		m_modWidth    = 39.37;
		m_Z0          = -34.094;   

		m_layerWidth  =   36.99;
		m_layerHeight =    2.446;

		m_logWidth    =  3.00;
		m_logLength   = 36.99;
		m_logHeight   =  2.10;
		m_logGap      =  0.084;


	}
	return sc;
}


//##############################################
StatusCode CalGeometrySvc::finalize()
//##############################################
{
    return StatusCode::SUCCESS;
}



//----------------------------------------------------
//##############################################
detGeo CalGeometrySvc::getLayer(int ilayer, detGeo::axis a)
//##############################################
{
	double xpos = 0.;
	double ypos = 0.;
	double zfar = 0.;
	if (a == detGeo::Y) zfar = layerHeight();
	double zpos = Z0()+2.*ilayer*layerHeight()+zfar;

	double xsize = layerWidth();
	double ysize = layerWidth();
	double zsize = layerHeight();

	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*zsize);

	detGeo layer(ilayer,a,ilayer,P,S);
	return layer;
}
//##############################################
detGeo CalGeometrySvc::getLog(int ilayer, detGeo::axis a, int ilog, idents::ModuleId mod)
//##############################################
{

	double xmod = (mod.ix()-(m_nmodx+1)*0.5)*m_modWidth;
	double ymod = (mod.iy()-(m_nmody+1)*0.5)*m_modWidth;
	Vector modcenter(xmod,ymod,0);
	detGeo log = getLog(ilayer, a, ilog);
	log.setPosition(log.position()+modcenter);
	return log;
}
//##############################################
detGeo CalGeometrySvc::getLog(int ilayer, detGeo::axis a, int ilog)
//##############################################
{
	detGeo layer = getLayer(ilayer,a);
	double xpos = layer.position().x();
	double ypos = layer.position().y();
	double zpos = layer.position().z();
	
	int nlogs = numLogs();
	double posRef = -0.5*(nlogs-1)*(logWidth()+logGap());
	double pos = posRef+ilog*(logWidth()+logGap());
	if (fabs(pos) < 1e-5) pos =0.; 
	if (a == detGeo::X) xpos = pos;
	else ypos = pos;

	double xsize = logLength();
	double ysize = logLength();
	if (a == detGeo::X) xsize = logWidth();
	else ysize = logWidth();
	double zsize = logHeight();
	
	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*zsize);

	detGeo log(ilayer,a,ilog,P,S);

	return log;
}
