
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/SvcFactory.h"

#include "CalRecon/CalGeometrySvc.h"

static const SvcFactory<CalGeometrySvc> s_factory;
const ISvcFactory& CalGeometrySvcFactory = s_factory;


CalGeometrySvc::CalGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator)

:Service(name, pSvcLocator) 
{
	declareProperty("geometryType", m_geoType = 0);
	declareProperty("balloonFlight", m_balloonFlight = 0);
	
	
		
// testbeam geometry (default)


		m_nviews  = 2;
		m_nlayers = 4;
		m_nLogs   = 10;
		m_nmodx   = 1;
		m_nmody   = 1;

//   testbeam geometry 
		m_Z0          = -25.9602; 


        m_layerWidth  =   31.05;
		m_layerHeight =    2.614;

		m_logWidth    =  3.05;
		m_logLength   = 31.05;
		m_logHeight   =  2.35;
		m_logGap      =  0.06;
        m_latt        =  0.35;


	



}


//##############################################
StatusCode CalGeometrySvc::initialize()
//##############################################
{
    StatusCode sc = StatusCode::SUCCESS;
	Service::initialize();
	setProperties();
    
    //    balloon flight geometry
    if (m_geoType == 0 && m_balloonFlight == 1)
//            m_Z0 = -24.346; 
//   calorimeter - tracker distance corrected according to G.Godfrey information
            m_Z0 = -27.416; 

	if (m_geoType == 1)
	{
	
		// Flight geometry, compressed cell design 
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
        m_latt        =  0.35;

	}
	else if (m_geoType == 2)
	{
	
		// Flight geometry, carbon cell design 
		m_nviews  = 2;
		m_nlayers = 4;
		m_nLogs   = 12;
		m_nmodx   = 4;
		m_nmody   = 4;

//		m_modWidth    = 37.35;

//      after grid design modification m_modWidth became                
                m_modWidth    = 37.45;

//		m_Z0          = -32.1425;        
//     after tracker geometry modification m_Z0 became
//		m_Z0          = -32.7756;   

//      after grid design modification m_Z0 became                
                
                m_Z0          = -30.9176;   

		m_layerWidth  =   33.3;
		m_layerHeight =    2.139;

		m_logWidth    =  2.6714;
		m_logLength   = 33.3;
		m_logHeight   =  1.99;
		m_logGap      =  0.113;
        m_latt        =  0.35;

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
CalDetGeo CalGeometrySvc::getLayer(int ilayer, CalDetGeo::axis a)
//##############################################
{
	double xpos = 0.;
	double ypos = 0.;
	double zfar = 0.;
	if (a == CalDetGeo::Y) zfar = layerHeight();
	double zpos = Z0()+2.*ilayer*layerHeight()+zfar;

	double xsize = layerWidth();
	double ysize = layerWidth();
	double zsize = layerHeight();

	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*zsize);

	CalDetGeo layer(ilayer,a,ilayer,P,S);
	return layer;
}
//##############################################
CalDetGeo CalGeometrySvc::getLog(int ilayer, CalDetGeo::axis a, int ilog, idents::ModuleId mod)
//##############################################
{

	double xmod = (mod.ix()-(m_nmodx+1)*0.5)*m_modWidth;
	double ymod = (mod.iy()-(m_nmody+1)*0.5)*m_modWidth;
	Vector modcenter(xmod,ymod,0);
	CalDetGeo log = getLog(ilayer, a, ilog);
	log.setPosition(log.position()+modcenter);
	return log;
}
//##############################################
CalDetGeo CalGeometrySvc::getLog(int ilayer, CalDetGeo::axis a, int ilog)
//##############################################
{
	CalDetGeo layer = getLayer(ilayer,a);
	double xpos = layer.position().x();
	double ypos = layer.position().y();
	double zpos = layer.position().z();
	
	int nlogs = numLogs();
	double posRef = -0.5*(nlogs-1)*(logWidth()+logGap());
	double pos = posRef+ilog*(logWidth()+logGap());
	if (fabs(pos) < 1e-5) pos =0.; 
	if (a == CalDetGeo::X) xpos = pos;
	else ypos = pos;

	double xsize = logLength();
	double ysize = logLength();
	if (a == CalDetGeo::X) xsize = logWidth();
	else ysize = logWidth();
	double zsize = logHeight();
	
	Point P(xpos,ypos,zpos);
	Point S(0.5*xsize,0.5*ysize,0.5*zsize);

	CalDetGeo log(ilayer,a,ilog,P,S);

	return log;
}

// queryInterface
StatusCode  CalGeometrySvc::queryInterface (const IID& riid, void **ppvIF)
{
    if (IID_ICalGeometrySvc == riid) {
        *ppvIF = dynamic_cast<ICalGeometrySvc*> (this);
        return StatusCode::SUCCESS;
    }
    else {
        return Service::queryInterface (riid, ppvIF);
    }
}

// access the type of this service
const IID&  CalGeometrySvc::type () const {
    return IID_ICalGeometrySvc;
}

