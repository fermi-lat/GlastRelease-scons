
#ifndef CALGEOMETRYSVC_H
#define CALGEOMETRYSVC_H 1

#include "GaudiKernel/Service.h"
#include "CalRecon/ICalGeometrySvc.h"


//----------------------------------------------
//
//   CalGeometrySvc
//
//     It stores the calorimeter geometry.
//
//   A.Chekhtman, NRL, 3/30/01 
//  
//






//##########################################################
class CalGeometrySvc : public Service , virtual public ICalGeometrySvc
//##########################################################
{
public:

    
	//! Constructor of this form must be provided
    CalGeometrySvc(const std::string& name, ISvcLocator* pSvcLocator); 
    virtual ~CalGeometrySvc() {}
    
    StatusCode initialize();
    StatusCode finalize();



	 int numModulesX()    {return m_nmodx;}
	 int numModulesY()    {return m_nmody;}
	 int numViews()     {return m_nviews;}
	 int numLayers()    {return m_nlayers;}
	 int numLogs()      {return m_nLogs;}

	 double moduleWidth()  {return m_modWidth;}
	 double Z0()          {return m_Z0;}
	 double layerWidth()  {return m_layerWidth;}
	 double layerHeight() {return m_layerHeight;}

	 double logWidth()       {return m_logWidth;}
	 double logHeight()      {return m_logHeight;}
	 double logLength()      {return m_logLength;}
	 double logGap()         {return m_logGap;}

	 CalDetGeo getLayer(int ilayer, CalDetGeo::axis a);
	 CalDetGeo getLog(int ilayer, CalDetGeo::axis a, int ilog); 
	 CalDetGeo getLog(int ilayer, CalDetGeo::axis a, int ilog, idents::ModuleId mod);

	     /// queryInterface - for implementing a Service this is necessary
    StatusCode queryInterface(const IID& riid, void** ppvUnknown);

	    /// return the service type
    const IID& type() const;

	static const InterfaceID& interfaceID() { return ICalGeometrySvc::interfaceID(); }



private:

	 int m_geoType;	 

	 int m_nmodx;
	 int m_nmody;
	 int m_nviews;
	 int m_nlayers;
	 int m_nLogs;

	 double m_modWidth;
	
	 double m_Z0;

	 double m_layerWidth;
	 double m_layerHeight;

	 double m_logWidth;
	 double m_logLength;
	 double m_logHeight;
	 double m_logGap;

};

#endif
