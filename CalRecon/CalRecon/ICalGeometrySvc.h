#ifndef ICALGEOMETRYSVC_H
#define ICALGEOMETRYSVC_H 1

#include "CalRecon/CalDetGeo.h"
#include "idents/ModuleId.h"
#include "GaudiKernel/IInterface.h"

//----------------------------------------------
//
//   CalGeometrySvc
//
//     It stores the calorimeter geometry.
//
//   A.Chekhtman, NRL, 3/30/01 
//  
//




static const InterfaceID IID_ICalGeometrySvc(906, 1 , 0);


//##########################################################
class ICalGeometrySvc : virtual public IInterface 
//##########################################################
{
public:




	virtual int numModulesX()=0;
	virtual int numModulesY()=0;
	virtual int numViews()=0;
	virtual int numLayers()=0;
	virtual int numLogs()=0;

	virtual double moduleWidth()=0;
	virtual double Z0()=0;
	virtual double layerWidth()=0;
	virtual double layerHeight()=0;

	virtual double logWidth()=0;
	virtual double logHeight()=0;
	virtual double logLength()=0;
	virtual double logGap() =0;
    virtual double light_att() =0;

	virtual CalDetGeo getLayer(int ilayer, CalDetGeo::axis a) =0;
	virtual CalDetGeo getLog(int ilayer, CalDetGeo::axis a, int ilog)=0; 
	virtual CalDetGeo getLog(int ilayer, CalDetGeo::axis a, int ilog, idents::ModuleId mod) =0;
	static const InterfaceID& interfaceID() { return IID_ICalGeometrySvc; }


};

#endif
