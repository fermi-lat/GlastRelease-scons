
#ifndef CALLOGID_H
#define CALLOGID_H 1

#include "TkrRecon/detGeo.h"
#include "idents/ModuleId.h"
//############################
class CalLogID
//############################
{
public:
	// construct
	CalLogID():m_modId(0) {ini();}
	CalLogID(int ilayer, int iview, int icolumn);
	CalLogID(int ilayer, int iview, int icolumn, idents::ModuleId mod);
	CalLogID(int ilayer, detGeo::axis v, int icolumn);
	CalLogID(int ilayer, detGeo::axis v, int icolumn, idents::ModuleId mod);
	~CalLogID() {}

	// access
	int layer()         const {return m_layer;}
	detGeo::axis view() const {return m_view;}
	int column()        const {return m_column;}
	int logID()         const {return m_logID;}
	idents::ModuleId modId() const {return m_modId;}	

	// operations
	static int ID(int ilayer, detGeo::axis v, int icolumn);
	static int ID(int ilayer, detGeo::axis v, int icolumn, idents::ModuleId mod);
	static int layer(int logID);
	static detGeo::axis view(int logID);
	static int column(int logID);
	static idents::ModuleId modId(int logID);

private:

	void ini();

private:

	int m_layer;
	detGeo::axis m_view;
	int m_column;
	int m_logID;
	idents::ModuleId m_modId;
};

#endif
