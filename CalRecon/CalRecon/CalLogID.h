
#ifndef CALLOGID_H
#define CALLOGID_H 1

#include "TkrRecon/detGeo.h"

//############################
class CalLogID
//############################
{
public:
	// construct
	CalLogID() {ini();}
	CalLogID(int ilayer, int iview, int icolumn);
	CalLogID(int ilayer, detGeo::axis v, int icolumn);
	~CalLogID() {}

	// access
	int layer()         const {return m_layer;}
	detGeo::axis view() const {return m_view;}
	int column()        const {return m_column;}
	int logID()         const {return m_logID;}

	// operations
	static int ID(int ilayer, detGeo::axis v, int icolumn);
	static int layer(int logID);
	static detGeo::axis view(int logID);
	static int column(int logID);

private:

	void ini();

private:

	int m_layer;
	detGeo::axis m_view;
	int m_column;
	int m_logID;
};

#endif
