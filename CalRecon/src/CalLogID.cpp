
#include "CalRecon/CalLogID.h"
#include "CalRecon/calorimeterGeo.h"

//#######################################
CalLogID::CalLogID(int ilayer, int iview, int icolumn)
//#######################################
{
	ini();
	m_layer = ilayer;
	m_view = detGeo::makeAxis(iview);
	m_column = icolumn;

	m_logID = CalLogID::ID(ilayer, m_view, icolumn);
}
//#######################################
CalLogID::CalLogID(int ilayer, detGeo::axis v, int icolumn)
//#######################################
{
	ini();
	m_layer = ilayer;
	m_view = v;
	m_column = icolumn;

	m_logID = CalLogID::ID(ilayer, m_view, icolumn);
}

//#######################################
int CalLogID::ID(int ilayer, detGeo::axis v, int icolumn)
//#######################################
{
//	int ifactor = (v == detGeo::X? 1 : 0);
	int ID = 100*ilayer + 10*v + icolumn;
//	int idLog = 2*nLogs*(3-ilayer)+ifactor*nLogs + icolumn;
	return ID;
}
//#######################################
int CalLogID::layer(int logID)
//#######################################
{
	int ilayer = logID/100;
	return ilayer;
}
//#######################################
detGeo::axis CalLogID::view(int logID)
//#######################################
{
	int ilayer = logID/100;
	int iview = (logID/10)-10*ilayer;
	detGeo::axis v = detGeo::X;
	if (iview == 1) v = detGeo::Y;
	return v;
}
//#######################################
int CalLogID::column(int logID)
//#######################################
{
	int ilayer = logID/100;
	int iview =  (logID/10) - 10*ilayer;
	int icolumn = logID-100*ilayer-10*iview;
	return icolumn;
}
//------------ private ------------------
//#######################################
void CalLogID::ini()
//#######################################
{
	m_layer = -1;
	m_view = detGeo::PASSIVE;
	m_column = -1;
	m_logID = -1;
}
