
#include "CalRecon/CalLogID.h"
//#######################################
CalLogID::CalLogID(int ilayer, detGeo::axis v, int icolumn, idents::ModuleId mod)
//#######################################
:m_modId(mod)
{
	ini();
	m_layer = ilayer;
	m_view = v;
	m_column = icolumn;

	m_logID = CalLogID::ID(ilayer, m_view, icolumn, mod);

}
//#######################################
CalLogID::CalLogID(int ilayer, int iview, int icolumn, idents::ModuleId mod)
//#######################################
:m_modId(mod)
{
	ini();
	m_layer = ilayer;
	m_view = detGeo::makeAxis(iview);
	m_column = icolumn;

	m_logID = CalLogID::ID(ilayer, m_view, icolumn, mod);

}
//#######################################
CalLogID::CalLogID(int ilayer, int iview, int icolumn):m_modId(0)
//#######################################
{
	ini();
	m_layer = ilayer;
	m_view = detGeo::makeAxis(iview);
	m_column = icolumn;

	m_logID = CalLogID::ID(ilayer, m_view, icolumn);
}
//#######################################
CalLogID::CalLogID(int ilayer, detGeo::axis v, int icolumn):m_modId(0)
//#######################################
{
	ini();
	m_layer = ilayer;
	m_view = v;
	m_column = icolumn;

	m_logID = CalLogID::ID(ilayer, m_view, icolumn);
}
//#######################################
int CalLogID::ID(int ilayer, detGeo::axis v, int icolumn, idents::ModuleId mod)
//#######################################
{
	int id = ID(ilayer,v,icolumn)+10000*mod;
	return id;
}
//#######################################
int CalLogID::ID(int ilayer, detGeo::axis v, int icolumn)
//#######################################
{
//	int ifactor = (v == detGeo::X? 1 : 0);
	int ID = 1000*ilayer + 100*v + icolumn;
//	int idLog = 2*nLogs*(3-ilayer)+ifactor*nLogs + icolumn;
	return ID;
}
//#######################################
idents::ModuleId CalLogID::modId(int logID)
//#######################################
{
	idents::ModuleId mod (logID/10000);
	return mod;
}
//#######################################
int CalLogID::layer(int logID)
//#######################################
{
	int i = logID/10000;
	int ilayer = logID/1000 - 10*i;
	return ilayer;
}
//#######################################
detGeo::axis CalLogID::view(int logID)
//#######################################
{
	int i = logID/1000;
	int iview = (logID/100)-10*i;
	detGeo::axis v = detGeo::X;
	if (iview == 1) v = detGeo::Y;
	return v;
}
//#######################################
int CalLogID::column(int logID)
//#######################################
{
	int i = logID/100;
	int icolumn = logID-100*i;
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
