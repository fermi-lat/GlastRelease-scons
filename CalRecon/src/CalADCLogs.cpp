
#include "TkrRecon/CalADCLogs.h"
#include "TkrRecon/calorimeterGeo.h"
// #include "Event/messageManager.h"

//######################################
CalADCLog::CalADCLog(int ilayer, int iview, int icolumn):CalLogID(ilayer,iview,icolumn)
//######################################
{
	clear();
}
//######################################
void CalADCLog::clear()
//######################################
{
	for (int iside = 0 ; iside < CALNSIDES; iside++) {
		for (int irange = 0; irange < CALNRANGES ; irange++) 
			m_ADC[iside][irange] = 0.;
	}
}
//################################################
void CalADCLog::writeOut() const
//################################################
{

	std::cout << logID() << " " << layer() << view() << column() << " ";
	std::cout << " PED ";
	for (int iside = 0; iside < CALNSIDES; iside++) {
		for (int irange = 0; irange < CALNRANGES; irange++) 
			std::cout << m_ADC[iside][irange]<< " ";
	}
	std::cout<<"\n";
}
//--------------------------------------------------
//######################################
void CalADCLogs::ini()
//######################################
{
	int nLayers = calorimeterGeo::numLayers();
	int nLogs   = calorimeterGeo::numLogs();
	for (int ilayer = 0; ilayer < nLayers; ilayer++) {
		for (int ilog = 0; ilog < nLogs; ilog++) {
			m_List.push_back(new CalADCLog(ilayer,detGeo::X,ilog));
			m_List.push_back(new CalADCLog(ilayer,detGeo::Y,ilog));
		}
	}
}
//######################################
CalADCLogs::~CalADCLogs()
//######################################
{
	for (int ilog = 0; ilog < m_List.size(); ilog++) delete m_List[ilog];
}
//######################################
void CalADCLogs::clear()
//######################################
{
	for (int ilog = 0; ilog < m_List.size(); ilog++) {
		m_List[ilog]->clear();
	}
}
//######################################
CalADCLog* CalADCLogs::getLogID(int logID) const
//######################################
{
	int jlog = -1;
	for (int ilog = 0; ilog < m_List.size(); ilog++) {
		if (m_List[ilog]->logID() == logID) {
			jlog  = ilog;
			break;
		}
	}
	return m_List[jlog];
}

//######################################
void CalADCLogs::writeOut() const
//######################################
{
	if (m_List.size()<=0) return;


	std::cout << " CalADCLogs " << m_List.size() <<"\n";
	for (int i = 0; i < m_List.size();i++) {
		Log(i)->writeOut();
	}
}
