
#include "CalRecon/CalADCLogs.h"
//######################################
CalADCLog::CalADCLog(int ilayer, int iview, int icolumn)
:CalLogID(ilayer,iview,icolumn)
//######################################
{
	clear();
}

//######################################
CalADCLog::CalADCLog(int ilayer, int iview, int icolumn, idents::ModuleId mod)
:CalLogID(ilayer,iview,icolumn, mod)
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

	std::cout << logID() << " " << modId() << layer() << view() << column() << " ";
	std::cout << " PED ";
	for (int iside = 0; iside < CALNSIDES; iside++) {
		for (int irange = 0; irange < CALNRANGES; irange++) 
			std::cout << m_ADC[iside][irange]<< " ";
	}
	std::cout<<"\n";
}
//--------------------------------------------------
//######################################
void CalADCLogs::ini(int nModX, int nModY, int nLogs, int nLayers)
//######################################
{
	for (int iy = 1; iy <= nModY; iy++){		
		for (int ix = 1; ix <= nModX; ix++){
			idents::ModuleId mod(ix,iy);
			for (int ilayer = 0; ilayer < nLayers; ilayer++) {
				for (int ilog = 0; ilog < nLogs; ilog++) {
					m_List.push_back(new CalADCLog(ilayer,detGeo::Y,ilog,mod));
					m_List.push_back(new CalADCLog(ilayer,detGeo::X,ilog,mod));
				}
			}
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
