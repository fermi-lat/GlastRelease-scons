
#include "CalRecon/CalRecLogs.h"
// #include "Event/messageManager.h"
#include "CalRecon/calorimeterGeo.h"
// #include "Event/dataManager.h"
// #include "reconstruction/CsIClusters.h"

//----------------- CalRecLog ------------------

//################################################
CalRecLog::CalRecLog(int ilayer, detGeo::axis v, int icolumn) : 
CalADCLog(ilayer,v,icolumn)
//################################################
{
	clear();
}
//################################################
double CalRecLog::asymmetry() const
//################################################
{
	double asymmetry = 0.;
	if (m_pEnergy > 0 && m_nEnergy > 0)
		asymmetry = (m_pEnergy-m_nEnergy)/(m_pEnergy+m_nEnergy);
	return asymmetry;
}
//################################################
void CalRecLog::writeOut() const
//################################################
{

	std::cout << " ID " << logID() <<  "-" << layer() << view() << column() << " ";
	std::cout << " E-+ " << energy() << " " << negEnergy() << " " << posEnergy() ;
	std::cout << " Asy " << asymmetry();
	std::cout << " Pos " << position().x() << " " << position().y() << " " << position().z();
	std::cout << "\n";
}
//######################################################
void CalRecLog::draw(gui::DisplayRep& v) const
//######################################################
{
	double MINENE = 10.;
	double FACTOR = 0.1;
	double delta = 0.2*calorimeterGeo::logWidth();
	double x = position().x();
	double y = position().y();
	double z = position().z();
	double ene = energy();
	if (ene <= MINENE) return;

	v.setColor("black");
	v.moveTo(Point(x-delta,y,z));
	v.lineTo(Point(x+delta,y,z));
	v.moveTo(Point(x,y-delta,z));
	v.lineTo(Point(x,y+delta,z));

	v.moveTo(Point(x,y,z));
	v.lineTo(Point(x,y,z-log(FACTOR*ene)));
}
//------------- private --------------------------
//################################################
void CalRecLog::clear()
//################################################
{
	CalADCLog::clear();

	for (int irange = 0; irange < CALNRANGES; irange++) {
		m_negEnergy[irange] = 0;
		m_posEnergy[irange] = 0;
	}

	m_nEnergy = 0.;
	m_pEnergy = 0.;
	m_position = Point(0.,0.,0.);
}

//------------ private ---------------------------
//######################################
void CalRecLogs::ini()
//######################################
{
	int nLayers = calorimeterGeo::numLayers();
	int nLogs   = calorimeterGeo::numLogs();
	for (int ilayer = 0; ilayer < nLayers; ilayer++) {
		for (int ilog = 0; ilog < nLogs; ilog++) {
			m_List.push_back(new CalRecLog(ilayer,detGeo::X,ilog));
			m_List.push_back(new CalRecLog(ilayer,detGeo::Y,ilog));
		}
	}
}
//######################################
CalRecLogs::~CalRecLogs()
//######################################
{
	for (int ilog = 0; ilog < m_List.size(); ilog++) delete m_List[ilog];
}
//######################################
void CalRecLogs::clear()
//######################################
{
	for (int ilog = 0; ilog < m_List.size(); ilog++) {
		m_List[ilog]->clear();
	}
}
//######################################
CalRecLog* CalRecLogs::getLogID(int logID) const
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
void CalRecLogs::writeOut() const
//######################################
{
	if (m_List.size()<=0) return;

	std::cout << " CalRecLogs " << m_List.size() <<"\n";
	for (int i = 0; i < m_List.size();i++) {
		m_List[i]->writeOut();
	}
}
//######################################
void CalRecLogs::draw(gui::DisplayRep& v) const
//######################################
{
	if (m_List.size() <=0) return;
	for (int i = 0; i < m_List.size();i++) m_List[i]->draw(v);
}
