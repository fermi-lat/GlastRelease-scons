#include "CalRecon/CalCalibLogs.h"
// #include "Event/messageManager.h"
#include <fstream>

CalCalibLog::CalCalibLog(int ilayer, int iview, int icolumn)
:CalLogID(ilayer,iview,icolumn){}

//######################################
CalCalibLog* CalCalibLogs::getLogID(int logID) const
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

void CalCalibLogs::ini(int nLogs, int nLayers)
{
	
	for (int ilayer = 0; ilayer < nLayers; ilayer++) {
		for (int ilog = 0; ilog < nLogs; ilog++) {
			m_List.push_back(new CalCalibLog(ilayer,detGeo::X,ilog));
			m_List.push_back(new CalCalibLog(ilayer,detGeo::Y,ilog));
		}
	}
	
	m_GainFile = "";
	m_IntlinFile = "";
	m_RailFile = "";
	m_SlopeFile = "";
}

//######################################
CalCalibLogs::~CalCalibLogs()
//######################################
{
	for (int ilog = 0; ilog < m_List.size(); ilog++) delete m_List[ilog];
}

/*
//######################################
void CalCalibLogs::defineOption()
//######################################
{
	optionVI::defineOption("GainFileName",&m_GainFile);
	optionVI::defineOption("RailFileName",&m_RailFile);
	optionVI::defineOption("IntlinFileName",&m_IntlinFile);
	optionVI::defineOption("SlopeFileName",&m_SlopeFile);
}
*/

void CalCalibLogs::make()
{
	readGain();
	readSlope();
	readIntlin();
	readRail();
}
void CalCalibLogs::readIntlin()
{	
	std::cout << " Calorimeter linearity file : "+m_IntlinFile << std::endl;
	if (m_IntlinFile == "") return;
	
	std::ifstream file;
	file.open(m_IntlinFile.c_str());
	while(file.good()){
		int range;
		int side;
		int col;
		int layer;

		file >> range;
		file >> side;
		file >> col;
		file >> layer;
		if(!file.good())break;

		int view = 1 - layer%2; layer /= 2;
		  side = 1-side; layer = 3-layer;

		CalCalibLog* log = 
		getLogID(CalLogID::ID(layer,detGeo::makeAxis(view),col));
		log->readIntlin(file,side,range);
	}
	file.close();
}

void CalCalibLogs::readGain()
{	
	std::cout << " Calorimeter gain file : "+m_GainFile << std::endl;
	if (m_GainFile == "") return;
	
	std::ifstream file;
	file.open(m_GainFile.c_str());
	while(file.good()){
		int side;
		int col;
		int layer;

		file >> side;
		file >> col;
		file >> layer;
		if(!file.good())break;

		int view = 1 - layer%2; layer /= 2;
		 side = 1-side; layer = 3-layer;

		CalCalibLog* log = 
		getLogID(CalLogID::ID(layer,detGeo::makeAxis(view),col));
		log->readGain(file,side);
	}
	file.close();
}
void CalCalibLogs::readRail()
{	
	std::cout << " Calorimeter rails file : "+m_RailFile << std::endl;
	if (m_RailFile == "") return;
	
	std::ifstream file;
	file.open(m_RailFile.c_str());
	while(file.good()){
		int side;
		int col;
		int layer;

		file >> side;
		file >> col;
		file >> layer;
		if(!file.good())break;

		int view =1-layer%2; layer /= 2;
		 side = 1-side; layer = 3-layer;

		CalCalibLog* log = 
		getLogID(CalLogID::ID(layer,detGeo::makeAxis(view),col));
		log->readRail(file,side);
	}
	file.close();
}
void CalCalibLogs::readSlope()
{	
	std::cout << " Calorimeter  light asymmetry file : "+m_SlopeFile << std::endl;
	if (m_SlopeFile == "") return;
	
	std::ifstream file;
	file.open(m_SlopeFile.c_str());
	while(file.good()){
		int col;
		int layer;

		file >> col;
		file >> layer;
		if(!file.good())break;
		int view = 1-layer%2; layer /= 2;
	    layer = 3-layer;

		CalCalibLog* log = 
		getLogID(CalLogID::ID(layer,detGeo::makeAxis(view),col));
		log->readSlope(file);
	}
	file.close();
}
void CalCalibLog::readIntlin(std::istream& file, int side, int range)
{
	file >> m_brkpt[side][range];
	int i = 0;
	for ( ;i<3;i++)	file >> m_coefa[i][side][range] ;
	for (i=0;i<3;i++)	file >> m_coefb[i][side][range] ;
}

void CalCalibLog::readGain(std::istream& file, int side)
{	const double g[4] = {2.1, 2.1, 8.4, 8.4};
	double gain;
	for ( int i=0;i<CALNRANGES;i++){	file >>gain;
	m_gain[side][i] = gain>0? gain : g[i];
	}
}

void CalCalibLog::readRail(std::istream& file, int side)
{
	for ( int i=0;i<CALNRANGES;i++)	file >> m_rail[side][i];
}

void CalCalibLog::readSlope(std::istream& file)
{
	for ( int i=0;i<CALNRANGES;i++)	file >> m_slope[i];
}

double CalCalibLog::adc_to_MeV(double adc, int s, int r) const
{
	double a[3];
	double bp = getBrkpt(s,r);
	for( int i = 0; i<3; i++)
		a[i] = adc<bp ? getCoefa(i,s,r):getCoefb(i,s,r);
	return (a[0]+(a[1]+a[2]*adc)*adc)*getGain(s,r);
}
