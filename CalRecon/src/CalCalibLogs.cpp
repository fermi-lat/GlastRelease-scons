#include "CalRecon/CalCalibLogs.h"
#include <fstream>
#include "xml/IFile.h"

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
			m_List.push_back(new CalCalibLog(ilayer,CalDetGeo::X,ilog));
			m_List.push_back(new CalCalibLog(ilayer,CalDetGeo::Y,ilog));
		}
	}
	
	m_GainFile = "";
	m_IntlinFile = "";
	m_RailFile = "";
	m_SlopeFile = "";
    m_MuPeaksFile = "";
    m_ChargePeaksFile = "";
}

//######################################
CalCalibLogs::~CalCalibLogs()
//######################################
{
	for (int ilog = 0; ilog < m_List.size(); ilog++) delete m_List[ilog];
}


void CalCalibLogs::make()
{
	readGain();
	readSlope();
	readIntlin();
	readRail();
    readMuPeaks();
    readChargePeaks();
}

void CalCalibLogs::readMuPeaks()
{
    xml::IFile::extractEnvVar(&m_MuPeaksFile);
	std::cout << " Muon peaks file : "+m_MuPeaksFile << std::endl;
	if (m_MuPeaksFile == "") return;

    std::ifstream file;
    file.open(m_MuPeaksFile.c_str());
    while(file.good()){
        int side;
        int col;
        int layer;
        file >> side >> col >> layer;
        if(!file.good())break;

		int view = 1 - layer%2; layer /= 2;
		  side = 1-side; layer = 3-layer;
		CalCalibLog* log = 
		getLogID(CalLogID::ID(layer,CalDetGeo::makeAxis(view),col));
		log->readMuPeak(file,side);
    }
    file.close();
}

void CalCalibLogs::readChargePeaks()
{	
	xml::IFile::extractEnvVar(&m_ChargePeaksFile);
	std::cout << " Charge peaks file : "+m_ChargePeaksFile << std::endl;
	if (m_IntlinFile == "") return;

	std::ifstream file;
	file.open(m_ChargePeaksFile.c_str());
	while(file.good()){
		int range;
		int side;
		int col;
		int layer;
        int point;

		file >> side;
		file >> col;
		file >> layer;
		file >> range;
        file >> point;
		if(!file.good())break;

		int view = 1 - layer%2; layer /= 2;
		  side = 1-side; layer = 3-layer;
 //         if(view == 1) col = 9-col;

		CalCalibLog* log = 
		getLogID(CalLogID::ID(layer,CalDetGeo::makeAxis(view),col));
		log->readChargePeak(file,side,range);
	}
	file.close();
}
void CalCalibLogs::readIntlin()
{	
	xml::IFile::extractEnvVar(&m_IntlinFile);
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
		getLogID(CalLogID::ID(layer,CalDetGeo::makeAxis(view),col));
		log->readIntlin(file,side,range);
	}
	file.close();
}

void CalCalibLogs::readGain()
{	
	xml::IFile::extractEnvVar(&m_GainFile);
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
		getLogID(CalLogID::ID(layer,CalDetGeo::makeAxis(view),col));
		log->readGain(file,side);
	}
	file.close();
}
void CalCalibLogs::readRail()
{	
	xml::IFile::extractEnvVar(&m_RailFile);
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
		getLogID(CalLogID::ID(layer,CalDetGeo::makeAxis(view),col));
		log->readRail(file,side);
	}
	file.close();
}
void CalCalibLogs::readSlope()
{	
	xml::IFile::extractEnvVar(&m_SlopeFile);
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
		getLogID(CalLogID::ID(layer,CalDetGeo::makeAxis(view),col));
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
    double slope;
    file >> slope >> m_offset;
    for ( int i=0;i<CALNRANGES;i++) m_slope[i] = slope;
}

void CalCalibLog::readMuPeak(std::istream& file, int side)
{
		file >> m_mu_peak[side];
}
void CalCalibLog::readChargePeak(std::istream& file, int side,int range)
{
    float dac,adc;
    int good;
    file >> dac;
    file >> adc;
    file >> good;
    if(good == 1 ){
        std::map<float,float>& mcp = m_charge_peaks[side][range];
        mcp[4095-adc] = dac;
    }
}

double CalCalibLog::adc_to_MeV(double adc, int s, int r) const
{
	double a[3];
	double bp = getBrkpt(s,r);
	for( int i = 0; i<3; i++)
		a[i] = adc<bp ? getCoefa(i,s,r):getCoefb(i,s,r);
	return (a[0]+(a[1]+a[2]*adc)*adc)*getGain(s,r);
}

double CalCalibLog::adc_to_dac(double adc, int s, int r) const
{
    const std::map<float,float>& mcp = getChargePeaks(s,r);
    std::map<float,float>::const_iterator adc0 = mcp.lower_bound(adc);
    if(adc0 == mcp.begin())adc0++;
    std::map<float,float>::const_iterator adc1 = adc0--;
    if((*adc0).second == 0){adc0++;adc1++;}
    float a0 = (*adc0).first;
    float a1 = (*adc1).first;
    float d0 = (*adc0).second;
    float d1 = (*adc1).second;
    float ped = (*(mcp.begin())).first;
    float d = d0 + (adc-a0)*(d1-d0)/(a1-a0);
    return (d<2)? 0 : d;
}

double CalCalibLog::dac_to_MeV(double dac, int s, int r) const
{
     double adc_mu = getMuPeak(s);
     double offset = getOffset();
     double dac_mu = adc_to_dac(adc_mu,s,0)*(1.0-offset*(1-2*s));
     if(dac_mu < 1.0) dac_mu = 10.0;
     double ene = 14.0*dac/dac_mu;
     return (r < 2) ? ene : 60*ene; 
}