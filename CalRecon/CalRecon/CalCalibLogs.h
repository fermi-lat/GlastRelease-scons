#ifndef CALCALIBLOGS_H
#define CALCALIBLOGS_H 1

#include <string>
#include <map>
// #include "Event/serviceI.h"
// #include "Event/trsDataVI.h"

#include "CalRecon/CalBase.h"
#include "CalRecon/CalLogID.h"

class CalCalibLog : public CalLogID
{
public:

//constructor
CalCalibLog(int ilayer, int iview, int column);
~CalCalibLog() {}

double getBrkpt( int side, int range)const {return m_brkpt[side][range];}
double getCoefa( int i,int side,int range)const {return m_coefa[i][side][range];}
double getCoefb( int i,int side,int range)const {return m_coefb[i][side][range];}
double getGain(int side,int range)const {return m_gain[side][range];}
double getRail(int side,int range)const {return m_rail[side][range];}
double getSlope(int range)const {return m_slope[range];}
double getOffset()const {return m_offset;}
double getMuPeak(int side)const {return m_mu_peak[side];}
const std::map<float,float>& getChargePeaks(int side, int range)const 
                                    {return m_charge_peaks[side][range];}
void readIntlin(std::istream& file,int side,int range);
void readGain(std::istream& file, int side);
void readRail(std::istream& file, int side);
void readSlope(std::istream& file);
void readChargePeak(std::istream& file, int side, int range);
void readMuPeak(std::istream& file, int side );
double adc_to_MeV(double adc, int side, int range) const;
double adc_to_dac(double adc, int side, int range) const;
double dac_to_MeV(double dac, int side, int range) const;
private:

	double m_brkpt[CALNSIDES][CALNRANGES];
	double m_coefa[3][CALNSIDES][CALNRANGES];
	double m_coefb[3][CALNSIDES][CALNRANGES];
	double m_gain[CALNSIDES][CALNRANGES];
	double m_rail[CALNSIDES][CALNRANGES];
	double m_slope[CALNRANGES];
    double m_offset;
    std::map<float,float> m_charge_peaks[CALNSIDES][CALNRANGES];
    double m_mu_peak[CALNSIDES];
};

class CalCalibLogs 
{
public:

	CalCalibLogs(int nLogs = 10, int nLayers = 4) { ini(nLogs,nLayers);}
	~CalCalibLogs();

	void setFileNames(const std::string& IntlinFileName,
					  const std::string& GainFileName,
					  const std::string& RailFileName,
					  const std::string& SlopeFileName,
                      const std::string& ChargePeaksFileName,
                      const std::string& MuPeaksFileName)
	{

		m_IntlinFile = IntlinFileName;
		m_GainFile   = GainFileName;
		m_RailFile   = RailFileName;
		m_SlopeFile  = SlopeFileName;
        m_ChargePeaksFile = ChargePeaksFileName;
        m_MuPeaksFile = MuPeaksFileName;

	}

	
	void readGain();
	void readIntlin();
	void readRail();
	void readSlope();
    void readMuPeaks();
    void readChargePeaks();
	virtual void make();
	virtual void clear() {}
	virtual void ini(int nLogs, int nLayers);

	// access
	int num()                 const {return m_List.size();}
	CalCalibLog* Log(int i)     const {return m_List[i];}	
	CalCalibLog* getLogID(int i) const;

private:

	std::string m_IntlinFile;
	std::string m_GainFile;
	std::string m_RailFile;
	std::string m_SlopeFile;
    std::string m_MuPeaksFile;
    std::string m_ChargePeaksFile;
	std::vector<CalCalibLog*> m_List;

};
#endif
