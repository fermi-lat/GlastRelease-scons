#ifndef CALCALIBLOGS_H
#define CALCALIBLOGS_H 1

#include <string>
// #include "Event/serviceI.h"
// #include "Event/trsDataVI.h"

#include "TkrRecon/CalBase.h"
#include "TkrRecon/CalLogID.h"

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
void readIntlin(std::istream& file,int side,int range);
void readGain(std::istream& file, int side);
void readRail(std::istream& file, int side);
void readSlope(std::istream& file);
double adc_to_MeV(double adc, int side, int range) const;
private:

	double m_brkpt[CALNSIDES][CALNRANGES];
	double m_coefa[3][CALNSIDES][CALNRANGES];
	double m_coefb[3][CALNSIDES][CALNRANGES];
	double m_gain[CALNSIDES][CALNRANGES];
	double m_rail[CALNSIDES][CALNRANGES];
	double m_slope[CALNRANGES];
};

class CalCalibLogs 
{
public:

	CalCalibLogs() { ini();}
	~CalCalibLogs();

	void setFileNames(const std::string& IntlinFileName,
					  const std::string& GainFileName,
					  const std::string& RailFileName,
					  const std::string& SlopeFileName)
	{

		m_IntlinFile = IntlinFileName;
		m_GainFile   = GainFileName;
		m_RailFile   = RailFileName;
		m_SlopeFile  = SlopeFileName;

	}

	
	void readGain();
	void readIntlin();
	void readRail();
	void readSlope();
	virtual void make();
	virtual void clear() {}
	virtual void ini();

	// access
	int num()                 const {return m_List.size();}
	CalCalibLog* Log(int i)     const {return m_List[i];}	
	CalCalibLog* getLogID(int i) const;

private:

	std::string m_IntlinFile;
	std::string m_GainFile;
	std::string m_RailFile;
	std::string m_SlopeFile;
	std::vector<CalCalibLog*> m_List;

};
#endif
