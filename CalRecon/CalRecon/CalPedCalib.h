#ifndef CALPEDCALIB_H
#define CALPEDCALIB_H

#include <string>
// #include "Event/serviceI.h"
// #include "Event/trsDataVI.h"
#include "TkrRecon/CalBase.h"
#include "TkrRecon/CalADCLogs.h"


//###########################################
class CalPedCalib : public CalADCLogs
//###########################################
{
public:

	// construct
	CalPedCalib() {ini(); }
	~CalPedCalib() {} 

	void setFileName(const std::string& filename) { m_fileName = filename;}

	//! define options
//	virtual void defineOption();

	void setNegGain(CalBase::RANGE r, double v) {m_gain[CalBase::NEG][r] = v;}
	void setPosGain(CalBase::RANGE r, double v) {m_gain[CalBase::POS][r] = v;}

	// access
	double gain(CalBase::SIDE s, CalBase::RANGE r) const {return m_gain[s][r];}

	// operation
	virtual void clear();
	virtual void make();


protected:

	virtual void ini();

private:

	void readGain(std::istream& file);
	void readPedestal(std::istream& file);

private:

	std::string m_fileName;
	double m_gain[CALNSIDES][CALNRANGES];
};
#endif
