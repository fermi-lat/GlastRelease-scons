
#ifndef CALADCLOGS_H
#define CALADCLOGS_H 1

#include "CalRecon/CalBase.h"
#include "CalRecon/CalLogID.h"
#include "GaudiKernel/DataObject.h"

extern const CLID& CLID_CalADCLogs;

//##############################
class CalADCLog : public CalLogID
//##############################
{
public:

	// construct
	CalADCLog(int ilayer, int iview, int column);
	CalADCLog(int ilayer, int iview, int column, idents::ModuleId mod);
	~CalADCLog() {};
                
	void setNegADC(CalBase::RANGE r, double v)   {m_ADC[CalBase::NEG][r] = v;}
	void setPosADC(CalBase::RANGE r, double v)   {m_ADC[CalBase::POS][r] = v;}
	void addNegADC(CalBase::RANGE r, double v)   {m_ADC[CalBase::NEG][r] += v;}
	void addPosADC(CalBase::RANGE r, double v)   {m_ADC[CalBase::POS][r] += v;}

	// access                        
	double ADC(CalBase::SIDE s, CalBase::RANGE r) const {return m_ADC[s][r];}

	// operations
	void clear();
	void writeOut() const;

private:

	double m_ADC[CALNSIDES][CALNRANGES];
};

//##############################
class CalADCLogs : public DataObject
//##############################
{
public:

	// constructor
	CalADCLogs(int nmodX=1, int nmodY=1, int nlogs=10, int nlayers=4)
	{ini(nmodX,nmodY,nlogs,nlayers);}
	virtual ~CalADCLogs();

	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_CalADCLogs;}
	virtual const CLID& clID() const {return classID();}

	// access
	int num()                 const {return m_List.size();}
	CalADCLog* Log(int i)     const {return m_List[i];}	
	CalADCLog* getLogID(int i) const;

	// operation
	virtual void clear();
	virtual void make() {}
	virtual void writeOut() const;

protected:

	virtual void ini(int nmodX=1, int nmodY=1, int nlogs=10, int nlayers=4);

private:

	std::vector<CalADCLog*> m_List;
};
#endif
