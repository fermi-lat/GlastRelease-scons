#ifndef __CALRECLOGS_H
#define __CALRECLOGS_H 1

//#include <stdlib.h>
#include <iostream>
#include <vector>
// #include "Event/trsDataVI.h"
#include "TkrRecon/detGeo.h"
#include "TkrRecon/CalADCLogs.h"
#include "TkrRecon/CalBase.h"
#include "Gaudi/Kernel/DataObject.h"
#include "gui/DisplayRep.h"


extern const CLID& CLID_CalRecLogs;


//----------------------------------------------
//
//   CalRecLogs
//
//   Transient Storage Data
//----------------------------------------------
//   It contains the Rec data for tha calorimeter logs
//----------------------------------------------
//             J.A Hernando, Santa Cruz 02/29/00
//----------------------------------------------





//##############################################
class CalRecLog : public CalADCLog
//##############################################
{
public:

	//construct
	CalRecLog(int ilayer, detGeo::axis iview, int ilog);
	void setNegEnergy(CalBase::RANGE r, double e)  {m_negEnergy[r] = e;}
	void setPosEnergy(CalBase::RANGE r, double e)  {m_posEnergy[r] = e;}
	void setNegEnergy(double e)                    {m_nEnergy = e;}
	void setPosEnergy(double e)                    {m_pEnergy = e;}
	void setPosition(Point p)    {m_position = p;}
	void setBestRange(CalBase::RANGE r) {m_bestRange = r;}
	~CalRecLog() {};
	
	double negEnergy(CalBase::RANGE r)   const {return m_negEnergy[r];}
	double posEnergy(CalBase::RANGE r)   const {return m_posEnergy[r];}
	double negEnergy()                   const {return m_nEnergy;}
	double posEnergy()                   const {return m_pEnergy;}
	double energy()      const {return 0.5*(posEnergy()+negEnergy());}
	double energy(CalBase::RANGE r)      const {return 0.5*(negEnergy(r)+posEnergy(r));}
	double asymmetry()   const;
	Point  position()    const {return m_position;}
	CalBase::RANGE bestRange() const {return m_bestRange;}
	
	// operations
	void clear();
	void writeOut() const;
	void draw(gui::DisplayRep& v) const;

private:

	double m_negEnergy[CALNRANGES];
	double m_posEnergy[CALNRANGES];
	double m_nEnergy;
	double m_pEnergy;
	Point m_position;
	CalBase::RANGE m_bestRange;

};

//##############################################
class CalRecLogs : public DataObject
//##############################################
{
public:

	// constructor
	CalRecLogs()  {ini();}
	~CalRecLogs();


	// GAUDI members to be use by the converters
	static const CLID& classID() {return CLID_CalRecLogs;}
	virtual const CLID& clID() const {return classID();}
	
	
	
	
	void add(CalRecLog* log)   {m_List.push_back(log);}
	
	// access
	int num()                const {return m_List.size();}
	CalRecLog* Log(int i)    const {return m_List[i];}
	CalRecLog* getLogID(int logID) const;

	// operations
	virtual void clear();
	virtual void make() {};

	void writeOut() const;
	void update(gui::DisplayRep& v)  {draw(v);}

private:

	virtual void ini();

	void draw(gui::DisplayRep& v) const;

private:

	std::vector<CalRecLog*>  m_List;

};

#endif
