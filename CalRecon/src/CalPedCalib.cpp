
#include "CalRecon/CalPedCalib.h"
// #include "Event/messageManager.h"

#include <fstream>

//######################################
void CalPedCalib::ini()
//######################################
{
	CalADCLogs::ini();

	m_fileName = "";
	clear();
}
//######################################
void CalPedCalib::clear()
//######################################
{
	CalADCLogs::clear();
	for (int iside = 0 ; iside < CALNSIDES; iside++) {
		for (int irange = 0; irange < CALNRANGES ; irange++) 
			m_gain[iside][irange] = 0.;
	}
}
//######################################
void CalPedCalib::make()
//######################################
{
	
	std::cout << " Calorimeter pedestals file : "+m_fileName << std::endl;
	if (m_fileName == "") return;
	
	std::ifstream file ( m_fileName.c_str() );
	std::string name = "";
	while(!file.eof()) {
		file >> name;
		if (name == "gain") readGain(file);
		if (name == "ped") readPedestal(file);
	}
	file.close();
}

//################################################
void CalPedCalib::readGain(std::istream& file)
//################################################
{
	double gain = 0;

	int idLog = 0;
	int ilayer = 0;
	int iview = 0;
	int icolumn = 0;

	file >> idLog;
	file >> ilayer;
	file >> iview;
	file >> icolumn;
	
	file >> gain;
	setPosGain(CalBase::LEX, gain);
	file >> gain;
	setPosGain(CalBase::LE,  gain);
	file >> gain;
	setPosGain(CalBase::HEX, gain);
	file >> gain;
	setPosGain(CalBase::HE,  gain);

	file >> gain;
	setNegGain(CalBase::LEX, gain);
	file >> gain;
	setNegGain(CalBase::LE,  gain);
	file >> gain;
	setNegGain(CalBase::HEX, gain);
	file >> gain;
	setNegGain(CalBase::HE,  gain);
}

//################################################
void CalPedCalib::readPedestal(std::istream& file)
//################################################
{
	int idLog = 0;
	int ilayer = 0;
	int iview = 0;
	int icolumn = 0;

	file >> idLog;
	file >> ilayer;
	file >> iview;
	file >> icolumn;


	// JAH 04/15/00 new definition
 //	if (iview == 0) icolumn = 9-icolumn;



	CalADCLog* log = getLogID(CalLogID::ID(ilayer,detGeo::makeAxis(iview),icolumn));

	double pedLEX;
	double pedLE;
	double pedHEX;
	double pedHE;
	
	pedLEX = 0;
	pedLE  = 0;
	pedHEX = 0;
	pedHE  = 0;

	file >> pedLEX;
	file >> pedLE;
	file >> pedHEX;
	file >> pedHE;

	log->setPosADC(CalBase::LEX, pedLEX);
	log->setPosADC(CalBase::LE,  pedLE);
	log->setPosADC(CalBase::HEX, pedHEX);
	log->setPosADC(CalBase::HE,  pedHE);

	pedLEX = 0;
	pedLE  = 0;
	pedHEX = 0;
	pedHE  = 0;

	file >> pedLEX;
	file >> pedLE;
	file >> pedHEX;
	file >> pedHE;

	log->setNegADC(CalBase::LEX, pedLEX);
	log->setNegADC(CalBase::LE,  pedLE);
	log->setNegADC(CalBase::HEX, pedHEX);
	log->setNegADC(CalBase::HE,  pedHE);

}
/*
//################################################
void CalPedCalib::defineOption()
//################################################
{
	optionVI::defineOption("PedestalFileName",&m_fileName);
}
*/