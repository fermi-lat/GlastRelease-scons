// $Header$

#include "DMmanager.h"

#include "detModel/Management/Manager.h"
#include "detModel/Management/XercesBuilder.h"
#include "detModel/Sections/Volume.h"

#include <string>
#include <cassert>
#include <iomanip>

DMmanager::DMmanager( )
:m_dm(detModel::Manager::getPointer())
{
}

//! initialize with command-line args, supplying defaults
void DMmanager::init(int argc, char* argv[] ) 
{
    // set up args: filename and topvol; defaults are the standard flight.xml and LAT
    //
    std::string filename = (argc<2)? std::string( "-") : std::string(argv[1]);
    std::string topvol   = (argc<3)? std::string("LAT") : std::string(argv[2]);
    std::string mode("fastmc");
    init( filename, mode, topvol);
}

//! Note that if filename is "-" or blank, use xmlUtil default
void DMmanager::init(std::string filename, std::string mode, std::string topvol) 
{
    m_dm->setBuilder(new detModel::XercesBuilder);
    if( filename.empty() || filename == "-" ) filename = std::string(::getenv("XMLUTILROOT"))+"/xml/flight.xml" ;
    ::strncpy(nameBuffer, filename.c_str(), sizeof(nameBuffer));
    m_dm->setNameFile( nameBuffer);
    m_dm->setMode(m_mode=mode);
    m_dm->build(detModel::Manager::all);
    m_vol = m_dm->getGdd()->getVolumeByName(topvol);
    if( m_vol==0) {
        std::cerr << "Could not find a topvol " << topvol << std::endl;
        exit(-1);
    }
    
}

void DMmanager::accept( detModel::SectionsVisitor* v){ m_vol->AcceptNotRec(v);  }

void DMmanager::printSetup(std::ostream& out) { 
    out << " geometry setup: "<< std::endl
        << std::setw(15) << "document = " << m_dm->getNameFile() << std::endl
        << std::setw(15) << "mode = "     << m_mode << std::endl
        << std::setw(15) << "topvol = "   << std::string(m_vol->getName())  << std::endl;
}

DMmanager::~DMmanager() {        delete m_dm;    }

std::string DMmanager::topvol()const{ return m_vol->getName(); }
