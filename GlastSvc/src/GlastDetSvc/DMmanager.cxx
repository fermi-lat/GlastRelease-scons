// File and Version Information:
// $Header$
//
// Description:
// Wrapper class around detModel, to hide all the initialization and 
// provide defaults.

#include "DMmanager.h"

#include "detModel/Management/Manager.h"
#include "detModel/Management/XercesBuilder.h"
#include "detModel/Sections/Volume.h"
#include "detModel/Management/MaterialsVisitor.h"

#include "idents/VolumeIdentifier.h"

#include "CLHEP/Geometry/Transform3D.h"

#include "xml/IFile.h" //for extractEnvVar
#include <string>
#include <cassert>
#include <iomanip>

DMmanager::DMmanager( )
:m_dm(detModel::Manager::getPointer())
{
}

void DMmanager::init(int argc, char* argv[] ) 
{
    // Purpose and Method:  Set up args: filename and topvol; 
    //   defaults are the standard flight.xml and LAT
    
    std::string filename = (argc<2)? std::string( "-") : std::string(argv[1]);
    std::string topvol   = (argc<3)? std::string("LAT") : std::string(argv[2]);
    std::string mode("fastmc");
    init( filename, mode, topvol);
}

void DMmanager::init(std::string filename, std::string mode, std::string topvol) 
{
    m_dm->setBuilder(new detModel::XercesBuilder);
    if( filename.empty() || filename == "-" ) filename = "$(XMLGEODBSROOT)/xml/flight/flight.xml" ;
    xml::IFile::extractEnvVar(&filename);
    m_dm->setNameFile( filename);
    m_dm->build(detModel::Manager::all);
    m_vol = m_dm->getGdd()->getVolumeByName(topvol);
    if( m_vol==0) {
        std::cerr << "Could not find a topvol " << topvol << std::endl;
        exit(-1);
    }

    // We setup the IDmap with the finest choice mode
    std::string temp = m_dm->getFineChoice();
    m_idMap = new detModel::IDmapBuilder(topvol);
    if (temp!="0")
      m_idMap->setMode(temp);
    m_dm->startVisitor(m_idMap);    
    m_dm->setMode(m_mode=mode);
}

bool DMmanager::getNumericConstByName(std::string name, double* res)
{
  return m_dm->getNumericConstByName(name,res);
}

bool DMmanager::getNumericConstByName(std::string name, int* res)
{
  return m_dm->getNumericConstByName(name,res);
}

void DMmanager::accept( detModel::SectionsVisitor* v){ m_vol->AcceptNotRec(v);  }

void DMmanager::accept( detModel::MaterialsVisitor* v){ m_dm->startVisitor(v);  }

void DMmanager::printSetup(std::ostream& out) { 
    out << " geometry setup: "<< std::endl
        << std::setw(15) << "document = " << m_dm->getNameFile() << std::endl
        << std::setw(15) << "mode = "     << m_mode << std::endl
        << std::setw(15) << "topvol = "   << std::string(m_vol->getName())  << std::endl;

     m_idMap->summary(out);
}

DMmanager::~DMmanager() {        delete m_dm;    }

std::string DMmanager::topvol()const{ return m_vol->getName(); }

idents::VolumeIdentifier DMmanager::getIDPrefix() const
{
  return m_idMap->getIDPrefix();
}

bool DMmanager::getTransform3DByID(idents::VolumeIdentifier id, HepTransform3D* tr)
{
  return m_idMap->getTransform3DByID(id,tr);
}

const HepTransform3D& DMmanager::getTransform3DPrefix() const
{
  return m_idMap->getTransform3DPrefix();
}

bool DMmanager::getShapeByID(idents::VolumeIdentifier id,
                             std::string* s, 
                             std::vector<double>* params){

  return m_idMap->getShapeByID(id,s,params);
}
     
