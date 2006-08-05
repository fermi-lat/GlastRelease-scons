// do something useful another day
#include "AncillaryDataUtil/AncillaryGeometry.h"
#include "facilities/Util.h"
#include <iostream>
void TestGeometry()
{
  std::string m_geometryFilePath="$(ANCILLARYDATAUTILROOT)/data/Geometry_v0.dat";
  facilities::Util::expandEnvVar(&m_geometryFilePath); 
  std::cout<< "loading geometry from " << m_geometryFilePath << std::endl;
  AncillaryData::AncillaryGeometry   *m_geometry = new AncillaryData::AncillaryGeometry(m_geometryFilePath);
  m_geometry->print();
}



int main(){
  
  TestGeometry();
  return 0;
}


