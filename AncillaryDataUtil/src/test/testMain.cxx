// do something useful another day
#include "AncillaryDataUtil/AncillaryGeometry.h"
#include <iostream>
void TestGeometry()
{
  std::string m_geometryFilePath="$(ANCILLARYDATAUTILROOT)/data/Geometry_v0.dat";
  std::string m_rcReportSample="$(ANCILLARYDATAUTILROOT)/data/rcReport.out";
  std::cout<< "loading geometry from " << m_geometryFilePath << std::endl;
  std::cout<< "reading rcReport from " << m_rcReportSample << std::endl;
  AncillaryData::AncillaryGeometry   *m_geometry = new AncillaryData::AncillaryGeometry(m_geometryFilePath,m_rcReportSample);
  m_geometry->print();
}



int main(){
  
  TestGeometry();
  return 0;
}


