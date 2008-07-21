// do something useful another day
#include "AncillaryDataUtil/AncillaryGeometry.h"
#include "facilities/commonUtilities.h"
#include <iostream>
void TestGeometry()
{
  std::string m_geometryFilePath="$(ANCILLARYDATAUTILDATAPATH)/Geometry_v0.dat";
  std::string m_rcReportSample="$(ANCILLARYDATAUTILDATAPATH)/rcReport.out";
  std::cout<< "loading geometry from " << m_geometryFilePath << std::endl;
  std::cout<< "reading rcReport from " << m_rcReportSample << std::endl;
  AncillaryData::AncillaryGeometry   *m_geometry = new AncillaryData::AncillaryGeometry(m_geometryFilePath,m_rcReportSample);
  m_geometry->print();
}



int main(){
  facilities::commonUtilities::setupEnvironment();
  TestGeometry();
  return 0;
}


