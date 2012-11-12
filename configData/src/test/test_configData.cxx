
#include "configData/db/LatcDBImplFile.h"
#include "configData/db/TrgConfigDB.h"
#include "facilities/commonUtilities.h"
#include <vector>
#include <iostream>
#include <stdlib.h>


int main(int argc, char **argv){
  // std::string xmlfile=std::string(getenv("CONFIGDATAROOT"))+"/src/defaultTrgConfig.xml";
  facilities::commonUtilities::setupEnvironment();
  std::string xmlfile=
    facilities::commonUtilities::joinPath(facilities::commonUtilities::getXmlPath("configData"), "defaultTrgConfig.xml");
  std::cout<<"=================================="<<std::endl;
  std::cout<<"File "<<xmlfile<<std::endl;
  std::cout<<"=================================="<<std::endl<<std::endl;
  LatcDBImplFile lc;
  lc.setFilename(xmlfile);
  TrgConfigDB tcf(&lc);
  tcf.updateKey(1);
  tcf.printContrigurator(std::cout);
  std::cout << std::endl << "Test program completed successfully" << std::endl;
  exit(0);
}
