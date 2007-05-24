
#include "configData/db/LatcDBImplFile.h"
#include "configData/db/TrgConfigDB.h"
#include <vector>
#include <iostream>
#include <stdlib.h>


int main(int argc, char **argv){
  std::string xmlfile=std::string(getenv("CONFIGDATAROOT"))+"/src/defaultTrgConfig.xml";
  std::cout<<"=================================="<<std::endl;
  std::cout<<"File "<<xmlfile<<std::endl;
  std::cout<<"=================================="<<std::endl<<std::endl;
  LatcDBImplFile lc;
  lc.setFilename(xmlfile);
  TrgConfigDB tcf(&lc);
  tcf.updateKey(1);
  tcf.printContrigurator(std::cout);
}
