
#include "configData/db/LatcDBImplOld.h"
#include "configData/db/TrgConfigDB.h"
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace MOOT;

int main(int argc, char **argv){
  if (argc!=2){
    std::cout<<"Usage: dumpGemConfiguration latc_key"<<std::endl;
    exit(0);
  }
  unsigned key=(unsigned int)atoi(argv[1]);
  std::cout<<"=================================="<<std::endl;
  std::cout<<"GEM configuration for key "<<key<<std::endl;
  std::cout<<"=================================="<<std::endl<<std::endl;
  LatcDBImplOld lc;
  TrgConfigDB tcf(&lc);
  tcf.updateKey(1852);
  tcf.printContrigurator();
  //std::cout<<tcf<<std::endl;
}
