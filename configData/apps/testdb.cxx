
#include "mootCore/MootQuery.h"

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

using namespace MOOT;

int main(int argc, char **argv){
  unsigned key=(unsigned int)atoi(argv[1]);
  std::vector<FswInfo> info;
  MootQuery mq;
  mq.getConfigFswInfo(key,info);
  std::cout<<"Key "<<key<<std::endl;
  for (unsigned i=0;i<info.size();i++){
    std::cout<<info[i].getSrc()<<" "<<info[i].getFmxKey()<<std::endl;
  }
}
