
#include "configData/db/LatcDBImpl.h"
#include "configData/db/LatcDBImplOld.h"
#include "configData/db/LatcDBImplFile.h"
#include "configData/db/TrgConfigDB.h"
#include <vector>
#include <iostream>
#include <stdlib.h>
#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
#endif

using namespace MOOT;

int main(int argc, char **argv){
  bool fullconfig=false;
  bool allowMissing=false;
  int opt;
  int argn=1;
#ifdef WIN32
  while ( (opt = facilities::getopt(argc, argv, "fm")) != EOF ) {
#else
  while ( (opt = getopt(argc, argv, "fm")) != EOF ) {
#endif
    argn++;
    switch(opt){
    case 'f':
      fullconfig=true;
      break;
    case 'm':
      allowMissing=true;
      std::cout<<"Warning: Allowing partial configuration"<<std::endl;
      break;
    }
  }
  if (argc!=argn+1  ){
    std::cout<<"Usage: dumpGemConfiguration -f [full config] -m [allow partial config] latc_key or xml_file"<<std::endl;
    exit(0);
  }
  char *endptr;
  //unsigned key=(unsigned int)atoi(argv[1]);
  unsigned key=strtoul(argv[argn],&endptr,0);
  if (endptr[0]=='\0'){
    std::cout<<"=================================="<<std::endl;
    std::cout<<"GEM configuration for key "<<key<<std::endl;
    std::cout<<"=================================="<<std::endl<<std::endl;
    LatcDBImpl lc;
    TrgConfigDB tcf(&lc);
    tcf.allowMissing(allowMissing);
    tcf.updateKey(key);
    if (!fullconfig)tcf.printContrigurator(std::cout);
    else std::cout<<tcf<<std::endl;
  }else{
    std::cout<<"=================================="<<std::endl;
    std::cout<<"GEM configuration for file "<<argv[argn]<<std::endl;
    std::cout<<"=================================="<<std::endl<<std::endl;
    LatcDBImplFile lc;
    lc.setFilename(std::string(argv[argn]));
    TrgConfigDB tcf(&lc);
    tcf.allowMissing(allowMissing);
    tcf.updateKey(1);
    if (!fullconfig)tcf.printContrigurator(std::cout);
    else std::cout<<tcf<<std::endl;
  }
}
