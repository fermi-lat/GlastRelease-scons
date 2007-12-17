
#include "configData/db/LatcDBImpl.h"
#include "configData/db/LatcDBImplOld.h"
#include "configData/db/LatcDBImplFile.h"
#include "configData/db/TrgConfigDB.h"
#include <vector>
#include <iostream>
#include <fstream>
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
  std::ostream *out=&std::cout;
#ifdef WIN32
  while ( (opt = facilities::getopt(argc, argv, "fmo:")) != EOF ) {
#else
  while ( (opt = getopt(argc, argv, "fmo:")) != EOF ) {
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
    case 'o':
      out=new std::ofstream(optarg);
      argn++;
      break;
    }
  }
  if (argc!=argn+1  ){
    (*out)<<"Usage: dumpGemConfiguration -f [full config] -m [allow partial config] -o [filename] latc_key or xml_file"<<std::endl;
    exit(0);
  }
  char *endptr;
  //unsigned key=(unsigned int)atoi(argv[1]);
  unsigned key=strtoul(argv[argn],&endptr,0);
  if (endptr[0]=='\0'){
    (*out)<<"=================================="<<std::endl;
    (*out)<<"GEM configuration for key "<<key<<std::endl;
    (*out)<<"=================================="<<std::endl<<std::endl;
    LatcDBImpl lc;
    TrgConfigDB tcf(&lc);
    tcf.allowMissing(allowMissing);
    tcf.updateKey(key);
    if (!fullconfig)tcf.printContrigurator((*out));
    else (*out)<<tcf<<std::endl;
  }else{
    (*out)<<"=================================="<<std::endl;
    (*out)<<"GEM configuration for file "<<argv[argn]<<std::endl;
    (*out)<<"=================================="<<std::endl<<std::endl;
    LatcDBImplFile lc;
    lc.setFilename(std::string(argv[argn]));
    TrgConfigDB tcf(&lc);
    tcf.allowMissing(allowMissing);
    tcf.updateKey(1);
    if (!fullconfig)tcf.printContrigurator((*out));
    else (*out)<<tcf<<std::endl;
  }
}
