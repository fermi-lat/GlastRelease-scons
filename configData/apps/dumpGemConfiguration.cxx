
#include "configData/db/LatcDBImpl.h"
#include "configData/db/MootDBImpl.h"
#include "configData/db/LatcDBImplOld.h"
#include "configData/db/LatcDBImplFile.h"
#include "configData/db/TrgConfigDB.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cassert>
#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
using facilities::getopt;
#endif

using namespace MOOT;

int main(int argc, char **argv){
  bool fullconfig=false;
  bool allowMissing=false;
  bool mootKey=false;
  bool roi=false;
  bool gem=false;
  int opt;
  int argn=1;
  std::ostream *out=&std::cout;

  while ( (opt = getopt(argc, argv, "fMmrgo:")) != EOF ) {
    argn++;
    switch(opt){
    case 'f':
      fullconfig=true;
      break;
    case 'm':
      allowMissing=true;
      std::cout<<"Warning: Allowing partial configuration"<<std::endl;
      break;
    case 'M': // mootKey
      mootKey=true;
      break;
    case 'r': // print ROI config only
      roi=true;
      break;
    case 'g': // print GEM config only
      gem=true;
      break;
    case 'o':
      out=new std::ofstream(optarg);
      argn++;
      break;
    }
  }
  if (argc!=argn+1  ){
    std::cout << "jhpdebug " << argc << std::endl;
    (*out)<<"Usage: dumpGemConfiguration -f [full config] -m [allow partial config] -M [MOOT key] -o [filename] latc_key or xml_file"<<std::endl;
    exit(0);
  }
  char *endptr;
  //unsigned key=(unsigned int)atoi(argv[1]);
  unsigned key=strtoul(argv[argn],&endptr,0);
  TrgConfigDB *tcf;
  if (endptr[0]=='\0'){
    (*out)<<"=================================="<<std::endl;
    (*out)<<"GEM configuration for key "<<key<<std::endl;
    (*out)<<"=================================="<<std::endl<<std::endl;
    LatcDB *db;
    if (mootKey) {
      db = new MootDBImpl;
    } else {
      db = new LatcDBImpl;
    }
    tcf=new TrgConfigDB (db);
  }else{
    (*out)<<"=================================="<<std::endl;
    (*out)<<"GEM configuration for file "<<argv[argn]<<std::endl;
    (*out)<<"=================================="<<std::endl<<std::endl;
    LatcDBImplFile *lc=new LatcDBImplFile;
    lc->setFilename(std::string(argv[argn]));
    tcf=new TrgConfigDB(lc);
    key=1;
  }
  tcf->allowMissing(allowMissing);
  tcf->updateKey(key);
  if (!fullconfig && !roi)tcf->printContrigurator((*out));
  else if(fullconfig && ((!roi && !gem) || (roi&&gem)))(*out)<<(*tcf)<<std::endl;
  else if (roi && !gem)(*out)<<(*tcf->roi())<<std::endl;
  else if (gem && !roi){
    (*out)<<"GEM configuration: "<<std::endl;
    (*out)<<(*tcf->configuration())<<std::endl;
    (*out)<<(*tcf->windowParams())<<std::endl;
    (*out)<<(*tcf->periodicTrigger())<<std::endl;
    (*out)<<(*tcf->lut())<<std::endl;
    (*out)<<(*tcf->trgEngine())<<std::endl;
    (*out)<<(*tcf->disabledChannels())<<std::endl;
  }else{
    std::cout<<"Bad combination of options -f, -r, -g"<<std::endl;
    assert(0);
  }
}
