
//
// stl
#include <cstdio>
#include <string>
#include <iostream>
#include "facilities/Util.h"
#ifndef WIN32
#include <unistd.h>
#else
#include "facilities/XGetopt.h"
using facilities::getopt;
#endif
//

// local
#include "configData/db/MootVoteParser.h"

void usage() {
  std::cout << "findConfigs.exe [options]" << std::endl
	    << "\t\tOptions:" << std::endl
	    << "\t-p <Precinct:Alias>  \tMatch Precinct and Alias" << std::endl
	    << "\t-a <ancil_key>       \tMatch ancillary file key" << std::endl
	    << "\t-v <Precinct:VoteTag>\tMatch XML tag in Vote file" << std::endl
	    << "\t-O                   \tOR mode match, any parts of query may match" << std::endl
	    << "\t-b                   \tBrief printout" << std::endl
	    << std::endl;
}

bool makeVoteAlias(char* fromString, std::list<configData::VoteAlias>& toList) {
  std::vector<std::string> tokens;
  facilities::Util::stringTokenize(fromString,":",tokens);
  if ( tokens.size() != 2 ) return false;
  configData::PrecinctType pt = configData::MootVoteParser::precinctType(tokens[0]);
  if ( pt == configData::Unknown ) {
    return false;
  }
  toList.push_back(configData::VoteAlias(pt,tokens[1]));
  return true;
}

bool splitToken(const std::string& str, std::string& name, std::string& value) {
  std::vector<std::string> tokens;
  facilities::Util::stringTokenize(str,"=",tokens);
  switch ( tokens.size() ) {
  case 2:
    value = tokens[1];
    // no break
  case 1:
    name = tokens[0];
    break;
  default:
    return false;
  }
  return true;
}

bool makeVoteTag(char* fromString, std::list<configData::VoteTag>& toList) {
  std::vector<std::string> tokens;
  facilities::Util::stringTokenize(fromString,":",tokens);
  if ( tokens.size() != 2 ) return false;
  configData::PrecinctType pt = configData::MootVoteParser::precinctType(tokens[0]);
  if ( pt == configData::Unknown ) {
    return false;
  }
  std::vector<std::string> tags;
  facilities::Util::stringTokenize(tokens[1],",",tags);
  if ( tags.size() == 0 ) return false;
  std::string name;
  std::string value;
  if ( ! splitToken(tags[0],name,value) ) return false;
  toList.push_back(configData::VoteTag(pt,name,value));
  configData::VoteTag& theTag = toList.back();
  for ( unsigned i(1); i < tags.size(); i++ ) {
    if ( ! splitToken(tags[i],name,value) ) return false;
    theTag.m_tagAttr[name] = value;
  }
  return true;
}

bool makeAncilKey(char* fromString, std::list<unsigned>& toList) {
  try { 
    unsigned val = facilities::Util::stringToUnsigned(fromString);
    toList.push_back(val);
  } catch ( std::exception& e ) {
    return false;
  }
  return true;
}


//
int main(int argn, char** argc) {
  
  // parse options
  //char* endPtr;  
  std::list<configData::VoteAlias> vaList;
  std::list<configData::VoteTag>   vtList;
  std::list<unsigned>  akList;
  bool matchAll(true);
  bool brief(false);

  int opt;
  while ( (opt = getopt(argn, argc, "hObp:a:v:")) != EOF ) {
    switch (opt) {
    case 'h':   // help      
      usage();
      return 1;
    case 'O':
      matchAll = false;
      break;
    case 'b':
      brief = true;
      break;
    case 'p':
      if ( ! makeVoteAlias( optarg, vaList ) ) {
	usage();
	return 2;
      }
      break;
    case 'a':
      if ( ! makeAncilKey( optarg, akList ) ) {
	usage();
	return 2;
      }
      break;
    case 'v':   //  output
      if ( ! makeVoteTag( optarg, vtList ) ) {
	usage();
	return 2;
      }
      break;
    default:
      std::cout << opt << " not parsable..." << std::endl;
      usage();
      return 3;  
    }
  }

  if ( optind != argn ) {
    usage();
    return 3;  
  }

  configData::MootVoteParser mvp;
  std::set<unsigned> keys;
  bool started(false);
  if ( vaList.size() != 0 ) {
    std::set<unsigned> tmpKeys;
    if ( ! mvp.getMasterKeysWithVoteAliases(tmpKeys,vaList,matchAll) ) return 4;
    if ( ! started ) {
      keys = tmpKeys;
      started = true;
    } else {
      if ( matchAll ) {
	mvp.mergeSetsByISect(keys,tmpKeys);
      } else {
	mvp.mergeSetsByUnion(keys,tmpKeys);
      }
    }
  }
  if ( vtList.size() != 0 ) {
    std::set<unsigned> tmpKeys;
    if ( ! mvp.getMasterKeysWithVoteContent(tmpKeys,vtList,matchAll) ) return 4;
    if ( ! started ) {
      keys = tmpKeys;
      started = true;
    } else {
      if ( matchAll ) {
	mvp.mergeSetsByISect(keys,tmpKeys);
      } else {
	mvp.mergeSetsByUnion(keys,tmpKeys);
      }
    }
  }

  if ( akList.size() != 0 ) {
    std::set<unsigned> tmpKeys;
    if ( ! mvp.getMasterKeysWithAncilaryFiles(tmpKeys,akList,matchAll) ) return 4;
    if ( ! started ) {
      keys = tmpKeys;
      started = true;
    } else {
      if ( matchAll ) {
	mvp.mergeSetsByISect(keys,tmpKeys);
      } else {
	mvp.mergeSetsByUnion(keys,tmpKeys);
      }
    }
  }

  if ( ! mvp.printConfigs(std::cout,keys,brief) ) return 4;

  return 0;
}







