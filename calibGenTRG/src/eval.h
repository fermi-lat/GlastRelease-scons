// Read in python container ASCII files and return STL containers
// 
// Martin Kocian, SLAC, GLAST
// 
// 8/31/06
//
// 
#ifndef EVAL_H 
#define EVAL_H

#include <map>
#include <vector>
#include <string>
using std::map;
using std::vector;
using std::string;

class eval{
 public:
  eval(){}
  map<string,string> readdict(const string);
  const string unquote(const string sl);
  const vector<string> readlistortuple(const string l);
  
  int find(const string,char,int);
};

#endif
