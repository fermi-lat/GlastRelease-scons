// $Header$
/**  @file test_strips.cxx
   Sample program to exercise low-level calibration strip services.  
*/
#include "MyObject.h"

#include <string>
#include <iostream>
#include <fstream>


int main(int argc, char* argv[]) {
  using calibUtil::StripSrv;

  std::string name(argv[1]);
  
  StripSrv ssObj(name);

  std::vector<StripSrv::towerRC> tv1;
  ssObj.getBadTowers(tv1);   
  cout << "ROW IS" << tv1[0].row;
  cout << "COL IS" << tv1[0].col;

  StripSrv::towerRC trc;
  trc.row = 3;
  trc.col = 2; 

  cout << "bad is" <<  ssObj.nBad(trc,2,StripSrv::TOP) << endl;

  std::vector<unsigned int> stripList;
  ssObj.getBad(trc,0,StripSrv::TOP,stripList);
  cout << "strip number is" << stripList[0];

  cout<< "BAD TYPE IS" << ssObj.getBadType() << std::endl;

  cout << "instrument name is" << ssObj.getCalType() << endl;

  MyObject cli;
  ssObj.traverseInfo(&cli);

  return(0);

} /* end of main */









