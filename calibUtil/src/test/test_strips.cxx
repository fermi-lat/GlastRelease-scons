// $Header$
/**  @file test_strips.cxx
   Sample program to exercise low-level calibration strip services.  
*/
#include "MyObject.h"

#include <string>
#include <iostream>
#include <fstream>


int main(int, char* argv[]) {
  using calibUtil::StripSrv;

  std::string name(argv[1]);
  
  StripSrv ssObj(name);

  std::vector<StripSrv::towerRC> trc;
  ssObj.getBadTowers(trc);   
  for (unsigned int i = 0; i < trc.size(); i++ ) {
    std::cout << "Tower id is (" << trc[i].row << ", " << trc[i].col;
    std::cout << ")" << std::endl;
    std::cout << "# very bad is " << ssObj.nVeryBad(trc[i]) << std::endl;
    std::cout << "# bad is " << ssObj.nBad(trc[i]) << std::endl;
  }
  std::cout << std::endl;

  std::cout << "For first tower, # very bad for each tray: " << std::endl;
  for (unsigned int iTray = 0; iTray <= 18; iTray++)  {
    std::cout << "  Tray " << iTray <<  " #very bad   " 
              << ssObj.nVeryBad(trc[0], iTray)  << " #bad, TOP  " 
              << ssObj.nBad(trc[0], iTray, StripSrv::TOP) << std::endl;
  }

  
  StripSrv::StripCol  strips;
  unsigned int        trayNum = 0;
  
  ssObj.getBad(trc[0], trayNum ,StripSrv::TOP, strips);


  std::cout<< "BAD TYPE IS" << (int) ssObj.getBadType() << std::endl;

  std::cout << "instrument name is" << ssObj.getCalType() << std::endl;

  MyObject cli;
  ssObj.traverseInfo(&cli);

  return(0);

} /* end of main */


calibUtil::StripSrv::eRet
 MyObject::readData(calibUtil::StripSrv::towerRC towerId, 
                   unsigned int trayNum, 
                   calibUtil::StripSrv::eUnilayer uni, 
                   calibUtil::StripSrv::eBadness howBad,
                   const calibUtil::StripSrv::StripCol* const strips) {
  using calibUtil::StripSrv;

  std::cout << "In readData for tower (" << towerId.row << ", " << towerId.col;
  std::cout << ") tray #" << trayNum << "  unilayer ";
  if (uni == StripSrv::TOP) std::cout << "TOP ";
  if (uni == StripSrv::BOT) std::cout << "BOT ";
  std::cout << "  Badness: ";
  if (howBad == StripSrv::VERYBAD) std::cout << " VERYBAD";
  if (howBad == StripSrv::BAD) std::cout << " JUST BAD";

  std::cout << std::endl;
  std::cout << "Bad strip count:  " << strips->size() << std::endl;
  std::cout << std::endl;
  
  return StripSrv::CONT;
}







