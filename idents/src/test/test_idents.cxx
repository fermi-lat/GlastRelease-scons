// $Header$


#include "idents/VolumeIdentifier.h"
#include "idents/AcdId.h"
#include "idents/CalXtalId.h"
#include "idents/TkrId.h"
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>

int main() 
{
  idents::VolumeIdentifier id1, id2, id3;
  id1.append(1); id1.append(1);
  id2.append(1); id2.append(1); id2.append(0);id2.append(0);
  id3.append(1);id3.append(1);id3.append(0);

  std::vector< idents::VolumeIdentifier > idVect;
  idVect.push_back(id1);
  idVect.push_back(id3);
  idVect.push_back(id2);

  std::sort(idVect.begin(), idVect.end());
  
  std::vector< idents::VolumeIdentifier >::iterator i;
  std::cout << "Sorted vector of volume identifiers: "<< std::endl;
  for (i=idVect.begin(); i != idVect.end(); i++)
    {
    std::cout << i->name() << std::endl;
    }
  std::cout << std::endl;

  std::map<idents::VolumeIdentifier,double> idMap;
  idMap[id1] = 1.5;
  idMap[id2] = 2.5;
  idMap[id3] = 3.5;

  std::map<idents::VolumeIdentifier,double>::iterator m;
  std::cout << "Map elements: " << std::endl;
  for (m = idMap.begin(); m != idMap.end();m++)
    {
      std::cout << "Id = " << (*m).first.name();
      std::cout << "   Value = " << (*m).second << std::endl;
    }

  // Test the AcdIds
  idents::VolumeIdentifier acdVolId;
  acdVolId.append(1); acdVolId.append(0); acdVolId.append(41);
  acdVolId.append(1); acdVolId.append(2);
  try {
    idents::AcdId acdId(acdVolId);
    std::cout << "AcdVolId = " << acdVolId.name() << std::endl;
    std::cout << "AcdId = " << acdId.id() << std::endl;
    std::cout << acdId.tile() << " " << acdId.ribbon() << " " << acdId.ribbonNum() 
      << " " << acdId.ribbonOrientation() << std::endl;
  
    std::cout << "Value of acdVolId.isAcd() is " << acdVolId.isAcd() 
            << std::endl;
  }
  catch (std::invalid_argument ex) {
    std::cout << 
      "Caught expected invalid_argument exception in AcdId constructor " 
              << "with input acdVolId" << std::endl;
  }


  try {
    idents::VolumeIdentifier acdVolId2;
    acdVolId2.append(1); acdVolId2.append(0); acdVolId2.append(41);
    acdVolId2.append(0); acdVolId2.append(3);
    idents::AcdId acdId2(acdVolId2);

    std::cout << "AcdVolId2 = " << acdVolId2.name() << std::endl;
    std::cout << "AcdId2 = " << acdId2.id() << std::endl;
    std::cout << acdId2.tile() << " " << acdId2.ribbon() << " " << acdId2.ribbonNum() 
      << " " << acdId2.ribbonOrientation() << std::endl;

    std::cout << "Value of acdVolId2.isAcd() is " << acdVolId.isAcd() 
            << std::endl;
    } catch(std::invalid_argument ex){
        std::cout << 
        "Caught expected invalid_argument exception in AcdId constructor " 
              << "with input acdVolId2" << std::endl;
  }

  try {
      idents::VolumeIdentifier acdVolId3;
      acdVolId3.append(1); acdVolId3.append(1); acdVolId3.append(40);
      acdVolId3.append(0); acdVolId3.append(3);
      idents::AcdId acdId3(acdVolId3);
    
      std::cout << "AcdVolId3 = " << acdVolId3.name() << std::endl;
      std::cout << "AcdId3 = " << acdId3.id() << std::endl;
      std::cout << acdId3.tile() << " " << acdId3.ribbon() << " " << acdId3.face() 
          << " " << acdId3.row() << " " << acdId3.column() << std::endl;

      std::cout << "Value of acdVolId3.isAcd() is " << acdVolId.isAcd() 
            << std::endl;

      std::cout << "Value of acdVolId3.isCal() is " << acdVolId.isCal() 
            << std::endl;

      std::cout << "Value of acdVolId3.isTkr() is " << acdVolId.isTkr() 
            << std::endl;
  } catch (std::invalid_argument ex) {
    std::cout << 
      "Caught expected invalid_argument exception in AcdId constructor " 
              << "with input acdVolId3" << std::endl;
  }
  

  idents::VolumeIdentifier vIdCal, vIdBad;
  unsigned yNum = 1, xNum = 2, layer = 5, column = 10;
  vIdCal.append(0);        // LATobject = towers
  vIdCal.append(yNum); 
  vIdCal.append(xNum);        // xnum
  vIdCal.append(0);        // towerObject = Calorimeter
  vIdCal.append(layer);        // layer
  vIdCal.append(1);        // orientation. In real LAT layer 5 is y-measuring
  vIdCal.append(column);       // column

  std::cout << "Value of vIdCal.isAcd() is " << vIdCal.isAcd() 
            << std::endl;
  std::cout << "Value of vIdCal.isCal() is " << vIdCal.isCal() 
            << std::endl;
  std::cout << "Value of vIdCal.isTkr() is " << vIdCal.isTkr() 
            << std::endl;

  vIdBad.append(0);        // LATobject = towers
  vIdBad.append(yNum); 
  vIdBad.append(xNum);        // xnum
  vIdBad.append(1);        // towerObject = tracker
  vIdBad.append(layer);        // layer
  vIdBad.append(0);        // orientation
  vIdBad.append(column);       // column

  std::cout << "Value of vIdBad.isAcd() is " << vIdBad.isAcd() 
            << std::endl;
  std::cout << "Value of vIdBad.isCal() is " << vIdBad.isCal() 
            << std::endl;
  std::cout << "Value of vIdBad.isTkr() is " << vIdBad.isTkr() 
            << std::endl;


  idents::CalXtalId xtalOk(vIdCal);
  if ((xtalOk.getTower() != (int) (4*yNum  + xNum) ) ||
      (xtalOk.getLayer() != (int) layer)             ||
      (xtalOk.getColumn() != (int) column) ) {
    std::string badConstructor("CalXtalId constructor from vIdCal failed");
    throw std::logic_error(badConstructor);
  }
  else {
    std::cout << "Built CalXtalId properly from vIdCal" << std::endl;
    idents::VolumeIdentifier* loopVid = xtalOk.makeVolumeId();
    std::cout << "Original vIdCal has name " << vIdCal.name() << std::endl;
    std::cout << "loop vid has name " << loopVid->name() << std::endl;
  }
  
  try {
    idents::CalXtalId xtalBad(vIdBad);
    std::string 
      falsePositive("CalXtalId constructor for vIdBad failed to fail!");
    throw std::logic_error(falsePositive);
  }
  catch (std::invalid_argument ex) {
    std::cout << 
      "Caught expected invalid_argument exception in CalXtalId constructor " 
              << "with input vIdBad" << std::endl;
  }

  // Tests for TkrId class. Shouldn't be able to use valid acd id as input:
  try {
    idents::TkrId xtalBad(acdVolId);
    std::string 
      falsePositive("TkrId constructor for acdVolId failed to fail!");
    throw std::logic_error(falsePositive);
  }
  catch (std::invalid_argument ex) {
    std::cout << 
      "Caught expected invalid_argument exception in TkrId constructor " 
              << "with input acdVolId " << std::endl;
  }

  // Make a sensible tkr VolumeIdentifier
  idents::VolumeIdentifier vIdTkr, vIdTkrTrunc;
  unsigned tray = 3;
  unsigned view = 1;
  unsigned botTop = 0;  // bottom
  unsigned ladder = 1;
  unsigned wafer = 2;

  bool isTop = (botTop == 1);
  idents::TkrId noView(xNum, yNum, tray, isTop);
  idents::TkrId withViewX(xNum, yNum, tray, isTop, idents::TkrId::eMeasureX);
  idents::TkrId withViewY(xNum, yNum, tray, isTop, idents::TkrId::eMeasureY);
  idents::TkrId withViewBad(xNum, yNum, tray, isTop, 17);

  std::cout << "For noView validity of view is " << noView.hasView() 
            << std::endl;
  std::cout << "For withViewX validity of view is " << withViewX.hasView() 
            << " and view is " << withViewX.getView() << std::endl;
  std::cout << "For withViewY validity of view is " << withViewY.hasView() 
            << " and view is " << withViewY.getView() << std::endl;

  std::cout << "For withViewBad validity of view is " << withViewBad.hasView() 
            << std::endl;




  vIdTkr.append(0); vIdTkrTrunc.append(0);        // LATobject = towers
  vIdTkr.append(yNum); vIdTkrTrunc.append(yNum);
  vIdTkr.append(xNum); vIdTkrTrunc.append(xNum);
  vIdTkr.append(1);    vIdTkrTrunc.append(1);    // towerObject = Tracker
  vIdTkr.append(tray); vIdTkrTrunc.append(tray);
  vIdTkr.append(view); vIdTkrTrunc.append(view);


  vIdTkr.append(botTop);
  vIdTkr.append(ladder);
  vIdTkr.append(wafer);

  try {
    idents::TkrId tId(vIdTkr);

    std::cout << "Made TkrId with TowerY, TowerX = " 
              << yNum << ", " << xNum << std::endl;
    std::cout << " tray = " << tray << ", view = " << view << std::endl;
    std::cout << " (botTop, ladder, wafer) = (" << botTop << ", ";
    std::cout << ladder << ", " << wafer << ")" << std::endl;

    std::cout << "Getting tower Y and x : "  
              <<  tId.getTowerY() << ", " << tId.getTowerX() << std::endl;
    std::cout << "Getting tray and view : "
              << tId.getTray() << ", " << tId.getView() << std::endl;
    std::cout << "Getting botTop, ladder, wafer : "
              << tId.getBotTop() << ", " << tId.getLadder() 
              <<   ", " << tId.getWafer() << std::endl;
  }
  catch (std::invalid_argument ex) {
    std::cout << "Failed to construct TkrId from vIdTkr " << std::endl;
  }

  try {
    idents::TkrId tIdTrunc(vIdTkrTrunc);
    std::cout << "tIdTrunc is missing botTop,ladder,wafer " << std::endl;
    std::cout << "vIdTkrTrunc size is " << vIdTkrTrunc.size() << std::endl;
    unsigned val;
    std::cout << "Fetching fields from corresponding TkrId.. " << std::endl;
    val = tIdTrunc.getTowerY(); std::cout << val << " "; std::cout.flush();
    val = tIdTrunc.getTowerX(); std::cout << val << " "; std::cout.flush();
    val = tIdTrunc.getTray(); std::cout << val << " "; std::cout.flush();
    val = tIdTrunc.getView(); std::cout << val << " "; std::cout.flush();
    val = tIdTrunc.getBotTop(); std::cout << val << " "; std::cout.flush();
    val = tIdTrunc.getLadder(); std::cout << val << " "; std::cout.flush();
    val = tIdTrunc.getWafer(); std::cout << val << " " << std::endl; 
    std::cout.flush();
  }
  catch (std::domain_error ex) {
    std::string what(ex.what());
    std::cout << std::endl;
    std::cout << "Examining tIdTrunc fields, caught exception " 
              << what << std::endl;
  }
  catch (std::invalid_argument ex) {
    std::cout << "Failed to construct TkrId from vIdTrunc " << std::endl;
  }

  return 0;
}
