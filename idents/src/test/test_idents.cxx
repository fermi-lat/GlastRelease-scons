// $Header$


#include "idents/VolumeIdentifier.h"
#include "idents/AcdId.h"
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

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
  idents::AcdId acdId(acdVolId);

  std::cout << "AcdVolId = " << acdVolId.name() << std::endl;
  std::cout << "AcdId = " << acdId.id() << std::endl;
  std::cout << acdId.tile() << " " << acdId.ribbon() << " " << acdId.ribbonNum() 
      << " " << acdId.ribbonOrientation() << std::endl;
  
  idents::VolumeIdentifier acdVolId2;
  acdVolId2.append(1); acdVolId2.append(0); acdVolId2.append(41);
  acdVolId2.append(0); acdVolId2.append(3);
  idents::AcdId acdId2(acdVolId2);

  std::cout << "AcdVolId2 = " << acdVolId2.name() << std::endl;
  std::cout << "AcdId2 = " << acdId2.id() << std::endl;
  std::cout << acdId2.tile() << " " << acdId2.ribbon() << " " << acdId2.ribbonNum() 
      << " " << acdId2.ribbonOrientation() << std::endl;

  idents::VolumeIdentifier acdVolId3;
  acdVolId3.append(1); acdVolId3.append(1); acdVolId3.append(40);
  acdVolId3.append(0); acdVolId3.append(3);
  idents::AcdId acdId3(acdVolId3);

  std::cout << "AcdVolId3 = " << acdVolId3.name() << std::endl;
  std::cout << "AcdId3 = " << acdId3.id() << std::endl;
  std::cout << acdId3.tile() << " " << acdId3.ribbon() << " " << acdId3.face() 
      << " " << acdId3.row() << " " << acdId3.column() << std::endl;

  return 0;
}
