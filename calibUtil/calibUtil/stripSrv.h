#include "xml/XmlParser.h"
#include "xml/Dom.h"
#include "dom/DOM_Element.hpp"
#include "dom/DOM_NodeList.hpp"
#include "dom/DOM_TreeWalker.hpp"

#include <string>
#include <iostream>
#include <fstream>

namespace stripSrv {

  class stripSrv {
  public: 

    typedef unsigned int towerRC [2];

    // return the status (Hot or Dead) of the strip
    std::string getStripsStatus();

    // List all towers with bad strips 
    vector<towerRC> getBadTowers();

    // count for very bad strips for the tower specified 
    unsigned int countVeryBadStrips(towerRC towerId);
    
    // count for  bad strips (including very bad) for the tower specified 
    unsigned int countBadStrips(towerRC towerId);
    
    // count for very bad strips for the tower and tray specified 
    unsigned int countVeryBadStrips(towerRC towerId, unsigned int trayNum);
    
    // count for  bad strips (including very bad) for the tower and tray 
    // specified 
    unsigned int countBadStrips(towerRC towerId, unsigned int trayNum);
    
    // count for very bad strips for the tower,tray and unilayer  specified 
    unsigned int countVeryBadStrips(towerRC towerId, unsigned int trayNum, std::string uniLayer);
    
    // count for  bad strips (including very bad) for the tower, tray 
    // and unilayer specified 
    unsigned int countBadStrips(towerRC towerId, unsigned int trayNum, std::string uniLayer);

    // List all very bad strips with the tower,tray and unilayer  
    vector<unsigned int> getVeryBadStrips(towerRC towerId, unsigned int trayNum, std::string uniLayer);
    
    // List  bad strips (including very bad) with the tower, tray 
    // and unilayer specified 
    vector<unsigned int> getBadStrips(towerRC towerId, unsigned int trayNum, std::string uniLayer);
    
  };
}// end of namespace stripSrv













