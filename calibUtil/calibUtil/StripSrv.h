// $Header$ 
#ifndef CALIBUTIL_STRIPSRV_H
#define CALIBUTIL_STRIPSRV_H

#include <string>
#include <iostream>
#include <vector>

namespace calibUtil {

  typedef unsigned int towerRC [2]; 

  class StripSrv {
  public: 

    /// constructor. Initializes strip service by creating a DOM structure 
    /// out the XML file and filling in the internal data structures 
    StripSrv(std::string xmlFileName);
  
    /// return the status (Hot or Dead) of the strip
    std::string getBadType();

    /// List all towers with bad strips 
    std::vector<towerRC> getBadTowers();

    /// count for very bad strips for the tower specified 
    unsigned int countVeryBad(towerRC towerId);
    
    /// count for  bad strips (including very bad) for the tower specified 
    unsigned int countBad(towerRC towerId);
    
    /// count for very bad strips for the tower and tray specified 
    unsigned int countVeryBad(towerRC towerId, unsigned int trayNum);
    
    /// count for  bad strips (including very bad) for the tower and tray 
    /// specified 
    unsigned int countBad(towerRC towerId, unsigned int trayNum);
    
    /// count for very bad strips for the tower,tray and unilayer  specified 
    unsigned int countVeryBad(towerRC towerId, unsigned int trayNum, std::string uniLayer);
    
    /// count for  bad strips (including very bad) for the tower, tray 
    /// and unilayer specified 
    unsigned int countBad(towerRC towerId, unsigned int trayNum, std::string uniLayer);

    /// List all very bad strips with the tower,tray and unilayer  
    std::vector<unsigned int> getVeryBad(towerRC towerId, unsigned int trayNum, std::string uniLayer);
    
    /// List  bad strips (including very bad) with the tower, tray 
    /// and unilayer specified 
    std::vector<unsigned int> getBad(towerRC towerId, unsigned int trayNum, std::string uniLayer);
    
    
  private:

    typedef struct UniLayer {
      std::string  stripType;
      std::vector<unsigned int> stripList;
    };
  
    typedef struct Tray {
      unsigned int trayId;
      UniLayer* topLayer;
      UniLayer* botLayer;
    };
  
    typedef struct Tower {
      towerRC towerId;
      std::vector<Tray> trayList;
    };
   
    std::vector<Tower> towerList;
  };

}/// end of namespace calibUtil

#endif











