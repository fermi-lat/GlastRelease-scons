// $Header$ 
#ifndef CALIBUTIL_STRIPSRV_H
#define CALIBUTIL_STRIPSRV_H

#include <string>
#include <iostream>
#include <vector>
#include "calibUtil/GenericSrv.h"

namespace calibUtil {

  typedef struct stowerRC { unsigned int row; unsigned int col;} towerRC; 

  class StripSrv {

  public: 

    /// this function takes in a stripList in string format and 
    /// fills a vector with corresponding strip numbers
    strToNum(std::string strips, std::vector<unsigned int> &v);

    /// constructor. Initializes strip service by creating a DOM structure 
    /// out the XML file and filling in the internal data structures 
    StripSrv(std::string xmlFileName);
  
    /// return the status (Hot or Dead) of the strip
    std::string getBadType();

    /// List all towers with bad strips 
    const std::vector<towerRC> getBadTowers();

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
    
    
    /// methods giving access to generic data

    /// Get instrument name
    std::string getInst();
    
    /// Get timestamp
    std::string getTimestamp();

    /// Get calibration type
    std::string getCalType();
    
    /// Get format Version
    std::string getFmtVer();

  private:

    typedef struct sUniLayer {
      std::string  stripType;
      std::vector<unsigned int> stripCol;
    }UniLayer;
  
    typedef struct sTray {
      unsigned int id;
      UniLayer *topLayer;
      UniLayer *botLayer;
    }Tray;

    typedef struct sTower {
      towerRC id;
      std::vector<Tray> trayCol;
    }Tower;
   
    std::vector<Tower> towerCol;
    std::string badType;

    /// Contained object to store generic data
    GenericSrv genSrv;   
    

  };

}/// end of namespace calibUtil

#endif











