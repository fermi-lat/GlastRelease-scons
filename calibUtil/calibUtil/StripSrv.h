// $Header$ 
#ifndef CALIBUTIL_STRIPSRV_H
#define CALIBUTIL_STRIPSRV_H

#include <string>
#include <iostream>
#include <vector>

namespace calibUtil {

  typedef struct stowerRC { unsigned int row; unsigned int col;} towerRC;
  enum uniL {TOP, BOT};

  class ClientObject;
  class GenericSrv;

  class StripSrv {

  private: 

    /// this function takes in a stripList in string format and 
    /// fills a vector with corresponding strip numbers
    void strToNum(std::string strips, std::vector<unsigned int> &v);

  public:

    /// constructor. Initializes strip service by creating a DOM structure 
    /// out of the XML file and filling in the internal data structures 
    StripSrv(std::string xmlFileName);

    /// destructor. Deallocates memory
    ~StripSrv();
  
    /// returns the status (Hot or Dead) of the strip
    std::string getBadType();

    /// lists all towers with bad strips 
    const std::vector<towerRC> getBadTowers();

    /// counts very bad strips for the tower specified 
    unsigned int countVeryBad(towerRC towerId);
    
    /// counts bad strips (including very bad) for the tower specified 
    unsigned int countBad(towerRC towerId);
    
    /// counts very bad strips for the tower and tray specified 
    unsigned int countVeryBad(towerRC towerId, unsigned int trayNum);
    
    /// counts bad strips (including very bad) for the tower and trayspecified 
    unsigned int countBad(towerRC towerId, unsigned int trayNum);
    
    /// counts very bad strips for the tower,tray and unilayer  specified 
    unsigned int countVeryBad(towerRC towerId, unsigned int trayNum, 
                              uniL uniLayer);
    
    /// counts bad strips (including very bad) for the tower, tray 
    /// and unilayer specified 
    unsigned int countBad(towerRC towerId, unsigned int trayNum, 
                          uniL uniLayer);

    /// lists all very bad strips with the tower,tray and unilayer  
    std::vector<unsigned int> getVeryBad(towerRC towerId, 
                                         unsigned int trayNum, uniL uniLayer);
    
    /// lists  bad strips (including very bad) with the tower, tray 
    /// and unilayer specified 
    std::vector<unsigned int> getBad(towerRC towerId, 
                                     unsigned int trayNum, uniL uniLayer);
    
    
    /// methods giving access to generic data

    /// Get instrument name
    std::string getInst();
    
    /// Get timestamp
    std::string getTimestamp();

    /// Get calibration type
    std::string getCalType();
    
    /// Get format Version
    std::string getFmtVer();

    /// call back method for client to access large data
    void StripSrv::traverseInfo(ClientObject *client);

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

    // object to store generic data
    GenericSrv *genSrv;   
    

  };

}// end of namespace calibUtil

#endif











