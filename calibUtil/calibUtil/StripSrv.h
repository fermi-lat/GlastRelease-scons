// $Header$ 
#ifndef CALIBUTIL_STRIPSRV_H
#define CALIBUTIL_STRIPSRV_H

#include <string>
#include <iostream>
#include <vector>

namespace calibUtil {


  class ClientObject;
  class GenericSrv;

  class StripSrv {


  public:

    enum eUnilayer {UNKNOWN_UNI, TOP, BOT};
    enum eBadType  {UNKNOWN_BADTYPE, DEAD, HOT};

    typedef struct stowerRC { unsigned int row; unsigned int col;} towerRC;

    /// constructor. Initializes strip service by creating a DOM structure 
    /// out of the XML file and filling in the internal data structures 
    StripSrv(std::string xmlFileName);

    /// destructor. Deallocates memory
    ~StripSrv();
  
    /// returns the status (Hot or Dead) of the strip
    //    std::string getBadType();
    eBadType getBadType() const;

    /// lists all towers with bad strips 
    // const std::vector<towerRC> 
    void getBadTowers(std::vector<towerRC>& towerIds) const;

    /// counts very bad strips for the tower specified 
    unsigned int nVeryBad(towerRC towerId) const;
    
    /// counts bad strips (including very bad) for the tower specified 
    unsigned int nBad(towerRC towerId) const;
    
    /// counts very bad strips for the tower and tray specified 
    unsigned int nVeryBad(towerRC towerId, unsigned int trayNum) const;
    
    /// counts bad strips (including very bad) for the tower and trayspecified 
    unsigned int nBad(towerRC towerId, unsigned int trayNum) const;
    
    /// counts very bad strips for the tower,tray and unilayer  specified 
    unsigned int nVeryBad(towerRC towerId, unsigned int trayNum, 
                          eUnilayer uni) const;
    
    /// counts bad strips (including very bad) for the tower, tray 
    /// and unilayer specified 
    unsigned int nBad(towerRC towerId, unsigned int trayNum, 
                      eUnilayer uni) const;

    /// lists all very bad strips for a particular tower,tray and unilayer  
    //    std::vector<unsigned int> getVeryBad(towerRC towerId, unsigned int trayNum,
    //                                   eUnilayer uni);
    void getVeryBad(towerRC towerId, unsigned int trayNum, eUnilayer uni,
                    std::vector<unsigned int>& strips) const;

    /// lists  bad strips (including very bad) with the tower, tray 
    /// and unilayer specified 
    void getBad(towerRC towerId, unsigned int trayNum, eUnilayer uni,
                std::vector<unsigned int>& strips) const; 

    // std::vector<unsigned int> getBad(towerRC towerId, unsigned int trayNum,
    //                                 eUnilayer uniLayer);
    
    
    /// methods giving access to generic data

    /// Get instrument name
    std::string getInst() const;
    
    /// Get timestamp
    std::string getTimestamp() const;

    /// Get calibration type
    std::string getCalType() const;
    
    /// Get format Version
    std::string getFmtVer() const;

    /// call back method for client to access large data
    void StripSrv::traverseInfo(ClientObject *client) const;

  private:

    typedef struct sUnilayer {
      std::string  stripType;
      std::vector<unsigned int> stripCol;
    } Unilayer;
  
    typedef struct sTray {
      unsigned int id;
      Unilayer *topLayer;
      Unilayer *botLayer;
    } Tray;

    typedef struct sTower {
      towerRC id;
      std::vector<Tray> trayCol;
    } Tower;
   
    std::vector<Tower> m_towers;
    eBadType m_badType;

    // object to store generic data
    GenericSrv *m_genSrv;   
    
    /// this function takes in a stripList in string format and 
    /// fills a vector with corresponding strip numbers
    void strToNum(std::string strips, std::vector<unsigned int> &v);

  };

}// end of namespace calibUtil

#endif











