// $Header$ 
#ifndef CALIBUTIL_STRIPSRV_H
#define CALIBUTIL_STRIPSRV_H

#include <string>
#include <iostream>
#include <vector>
#include "dom/DOM_Element.hpp"

namespace calibUtil {


  class ClientObject;
  class GenericSrv;

  class StripSrv {


  public:

    enum eUnilayer {UNKNOWN_UNI, TOP, BOT};
    enum eBadType  {UNKNOWN_BADTYPE, DEAD, HOT};
    enum eBadness {VERYBAD=0, BAD, NBADNESS};
    enum eRet     {CONT, USER_DONE, ERROR, DONE};

    /// Clients shoud use as return values for readData


    typedef struct stowerRC { 
      unsigned short int row; 
      unsigned short int col;
    } towerRC;
    typedef std::vector<unsigned short int> StripCol; 

    /// constructor. Initializes strip service by creating a DOM structure 
    /// out of the XML file and filling in the internal data structures 
    StripSrv(std::string xmlFileName);

    /// Constructor to be used when creating new calibrations
    StripSrv(eBadType badType, const GenericSrv& gen);

    /** Add additional strips to the collection.  They must all come from
        the same (uni)layer of Si and must all be of the same "badness"
        @return  number of strips added
    */
    bool addBad(eBadness howBad, towerRC& t, int tray, eUnilayer uni,
                const StripCol& list);

    /// destructor. Deallocates memory
    ~StripSrv();
  
    /// returns the status (Hot or Dead) of the strip
    //    std::string getBadType();
    eBadType getBadType() const;

    /// lists all towers with bad strips 
    void getBadTowers(std::vector<towerRC>& towerIds) const;

    /// counts very bad strips for the tower specified 
    unsigned int nVeryBad(const towerRC& towerId) const;
    
    /// counts bad strips (including very bad) for the tower specified 
    unsigned int nBad(const towerRC& towerId) const;
    
    /// counts very bad strips for the tower and tray specified 
    unsigned int nVeryBad(const towerRC& towerId, unsigned int trayNum) const;
    
    /// counts bad strips (including very bad) for the tower and trayspecified 
    unsigned int nBad(const towerRC& towerId, unsigned int trayNum) const;
    
    /// counts very bad strips for the tower,tray and unilayer  specified 
    unsigned int nVeryBad(const towerRC& towerId, unsigned int trayNum, 
                          eUnilayer uni) const;
    
    /// counts bad strips (including very bad) for the tower, tray 
    /// and unilayer specified 
    unsigned int nBad(const towerRC& towerId, unsigned int trayNum, 
                      eUnilayer uni) const;

    /// lists all very bad strips for a particular tower,tray and unilayer  
    //    std::vector<unsigned int> getVeryBad(towerRC towerId, unsigned int trayNum,
    //                                   eUnilayer uni);
    void getVeryBad(const towerRC& towerId, unsigned int trayNum, eUnilayer uni,
                    StripCol& strips) const;

    /// lists  bad strips (including very bad) with the tower, tray 
    /// and unilayer specified 
    void getBad(const towerRC& towerId, unsigned int trayNum, eUnilayer uni,
                StripCol& strips) const; 

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
    StripSrv::eRet StripSrv::traverseInfo(ClientObject *client) const;

    eRet writeXml(std::ostream* out);

  private:
    // A unilayer can have one or both badnesses ("bad" and "veryBad")
    typedef struct sUnilayer {
      StripCol badLists[NBADNESS];
    } Unilayer;
  
    typedef struct sTray {
      unsigned int id;
      Unilayer *top;  //    Unilayer *topLayer;
      Unilayer *bot;  // Unilayer *botLayer;
    } Tray;

    typedef struct sTower {
      towerRC id;
      std::vector<Tray> trays;
    } Tower;
   
    std::vector<Tower> m_towers;
    eBadType m_badType;

    enum eState {
      FROM_PERS,       // read in from persistent form; read-only
      BUILDING,        // in process of building the list
      WRITTEN };       // built and written to persistent form (read-only)

    eState m_state;

    // object to store generic data
    GenericSrv *m_genSrv;   

    // A some point need to define a "cuts" object or just add three
    // new fields
    // float  expected;  // expected percentage of hits for a "normal"
                         //  channel for this input sample
    // float  tight;     // max (min) % of hits for a channel to be
                         // considered not hot (not dead)
    // float  loose;     // max (min) % of hits for a channel to be
                         // considered not very hot (not very dead)

    
    /// this function takes in a stripList in string format and 
    /// fills a vector with corresponding strip numbers
    void strToNum(std::string strips, std::vector<unsigned short int> &v);

    ///  Handles all the messy of details of extracting information
    ///  about a single unilayer from the XML representation
    void fillUni(const DOM_Element& uniElt, Unilayer* uni);

    void fillStrips(const DOM_Element& badElt, StripCol& list);
    
    /// Internal utility, used when this object must be mutable (e.g. addBad)
    Tower* findTower(towerRC& towerId);
    /// const version of findTower, used by const public methods
    const Tower* findTower(const towerRC& towerId) const;

    /// Internal utility, used when this object must be mutable 
    Tray* findTray(Tower& tower, unsigned int trayNum);

    /// const version of findTray, used by const public methods
    const Tray* findTray(const Tower& tower, unsigned int trayNum) const;

    bool addBad(eBadness howBad,  Tower* tow, int tray,
                eUnilayer uni, const StripCol& list);


    //    ----> TO BE WRITTEN           <---
    bool addBad(eBadness howBad,  Tray* tr, eUnilayer uni, 
                const StripCol& list);

  };

}// end of namespace calibUtil

#endif











