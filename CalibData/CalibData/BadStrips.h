// $Header$
#ifndef CalibData_BadStrips_h
#define CalibData_BadStrips_h

#include "GaudiKernel/DataObject.h"

class calibUtil::StripSrv;

/** @class BadStrips
    @brief Base class for tracker hot, dead or "merged" (potentially
           including both hot and dead) strips.

           The 3 derived class (for now) add nothing new except
           implement getType method.
*/
namespace CalibData {
  class HeaderInfo;
  class ClientObject;

  class BadStrips : virtual public DataObject {
  public:
    enum eBadType {
      DEAD = 0,
      HOT  = 1,
      MERGED = 2 };

    /// Used by client callback in TraversInfo
    enum eUserRet {CONT, DONE, ERROR}
    
    enum eUnilayer {UNKNOWN_UNI, TOP, BOT};
      
    typedef struct sTowerRC {
      unsigned short int row;
      unsigned short int col;
    }      TowerRC;

    typedef std::vector<unsigned short int> StripCol;

    static const unsigned short int bDead     = 0x0800;
    static const unsigned short int bVeryDead = 0x1000;
    static const unsigned short int bHot      = 0x2000;
    static const unsigned short int bVeryHot  = 0x4000;
      
    /// Construct from calibUtil::StripSrv object
    BadStrips(calibUtil::StripSrv* utilStrips);
    ~BadStrips();

    BadType getType() = 0;
    
    /// Return count of bad (includes very bad) strips in specified tower
    int getNBad(const TowerRC& tower) const;

    /// Return count of bad (includes very bad) strips in (tower,tray)
    int getNBad(const TowerRC& tower, int tray) const;

    /// Return count of bad (includes very bad) strips in (tower,tray,unilayer)
    int getNBad(const TowerRC& tower, int tray, eUnilayer u) const;

    /// Return count of very bad strips in specified tower
    int getNVeryBad(const TowerRC& tower) const;

    /// Return count of very bad strips in (tower,tray)
    int getNVeryBad(const TowerRC& tower, int tray) const;

    /// Return count of very bad strips in (tower,tray,unilayer)
    int getNVeryBad(const TowerRC& tower, int tray, eUnilayer u) const;

    /// Return object holding universal calibration information
    HeaderInfo* getHeader() const;

    /// call back method for client to access all data
    bool traverseInfo(ClientObject *client) const;
    


  };     // end Badstrips class definition
}        // end CalibData namespace
#endif
