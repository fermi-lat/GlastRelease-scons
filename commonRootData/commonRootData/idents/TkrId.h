#ifndef commonRootData_TKRID_H
#define commonRootData_TKRID_H 1

#include <stdexcept>
#include "TObject.h"

namespace commonRootData {
    
/**
*
* @class   TkrId        
*  
* @brief ID class for Tkr Id, based on TkrId class in idents, excluding
*  all mention of VolumeIdentifier
*  Includes access methods for  all tracker-relevant fields.  
* @author  J. Bogart
*
* $Header$
*/

  class TkrId : public TObject {
  public:

    /** Identify top or bottom Silicon layer.
        Should have same values as identically-named constants in xml 
        description.
    */
    enum
      {
        eTKRSiBottom = 0,
        eTKRSiTop = 1
      };
        
    /** Identify x-measuring or y-measuring.
        Should have same values as identically-named constants in xml 
        description.
    */
    enum
      {
        eMeasureX = 0,
        eMeasureY = 1,
        eMeasureNone = 2,
        eNoValue = 100      // used for missing ladder, wafer
      };
    /** constructor including just enough information to identify 
        silicon unilayer (view, ladder, wafer all optional)
    */
    TkrId(unsigned towerX, unsigned towerY, unsigned tray, bool top, 
          int view=eMeasureNone, int ladder=eNoValue, int wafer=eNoValue): 
      TObject() {
      init(towerX, towerY, tray, top, view, ladder, wafer);
    }
            
    /// Unuseful default constructor
    TkrId() : TObject(), m_packedId(0) {};

    /// Copy constructor
    TkrId(const TkrId& other) : TObject(other), m_packedId(other.m_packedId) {};

    ~TkrId() {};

    void Clear(Option_t * ="") {m_packedId = 0;}

    bool isEqual(const TkrId& other) {
      return ((m_packedId == other.m_packedId));
    }
    
    //    const bool operator<(const TkrId& right) const 
    //    {return m_packedId < right.m_packedId;}
    

    void init(unsigned towerX, unsigned towerY, unsigned tray, bool top, 
         int view=eMeasureNone, int ladder=eNoValue, int wafer=eNoValue);
        
    void init(UInt_t packedId) {
      m_packedId = packedId;
    }

  private:
    enum {
      MASKTowerY = 0x3,
      MASKTowerX = 0x3,
      MASKTray   = 0x3f,
      MASKMeas   = 0x1,
      MASKBotTop = 0x1,
      MASKLadder = 0x3,
      MASKWafer  = 0x3
    };

    enum {
      SHIFTTowerY = 0,
      SHIFTTowerX = 2,
      SHIFTTray   = 4,
      SHIFTMeas   = 10,
      SHIFTBotTop = 11,
      SHIFTLadder = 12,
      SHIFTWafer  = 14
    };

    /*
#define TKRID_VALIDTowerY (unsigned) 0x40000000
#define TKRID_VALIDTowerX (unsigned) 0x20000000
#define TKRID_VALIDTray  (unsigned) 0x10000000
#define TKRID_VALIDMeas   (unsigned) 0x8000000
#define TKRID_VALIDBotTop (unsigned) 0x4000000
#define TKRID_VALIDLadder (unsigned) 0x2000000
#define TKRID_VALIDWafer  (unsigned) 0x1000000
    */
    enum {
      VALIDTowerY =  0x40000000,
      VALIDTowerX =  0x20000000,
      VALIDTray =    0x10000000,
      VALIDMeas =   0x08000000,
      VALIDBotTop = 0x04000000,
      VALIDLadder =  0x02000000,
      VALIDWafer =   0x01000000
    };


    enum {
      SHMASKTowerY = MASKTowerY << SHIFTTowerY,
      SHMASKTowerX = MASKTowerX << SHIFTTowerX,
      SHMASKTray = MASKTray << SHIFTTray,
      SHMASKMeas = MASKMeas << SHIFTMeas,
      SHMASKBotTop = MASKBotTop << SHIFTBotTop,
      SHMASKLadder = MASKLadder << SHIFTLadder,
      SHMASKWafer = MASKWafer << SHIFTWafer
    };
  public:

                        
    bool hasTowerX() const {return ((m_packedId & VALIDTowerX) != 0);}
    unsigned int getTowerX() const {
      if (!(hasTowerX())) throw std::domain_error("No TowerX field");
      return (m_packedId & SHMASKTowerX) >> SHIFTTowerX;
    }

    bool hasTowerY() const {return ((m_packedId & VALIDTowerY) != 0);}
    unsigned int getTowerY() const  {
      if (!(hasTowerY())) throw std::domain_error("No TowerY field");
      return (m_packedId & SHMASKTowerY) >> SHIFTTowerY;
    }


    bool hasTray() const {return ((m_packedId & VALIDTray) != 0);}
    unsigned int getTray() const {
      if (!(hasTray())) throw std::domain_error("No Tray field");
      return (m_packedId & SHMASKTray) >> SHIFTTray;
    }


    bool hasBotTop() const {return ((m_packedId & VALIDBotTop) != 0);}
    unsigned int getBotTop() const {
      if (!(hasBotTop())) throw std::domain_error("No BotTop field");
      return (m_packedId & SHMASKBotTop) >> SHIFTBotTop;
    }


    bool hasView() const {return ((m_packedId & VALIDMeas) != 0);}
    unsigned int getView() const {
      if (!(hasView())) throw std::domain_error("No View field");
      return (m_packedId & SHMASKMeas) >> SHIFTMeas;
    }


    bool hasLadder() const {return ((m_packedId & VALIDLadder) != 0);}
    unsigned int getLadder() const {
      if (!(hasLadder())) throw std::domain_error("No Ladder field");
      return (m_packedId & SHMASKLadder) >> SHIFTLadder;
    }

    
    bool hasWafer() const {return ((m_packedId & VALIDWafer) != 0);}
    unsigned int getWafer() const {
      if (!(hasWafer())) throw std::domain_error("No Wafer field");
      return (m_packedId & SHMASKWafer) >> SHIFTWafer;
    }

    // Do we need this?
    //Access Methods for Tkr reconstruction semantics:
    unsigned int getLayer() const {return ( getPlane())/2;}
    unsigned int getPlane() const {return 2*getTray() + getBotTop() - 1;}

    /// Number of valid fields, never more than 7
    //    unsigned int getSize();

    /// conversion to int
    operator int() const {return m_packedId;}
            
  private:

    /// Bitmask containing tracker information
    //    unsigned short int m_packedId;
    UInt_t  m_packedId;

    ClassDef(TkrId,1) // Tkr Identifier
  };

} // namespace commonRootData
#endif    // commonRootData_TKRID_H
