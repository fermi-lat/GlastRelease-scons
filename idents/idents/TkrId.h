#ifndef idents_TKRID_H
#define idents_TKRID_H 1

#include <stdexcept>

namespace idents {
    
/**
*
* @class   TkrId        
*  
* @brief ID class for Tkr Id. Construct from a VolumeIdentifier.
*  Includes access methods for  all tracker-relevant fields.  
*  Will break if manner of assembling geometry in xmlGeoDbs changes
*  sufficiently; that is, if fields of interest are moved.
* @author  J. Bogart
*
* $Header$
*/
  class VolumeIdentifier;

  class TkrId {
  public:

    /** constructor from VolumeIdentifier.  Throws std::invalid_argument
        if @vId is not a valid tkr identifier
        @param vId  Volume identifier, containing complete geometry
        id of volume.  Precise form depends on description
        and nesting of volumes in xml geometry description
    */
    TkrId(const VolumeIdentifier& vId);
    /** constructor including just enough information to identify 
        silicon plane (view optional)     */
    TkrId(unsigned towerX, unsigned towerY, unsigned tray, bool top, 
          int view=eMeasureNone);

    TkrId(const TkrId& id) : m_packedId(id.m_packedId) {}
    TkrId() : m_packedId(0) {};

    ~TkrId() {};

    void copy(const TkrId& id) {m_packedId = id.m_packedId;}
    bool isEqual(const TkrId& other) {
      return ((m_packedId == other.m_packedId));
    }

    const bool operator<(const TkrId& right) const {return m_packedId < right.m_packedId;}
    
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
        eMeasureNone = 2
      };
        
  private:
    enum {
      MASKTowerY = 0x3,
      MASKTowerX = 0x3,
      MASKTray   = 0x1f,
      MASKMeas   = 0x1,
      MASKBotTop = 0x1,
      MASKLadder = 0x3,
      MASKWafer  = 0x3
    };

    enum {
      SHIFTTowerY = 0,
      SHIFTTowerX = 2,
      SHIFTTray   = 4,
      SHIFTMeas   = 9,
      SHIFTBotTop = 10,
      SHIFTLadder = 11,
      SHIFTWafer  = 13
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

    //Access Methods for Tkr reconstruction semantics:
    //unsigned int getLayer() const {return ( getPlane())/2;}
    //unsigned int getPlane() const {return 2*getTray() + getBotTop() - 1;}

    /// Number of valid fields, never more than 7
    //    unsigned int getSize();

    /*
        /// conversion to int
        operator int() const {return m_packedId;}
            
        /// operator to put CalXtalId to output stream
        friend std::ostream& operator<< (std::ostream &stream,
                                             CalXtalId XtalId);
        
        /// operator to get CalXtalId from input stream
        friend std::istream& operator>> (std::istream &stream,
                                            CalXtalId &XtalId);
            
        
        /// function to write to output stream (used by operator <<)
        void write( std::ostream& stream) const;

        /// function to read from input stream (used by operator >>)
        void read( std::istream& stream);
    */      
    void write(std::ostream &stream) const;

  private:

    /// Does the actual work; extracted here since gcc doesn't let
    /// debugger see symbols
    void constructorGuts(const VolumeIdentifier& vId);      

    /// Bitmask containing tracker information
    //    unsigned short int m_packedId;
    unsigned long m_packedId;

    /// Parallel bitmask indicating which fields in m_packedId have been set
    //    unsigned short int m_validFields;

  };
    

    // definition of operator <<
    inline std::ostream& operator<<(std::ostream &stream,
                                      TkrId id)
    { id.write(stream);return stream;}
    

    // definition of operator >>
  /*
    inline std::istream& operator>> (std::istream &stream,
                                        CalXtalId &XtalId)
    {XtalId.read(stream);return stream;}
  */
    
} // namespace idents
#endif    // idents_TKRID_H
