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
            
    ~TkrId() {};

    bool isEqual(const TkrId& other) {
      return ((m_packedId == other.m_packedId) && 
              (m_validFields == other.m_validFields));
    }
    
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
        eMeasureY = 1
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

                        
    bool hasTowerX() {return ((m_validFields & SHMASKTowerX) != 0);}
    unsigned int getTowerX() {
      if (!(hasTowerX())) throw std::domain_error("No TowerX field");
      return (m_packedId & SHMASKTowerX) >> SHIFTTowerX;
    }

    bool hasTowerY() {return ((m_validFields & SHMASKTowerY) != 0);}
    unsigned int getTowerY()  {
      if (!(hasTowerY())) throw std::domain_error("No TowerY field");
      return (m_packedId & SHMASKTowerY) >> SHIFTTowerY;
    }


    bool hasTray() {return ((m_validFields & SHMASKTray) != 0);}
    unsigned int getTray() {
      if (!(hasTray())) throw std::domain_error("No Tray field");
      return (m_packedId & SHMASKTray) >> SHIFTTray;
    }


    bool hasBotTop() {return ((m_validFields & SHMASKBotTop) != 0);}
    unsigned int getBotTop() {
      if (!(hasBotTop())) throw std::domain_error("No BotTop field");
      return (m_packedId & SHMASKBotTop) >> SHIFTBotTop;
    }


    bool hasView() {return ((m_validFields & SHMASKMeas) != 0);}
    unsigned int getView() {
      if (!(hasView())) throw std::domain_error("No View field");
      return (m_packedId & SHMASKMeas) >> SHIFTMeas;
    }


    bool hasLadder() {return ((m_validFields & SHMASKLadder) != 0);}
    unsigned int getLadder() {
      if (!(hasLadder())) throw std::domain_error("No Ladder field");
      return (m_packedId & SHMASKLadder) >> SHIFTLadder;
    }

    
    bool hasWafer() {return ((m_validFields & SHMASKWafer) != 0);}
    unsigned int getWafer() {
      if (!(hasWafer())) throw std::domain_error("No Wafer field");
      return (m_packedId & SHMASKWafer) >> SHIFTWafer;
    }


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
    /// Unuseful default constructor
    TkrId() : m_packedId(0), m_validFields(0) {};

    /// Does the actual work; extracted here since gcc doesn't let
    /// debugger see symbols
    void constructorGuts(const VolumeIdentifier& vId);      

    /// Bitmask containing tracker information
    unsigned short int m_packedId;

    /// Parallel bitmask indicating which fields in m_packedId have been set
    unsigned short int m_validFields;

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
