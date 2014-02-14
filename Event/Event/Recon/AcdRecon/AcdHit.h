#ifndef Event_ACDHIT_H
#define Event_ACDHIT_H

#include <vector>

#include "Event/TopLevel/Definitions.h"
#include "idents/AcdId.h"

#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/DataObject.h"
#include "GaudiKernel/IInterface.h"

#include "Event/Digi/AcdDigi.h"

class MsgStream;

static const CLID& CLID_AcdHitCol = InterfaceID("AcdHitCol", 1, 0);

/**
 *  @class Event::AcdHit
 *
 *  @brief This ROOT object stores the calibrated Acd information.  
 *
 *  There is a 1-1 correspondance between AcdHits and AcdDigis.  The hits also 
 *  add the pulse height information expressed in terms of MIPs.  To preserve the
 *  correspondance with AcdDigi no correlation is made with tracking.  This means 
 *  we don't correct for pathlength through the element.
 *  
 *  The main access functions are:
 *  - Functions equivalent to AcdDigi
 *    - const idents::AcdId& getAcdId()  
 *      - which returns the ID of the hit element
 *    - unsigned short getPha(PmtId id) 
 *      - which returns the raw PHA value for either PMT
 *    - unsigned short getFlags(PmtId id)
 *      - which returns a mask with the status bits associated with either PMT
 *         Those status bits are:
 *         -      PMT_ACCEPT_BIT = 0,              // channel is above zero suppresion threshold (applied to digital data)
 *         -      PMT_VETO_BIT = 1,                // channel fired veto discriminator (applied to analog data)
 *         -      PMT_RANGE_BIT = 2,               // channel was read out in high range
 *         -      PMT_CNO_BIT = 3,                 // could be any channel on that GARC
 *         -      PMT_ODD_PARITY_ERROR_BIT = 4,    // PHA data transmission had parity error
 *         -      PMT_HEADER_PARITY_ERROR_BIT = 5, // header data transmission had parity error 
 *         -      PMT_DEAD_BIT = 6,                // PMT was flagged as dead by offline calibration
 *         -      PMT_HOT_BIT = 7                  // PMT was flagged as hot by offline calibration
 *
 *  - Functions with calibrated data
 *    - float(PmtId id) const
 *      - which returns the calibrated pulse height expressed in terms of MIPS of either PMT
 *    - float mips( ) 
 *      - which returns the average calibrated pulse height expressed in terms of MIPS of both PMTs
 *
 *  \author Eric Charles
 *
 * $Header$
 **/

namespace Event
{
  
  class AcdHit : virtual public ContainedObject {
    
  public:

    // copied from AcdDigi
    typedef enum {
      A = 0,
      B = 1,
      nPmt
    } PmtId;
    
    // copied from AcdDigi
    typedef enum {
      NOERROR = 0,
      ERROR = 1
    } ParityError;
    
    // order of bits in bitMasks
    typedef enum {
      PMT_ACCEPT_BIT = 0,              // channel is above zero suppresion threshold (applied to digital data)
      PMT_VETO_BIT = 1,                // channel fired veto discriminator (applied to analog data)
      PMT_RANGE_BIT = 2,               // channel was read out in high range
      PMT_CNO_BIT = 3,                 // could be any channel on that GARC
      PMT_ODD_PARITY_ERROR_BIT = 4,    // PHA data transmission had parity error
      PMT_HEADER_PARITY_ERROR_BIT = 5, // header data transmission had parity error 
      PMT_DEAD_BIT = 6,                // PMT was flagged as dead by offline calibration
      PMT_HOT_BIT = 7                  // PMT was flagged as hot by offline calibration
    } FlagBits;
    
    // masks used to test conditions
    typedef enum { 
      PMT_ACCEPT_MASK = 0x1,              // channel is above zero suppresion threshold (applied to digital data)
      PMT_VETO_MASK = 0x2,                // channel fired veto discriminator (applied to analog data)
      PMT_RANGE_MASK = 0x4,               // just a trick to split the errors into the higher byte
      PMT_CNO_MASK = 0x8,                 // could be any channel on that GARC
      PMT_ODD_PARITY_ERROR_MASK = 0x10,   // PHA data transmission had parity error
      PMT_HEADER_PARITY_ERROR_MASK = 0x20,// header data transmission had parity error 
      PMT_DEAD_MASK = 0x40,               // PMT was flagged as dead by offline calibration
      PMT_HOT_MASK = 0x80,                // PMT was flagged as hot by offline calibration
      PMT_ANY_ERROR_MASK = 0xF0           // PMT has any error
    } FlagMasks;

  public:
    
    /// Default constructor
    AcdHit(){
      ini();
    }

    /// Constructor for use in reconstruction, takes digi and calibrated values
    AcdHit(const Event::AcdDigi& digi, float mipsPmtA, float mipsPmtB);

    /// Constructor for use in persistent -> transient conversion, 
    /// Takes arguements as they are stored in ROOT
    AcdHit(const idents::AcdId&, unsigned short flagsA, unsigned short flagsB, 
           unsigned short phaA, unsigned short phaB,
           float mipsPmtA, float mipsPmtB);
    
    /// Destructor is trivial
    virtual ~AcdHit() {};
    
    /// Direct access to parameters

    /// set everything at once
    void set(const idents::AcdId&, unsigned short flagsA, unsigned short flagsB, 
             unsigned short phaA, unsigned short phaB,
             float mipsPmtA, float mipsPmtB);
    
    /// this is to allow us to kill accept map bits for periodic triggers in the overlays
    void correctAcceptMapBits(bool acceptA, bool acceptB);

    /// Returns the id of the tile or ribbon
    inline const idents::AcdId& getAcdId() const { return m_acdId; };

    /// Return the calibrated pulse height ( relative to MIP signal ) for one PMT
    inline float mips( PmtId id ) const {
      return m_mipsPmt[id];
    }
    
    /// Return the calibrated pusle height ( relative to MIP signal ) 
    ///  This is averaged between the two PMTs if they are both OK, 
    ///  otherwise only to good one is taken.
    ///  Zero values are _NOT_ suppresed
    inline float mips( ) const {
      unsigned int nVal(0);
      float val(0.);
      if ( ! getPmtError(A) ) { nVal++; val += m_mipsPmt[A]; }      
      if ( ! getPmtError(B) ) { nVal++; val += m_mipsPmt[B]; }
      switch ( nVal ) {
      case 0: val = -1.; break;
      case 1: break;
      case 2: val /= 2.; break;
      }
      return val;
    }


    /// Return the Energy for tiles
    inline float tileEnergy() const {
      if ( ! m_acdId.tile() ) return -1.;
      static const float MeVMipTile10 = 1.9;
      static const float MeVMipTile12 = 2.28;
      float MeVMip = m_acdId.top() && m_acdId.row() == 2 ? MeVMipTile12 : MeVMipTile10;
      return mips() * MeVMip;
    }

    /// Return the Energy for either PMT on ribbons
    inline float ribbonEnergy(PmtId id) const {
      static const float MeVMipRibbon = 0.5;
      if ( ! m_acdId.ribbon() ) return -1;
      return mips(id) * MeVMipRibbon;            
    }

    /// Returns the raw PHA value.  
    /// Note that you also should check the range bit if you want to use this
    inline unsigned short getPha(PmtId id) const { 
      return m_pha[id];
    }

    /// Returns true if pmt flagged as ghost
    inline bool getGhost(PmtId id) const {
      static const float GhostThreshold = 0.5;
      return (mips(id) > GhostThreshold) && !getHitMapBit(id);
    }

    /// Returns true if hit flagged as ghost by both PMTs
    inline bool getGhost( ) const {
      return ( getGhost(A) && getGhost(B) );
    }

    /// Returns true if trigger veto bit set for either PMT
    inline bool getTriggerVeto( ) const {
      return ( getHitMapBit(A) || getHitMapBit(B) );
    }

    /// Returns all the flags at once as a bit mask
    inline unsigned short getFlags(PmtId id) const { 
      return m_flags[id];
    }

    /// Denotes that the PMT was above accept (AKA zero-suppresion) threshold
    inline bool getAcceptMapBit(PmtId id) const { 
      return (m_flags[id] & PMT_ACCEPT_MASK) != 0;
    };

    /// Denotes that the PMT was above hit (veto) threshold
    inline bool getHitMapBit(PmtId id) const { 
      return (m_flags[id] & PMT_VETO_MASK) != 0;
    };

    /// Denotes that a parity error was seen in transmitting the PHA data for this channel
    inline ParityError getOddParityError(PmtId id) const {
      return (m_flags[id] & PMT_ODD_PARITY_ERROR_MASK) != 0 ? ERROR : NOERROR;
    }
    
    /// Denotes that a parity error was seen in transmitting the header data in this event
    inline ParityError getHeaderParityError(PmtId id) const { 
      return (m_flags[id] & PMT_HEADER_PARITY_ERROR_MASK) != 0 ? ERROR : NOERROR;
    };
    
    /// Denotes that the PMT was flagged as DEAD by the offline calibrations
    inline bool getPmtDead(PmtId id) const { 
      return (m_flags[id] & PMT_DEAD_MASK) != 0;
    };

    /// Denotes that the PMT was flagged as HOT by the offline calibrations
    inline bool getPmtHot(PmtId id) const { 
      return (m_flags[id] & PMT_HOT_MASK) != 0;
    };

    /// Denotes that there is an error associated w/ the PMT (could be Parity, HOT or DEAD)
    inline bool getPmtError(PmtId id) const { 
      return (m_flags[id] & PMT_ANY_ERROR_MASK) != 0;
    };

    /// Print out this structure on a stream
    virtual void writeOut(MsgStream& stream) const;
    
  protected:

    /// Reset all the values to their defaults
    virtual void ini();

    /// Parse out the flags from the AcdDigi
    void setFlags(const Event::AcdDigi& digi);
    
  private:
    
    /// ID of hit tile
    idents::AcdId m_acdId;

    /// Various flags about the PMT.  See the enum above for their definitions
    unsigned short m_flags[nPmt];

    /// The raw PHA values
    unsigned short m_pha[nPmt];
    
    /// The calibrated PHA values
    float m_mipsPmt[nPmt];  
    
  };

   
  /*! 
   * @class AcdHitCol
   *
   *  @brief TDS class to store the set of AcdHits
   *  
   * It inherits from DataObject
   * (to be a part of Gaudi TDS) and from std::vector of pointers
   * to AcdHit objects. Some methods just rename the std::vector
   * methods with the same functionality for backward compartibility.
   *
   * @author Eric Charles
   *
   * @todo replace this class by typedef to ObjectVector, make corresponding
   *       changes in AcdPha2MipTool
   */
    

  class AcdHitCol : public DataObject, public std::vector<AcdHit*> 
  {
  public:

    /// Default constructor.  Builds empty collection
    AcdHitCol() { clear();}

    /// "copy" constructor.  Take ownerships of a vector of AcdHits
    AcdHitCol(const std::vector<AcdHit*>& acdhits);
    
    /// destructor - deleting the hits pointed
    /// by the vector elements
    ~AcdHitCol() { delHits();}
            
    // GAUDI members to be use by the converters
    static const CLID& classID() {return CLID_AcdHitCol;}
    virtual const CLID& clID() const {return classID();}
    
    // take ownership
    void init(std::vector<AcdHit*>& other) {
      for ( std::vector<AcdHit*>::iterator itr = other.begin(); itr != other.end(); itr++ ) {
        push_back(*itr);
      }
    }

    /// Add a new hit
    void add(AcdHit* cl) {push_back(cl);}
    
    /// get the number of hits in collection
    int num()                  const {return size();}
    
    /// get pointer to the hit with given number 
    AcdHit * getHit(int i) const {return operator[](i);}
    
    /// delete all hits pointed by the vector elements
    void delHits();
    
    /// write information for all hits to the ascii file 
    /// for debugging purposes
    virtual void writeOut(MsgStream& stream) const;
    
  protected:
    
    /// does the same function as clear() 
    virtual void ini();
        
  };

}

#endif
