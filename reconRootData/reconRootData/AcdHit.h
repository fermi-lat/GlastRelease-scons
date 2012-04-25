#ifndef ROOT_ACDHIT_H
#define ROOT_ACDHIT_H

#include "TObject.h"
/** 
 *  @class AcdHit
 *
 *  @brief ROOT object which stores calibrated ACD information.  
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
 * \author Eric Charles
 *
 * $Header$
 */

#include "commonRootData/idents/AcdId.h"

class AcdHit : public TObject {

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
  AcdHit();

  /// Constructor for use transient->persistent conversion, 
  /// Takes arguements as they are stored
  AcdHit(const AcdId& id, 
         UShort_t flagsA, UShort_t flagsB, 
         UShort_t phaA, UShort_t phaB, 
         Float_t mipsA, Float_t mipsB)
    :TObject() {
    set(id,flagsA,flagsB,phaA,phaB,mipsA,mipsB);
  }

  /// Destructor is trivial  
  virtual ~AcdHit() {
  }

  /// set everything at once
  void set(const AcdId& id, 
           UShort_t flagsA, UShort_t flagsB, 
           UShort_t phaA, UShort_t phaB, 
           Float_t mipsA, Float_t mipsB){
    m_acdId = id;
    m_flags[0] = flagsA; m_flags[1] = flagsB;
    m_pha[0] = phaA; m_pha[1] = phaB;
    m_mipsPmt[0] = mipsA; m_mipsPmt[1] = mipsB;
  }

  /// Returns the id of the tile or ribbon
  inline const AcdId& getId() const { return m_acdId; }

  /// Return the calibrated pulse height ( relative to MIP signal ) for one PMT
  inline Float_t getMips( PmtId id ) const {
    return m_mipsPmt[id];
  }
  
  /// Return the calibrated pusle height ( relative to MIP signal ) 
  ///  This is averaged between the two PMTs if they are both OK, 
  ///  otherwise only to good one is taken.
  ///  Zero values are _NOT_ suppresed
  inline Float_t getMips( ) const {
    unsigned int nVal(0);
    float val(0.);
    if ( ! getPmtDead(A) && ! getPmtHot(A) ) { nVal++; val += m_mipsPmt[A]; }      
    if ( ! getPmtDead(B) && ! getPmtHot(B) ) { nVal++; val += m_mipsPmt[B]; }
    switch ( nVal ) {
    case 0: val = -1.; break;
    case 1: break;
    case 2: val /= 2.; break;
    }
    return val;
  }

  /// Return the Energy for tiles
  inline Float_t getTileEnergy() const {
    if ( ! m_acdId.isTile() ) return -1.;
    static const float MeVMipTile10 = 1.9;
    static const float MeVMipTile12 = 2.28;
    float MeVMip = m_acdId.isTop() && m_acdId.getRow() == 2 ? MeVMipTile12 : MeVMipTile10;
    return getMips() * MeVMip;
  }

  /// Return the Energy for either PMT on ribbons
  inline Float_t getRibbonEnergy(PmtId id) const {
    static const float MeVMipRibbon = 0.5;
    if ( ! m_acdId.isRibbon() ) return -1;
    return getMips(id) * MeVMipRibbon;            
  }

  /// Return the PHA value
  inline UShort_t getPha(PmtId id) const { 
    return m_pha[id];
  }

  /// Return the flags bits (all at once)
  inline UShort_t getFlags(PmtId id) const { 
    return m_flags[id];
  }

  /// Returns true if pmt flagged as ghost
  inline Bool_t getGhost(PmtId id) const {
    static const Float_t GhostThreshold = 0.5;
    return (getMips(id) > GhostThreshold) && !getHitMapBit(id);
  }
  
  /// Returns true if hit flagged as ghost
  inline bool getGhost( ) const {
    return ( getGhost(A) && getGhost(B) );
  }

  /// Returns true if trigger veto bit set for either PMT
  inline bool getTriggerVeto( ) const {
    return ( getHitMapBit(A) || getHitMapBit(B) );
  }

  /// Denotes that the PMT was above accept (AKA zero-suppresion) threshold
  inline Bool_t getAcceptMapBit(PmtId id) const { 
    return (m_flags[id] & PMT_ACCEPT_MASK) != 0;
  };
  
  /// Denotes that the PMT was above hit (veto) threshold
  inline Bool_t getHitMapBit(PmtId id) const { 
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
  
  /// Denotes that the PMT was DEAD
  inline Bool_t getPmtDead(PmtId id) const { 
    return (m_flags[id] & PMT_DEAD_MASK) != 0;
  };
  
  /// Denotes that the PMT was HOT
  inline Bool_t getPmtHot(PmtId id) const { 
    return (m_flags[id] & PMT_HOT_MASK) != 0;
  };
  
  /// Denotes that there is an error associated w/ the PMT (could be Parity, HOT or DEAD)
  inline bool getPmtError(PmtId id) const { 
    return (m_flags[id] & PMT_ANY_ERROR_MASK) != 0;
  };
  
  /// Print out this structure
  virtual void Print(Option_t* option) const;

protected:

  /// Reset all the values to their defaults
  void clear() {
    m_acdId = AcdId();
    m_flags[0] = m_flags[1] = 0;
    m_pha[0] = m_pha[1] = 0;
    m_mipsPmt[0] = m_mipsPmt[1] = 0.;
  }
  
private:
  
  /// ID of hit tile
  AcdId m_acdId;
  
  /// Various flags about the PMT, see the enum above for their definitions
  UShort_t m_flags[nPmt];
  
  /// The raw PHA values
  UShort_t m_pha[nPmt];
  
  /// The calibrated PHA values
  Float_t m_mipsPmt[nPmt];
  
  ClassDef(AcdHit,1)

};
#endif
