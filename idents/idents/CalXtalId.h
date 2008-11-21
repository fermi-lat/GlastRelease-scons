#ifndef idents_CALXTALID_H
#define idents_CALXTALID_H 1


// Include files
#include <iostream>


namespace idents {

/**
*
* @class   CalXtalId        
*  
* @brief ID class for CAL Xtals
*                                
* Stores xtal id as a packed word of 32 bits.
* Provides access methods to get unpacked ID, i.e. tower, layer, and column.
*                
* Serialize methods give i/o for packed ID.
* Extractor/inserter give i/o for unpacked ID.
*
*    Retrieve packed ID or unpacked tower, layer, and column 
*        inline int getPackedId() const 
*        void getUnpackedId(short& tower, short& layer, short& column); 
*
* The ID may also identify face and range, but need not.  Boolean functions
* validFace() and validRange() indicate whether or not the ID contains this
* information
*
* The packedId looks like:
* @verbatim
*     15 14 13  12 11  10 9 8 7  6 5 4  3 2 1 0
*      _  _  _   _  _   _ _ _ _  _ _ _  _ _ _ _
*        Range   Face   Tower   Layer  Column
* BITS:  13-15  11-12    7-10    4-6    0-3
* @endverbatim
*
* where the top bit in the Range and Face fields is used to keep track
* of the validity of the remainder of the field.
* 
*    Retrieve tower, layer, and column numbers individually from packed ID 
*        inline short getTower() const 
*        inline short getLayer() const 
*        inline short getColumn() const 
*
*              
* @author  J. Eric Grove	                     
*
* $Header$
*/
  class VolumeIdentifier;
  
  class CalXtalId {
        
        
  public:
        
    /// Xtal ends are labeled by POSitive or NEGative face
    typedef enum
      {
        FACE_UNUSED = -1,
        POS = 0,
        NEG = 1
      } XtalFace;
        
        
    /// each Xtal end can report four energy ranges
    typedef enum
      {
        RANGE_UNUSED = -1,
        LEX8 = 0,
        LEX1 = 1,
        HEX8 = 2,
        HEX1 = 3
      } AdcRange;
        
    /// readout can be either best-of-four energy ranges
    /// or all energy ranges
    typedef enum
      {
        BESTRANGE = 0,
        ALLRANGE = 2
      } CalTrigMode;
        
    /// diode type  -  Large or Small
    typedef enum
      {
        LARGE = 0,
        SMALL = 1
      } DiodeType;
		
    /// Return diode type, given AdcRange
    static DiodeType rangeToDiode(AdcRange range) {
      if ((range == LEX8) || (range == LEX1)) return LARGE;
      else return SMALL;
    }

    /// constructor from packed Id
    CalXtalId(int packedId=0) :
      m_packedId(packedId) {};
        
    /// constructor from tower, layer and column numbers
    CalXtalId(short tower, short layer, short column, 
              short face=FACE_UNUSED, short range=RANGE_UNUSED)
    {
      packId(tower, layer, column, face, range);
    };
    /** constructor from VolumeId.  Throws std::invalid_argument
        if @vId is not a valid Cal identifier
        @param vId  Volume identifier, containing complete geometry
        id of volume.  Precise form depends on description
        and nesting of volumes in xml geometry description
        @param xNum number of towers in a row, needed to calculate
        single tower number from (row, column) representation
        in @a vid
    */
    CalXtalId(const VolumeIdentifier& vId, unsigned xNum=4);
            
    ~CalXtalId() {}; 

    VolumeIdentifier* makeVolumeId() const;
                        
    /// get packed ID  
    inline int getPackedId() const {return m_packedId;}

    /// get unpacked tower, layer, and column
    void getUnpackedId(short& tower, short& layer, short& column);
            
    /// get tower 
    inline short getTower() const 
    {return (m_packedId >> TOWER_SHIFT) & 0xf;}

    /// get layer 
    inline short getLayer() const 
    {return (m_packedId >> LAYER_SHIFT) & 0x7;}
        
    /// get column
    inline short getColumn() const 
    {return (m_packedId >> COLUMN_SHIFT) & 0xf;}

    bool validFace() const {
      return ( ( (1 << FACE_VALID_SHIFT) & m_packedId) != 0); }
    bool validRange() const {
      return ( ( (1 << RANGE_VALID_SHIFT) & m_packedId) != 0); }

    inline short getFace() const {
      if (validFace()) return (m_packedId >> FACE_SHIFT) & 0x1;
      else return FACE_UNUSED;
    }

    inline short getRange() const {
      if (validRange()) return (m_packedId >> RANGE_SHIFT) & 0x3;
      else return RANGE_UNUSED;
    }

    /// get measurement direction
    inline bool isX() const {return (getLayer()%2 == 0);};
            
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
            
  private:
    enum {
      COLUMN_SHIFT = 0,
      LAYER_SHIFT = 4,
      TOWER_SHIFT = 7,
      FACE_SHIFT = 11,
      FACE_VALID_SHIFT = 12,
      RANGE_SHIFT = 13,
      RANGE_VALID_SHIFT = 15
    };

    /// Packed word containing Xtal ID = (tower*8 + layer)*16 + column
    unsigned int m_packedId;
        
    /// private method to produce packed Id from tower, layer and column
    inline void packId(short tower, short layer, short column,
                       short face, short range) {
      m_packedId = (tower<<TOWER_SHIFT) + (layer<<LAYER_SHIFT) + 
        (column<<COLUMN_SHIFT);
      if (face != FACE_UNUSED) {
        m_packedId += ((1 << FACE_VALID_SHIFT) + (face << FACE_SHIFT));
      }
      if (range != RANGE_UNUSED) {
        m_packedId += ((1 << RANGE_VALID_SHIFT) + 
                       (range<<RANGE_SHIFT));
      }
    };
  };
    

  // definition of operator <<
  inline std::ostream& operator<<(std::ostream &stream,
                                  CalXtalId XtalId)
  { XtalId.write(stream);return stream;}
    

  // definition of operator >>
  inline std::istream& operator>> (std::istream &stream,
                                   CalXtalId &XtalId)
  {XtalId.read(stream);return stream;}
    
    
} // namespace idents
#endif    // idents_CALXTALID_H
