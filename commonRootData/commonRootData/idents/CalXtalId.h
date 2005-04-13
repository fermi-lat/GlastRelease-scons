#ifndef commonRootData_CalXtalId_H
#define commonRootData_CalXtalId_H 1

#include "TObject.h"

namespace commonRootData {

/* @class CalXtalId        
* @brief ID class for CAL logs   
*                             
* Based on CalXtalId in the idents package.
* Stores xtal id as a packed word (UInt_t) of 32 bits.
* Supports both packed log ID = (tower*8 + layer)*16 + column and
* unpacked ID, i.e. tower, layer, and column.
* 
* Retrieve packed ID or unpacked tower, layer, and column 
* inline int getPackedId() const 
* void getUnpackedId(short& tower, short& layer, short& column); 
*
* The ID may also identify face and range, but need not.  Boolean functions
* faceValid() and rangeValid() indicate whether or not the ID contains
* this information  
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
* @author  Heather Kelly
* $Header$
*/
class CalXtalId : public TObject {

public:
 
    /// @enum XtalFace Xtal ends are labeled by POSitive or NEGative face
    typedef enum
    {
        FACE_UNUSED = -1,
        POS = 0,
        NEG = 1
    } XtalFace;

    /// @enum AdcRange Each Xtal end can report four energy ranges
    typedef enum
    {
        RANGE_UNUSED = -1,
        LEX8 = 0,
        LEX1 = 1,
        HEX8 = 2,
        HEX1 = 3
    } AdcRange;
    
    /// @enum CalTrigMode Readout can be either best-of-four energy ranges or all energy ranges
    typedef enum
    {
        BESTRANGE = 0,
        ALLRANGE = 2
    } CalTrigMode;

    CalXtalId();

    CalXtalId(const CalXtalId& id):TObject(id) { m_packedId = id.m_packedId; };

    CalXtalId(UInt_t packedId);
    
    /// Default values of -1 for face and range signify these fields are to
    /// be marked invalid
    CalXtalId(Short_t tower, Short_t layer, Short_t column, 
              Short_t face=FACE_UNUSED, Short_t range=RANGE_UNUSED); 
    
    ~CalXtalId();

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void init(Short_t tower, Short_t layer, Short_t column, 
              Short_t face=FACE_UNUSED, Short_t range=RANGE_UNUSED);
    void init(UInt_t packedId);
    
    /// Returns the packed Id
    UInt_t getPackedId() const;
    /// Returns the unpacked tower, layer, and column
    void getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column);
    
    /// get tower, layer, and column numbers individually from packed ID
    Short_t getTower() const;
    Short_t getLayer() const;
    Short_t getColumn() const;
    Short_t getFace() const;
    Short_t getRange() const;
    
    bool validFace() const {
        return ( ( (1 << FACE_VALID_SHIFT) & m_packedId) != 0); }
    bool validRange() const {
        return ( ( (1 << RANGE_VALID_SHIFT) & m_packedId) != 0); }

    operator int() const { return m_packedId; };
            
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
    /// Packed word containing log ID = (tower*8 + layer)*16 + column
    UInt_t m_packedId;
    
    inline void packId(Short_t tower, Short_t layer, Short_t column,
                       Short_t face, Short_t range) 
    {m_packedId = (((tower<<3) + layer)<<4) + column;
    if (face != FACE_UNUSED) {
        m_packedId += ((1 << FACE_VALID_SHIFT) + (face << FACE_SHIFT));
    }
    if (range != RANGE_UNUSED) {
        m_packedId += ((1 << RANGE_VALID_SHIFT) + (range << RANGE_SHIFT));
    }

};
    

    ClassDef(CalXtalId,1) // Cal Xtal Identifier
};

}

#endif    // GlastEvent_LOGID_H
