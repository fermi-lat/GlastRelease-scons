#ifndef digiRootData_CalXtalId_H
#define digiRootData_CalXtalId_H 1

#include "TObject.h"


/* @class CalXtalId        
* @brief ID class for CAL logs                                
* Supports both packed log ID = (tower*8 + layer)*16 + column and
* unpacked ID, i.e. tower, layer, and column.
* Serialize methods give i/o for packed ID.
* Extractor/inserter give i/o for unpacked ID.
* Input methods ensure that packed and unpacked IDs stay in sync.
* 
* Retrieve packed ID or unpacked tower, layer, and column 
* inline int getPackedId() const 
* void getUnpackedId(short& tower, short& layer, short& column); 
*
* Retrieve tower, layer, and column numbers individually from packed ID 
* inline short getTower() const 
* inline short getLayer() const 
* inline short getColumn() const 
*             
* The packedId looks like:
* @verbatim
* 11 10 9 8 7  6 5 4  3 2 1 0
*  _  _ _ _ _  _ _ _  _ _ _ _
*      Tower   Layer  Column
* BITS: 7-10    4-6    0-3
* @endverbatim
* @author  J. Eric Grove
* $Header$
*/
class CalXtalId : public TObject {

public:
 
    /// Xtal ends are labeled by POSitive or NEGative face
    typedef enum
    {
        POS = 0,
        NEG
    } XtalFace;

    /// each Xtal end can report four energy ranges
    typedef enum
    {
        LEX8 = 0,
            LEX1,
            HEX8,
            HEX1
    } AdcRange;
    
    /// readout can be either best-of-four energy ranges or all energy ranges
    typedef enum
    {
        BESTRANGE = 0,
            ALLRANGE = 2
    } CalTrigMode;

    CalXtalId(UInt_t packedId=0);
    
    CalXtalId(Short_t tower, Short_t layer, Short_t column); 
    
    virtual ~CalXtalId();
    
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void init(Short_t tower, Short_t layer, Short_t column);
    void init(UInt_t packedId);
    
    /// get packed ID or unpacked tower, layer, and column
    UInt_t getPackedId() const;
    void getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column);
    
    /// get tower, layer, and column numbers individually from packed ID
    Short_t getTower() const;
    Short_t getLayer() const;
    Short_t getColumn() const;
    
    operator int() const { return m_packedId; };
            
private:
    
    /// Packed word containing log ID = (tower*8 + layer)*16 + column
    UInt_t m_packedId;
    
    inline void packId(Short_t tower, Short_t layer, Short_t column) 
    {m_packedId = (((tower<<3) + layer)<<4) + column;};
    

    ClassDef(CalXtalId,1) // Cal Xtal Identifier
};


#endif    // GlastEvent_LOGID_H
