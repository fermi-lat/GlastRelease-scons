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
* Stores packed Xtal ID = (tower*8 + layer)*16 + column .
* Provides access methods to get unpacked ID, i.e. tower, layer, and column.
*                
* Serialize methods give i/o for packed ID.
* Extractor/inserter give i/o for unpacked ID.
*
*    Retrieve packed ID or unpacked tower, layer, and column 
*        inline int getPackedId() const 
*        void getUnpackedId(short& tower, short& layer, short& column); 
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
    class CalXtalId {
        
        
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
        
        /// readout can be either best-of-four energy ranges
        /// or all energy ranges
        typedef enum
        {
            BESTRANGE = 0,
                ALLRANGE = 2
        } CalTrigMode;
        
        /// constructor from packed Id
        CalXtalId(int packedId=0) :
        m_packedId(packedId) {};
        
        /// constructor from tower, layer and column numbers
        CalXtalId(short tower, short layer, short column) 
        {
            packId(tower, layer, column);
        };
            
        ~CalXtalId() {};
                        
        /// get packed ID  
        inline int getPackedId() const {return m_packedId;}

        /// get unpacked tower, layer, and column
        void getUnpackedId(short& tower, short& layer, short& column);
            
        /// get tower 
        inline short getTower() const {return (m_packedId >> 0x7) & 0xf;}

        /// get layer 
        inline short getLayer() const {return (m_packedId >> 0x4) & 0x7;}
        
        /// get column
        inline short getColumn() const {return m_packedId & 0xf;}
            
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
        
        /// Packed word containing Xtal ID = (tower*8 + layer)*16 + column
        unsigned int m_packedId;
        
        /// private method to produce packed Id from tower, layer and column
        inline void packId(short tower, short layer, short column) 
        {m_packedId = (((tower<<3) + layer)<<4) + column;};
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
