#ifndef GlastEvent_LOGID_H
#define GlastEvent_LOGID_H 1


// Include files
#include <iostream>


/*!
//------------------------------------------------------------------------------
//
// ClassName:   CalLogId        
//  
// Description: ID class for CAL logs                                
//		Supports both packed log ID = (tower*8 + layer)*16 + column and
//			unpacked ID, i.e. tower, layer, and column.
//		Serialize methods give i/o for packed ID.
//		Extractor/inserter give i/o for unpacked ID.
//		Input methods ensure that packed and unpacked IDs stay in sync.
//
//    Retrieve packed ID or unpacked tower, layer, and column 
//        inline int getPackedId() const 
//        void getUnpackedId(short& tower, short& layer, short& column); 
// 
//    Retrieve tower, layer, and column numbers individually from packed ID 
//        inline short getTower() const 
//        inline short getLayer() const 
//        inline short getColumn() const 

//              
// Author:  J. Eric Grove	22 Mar 2001                      
//------------------------------------------------------------------------------
 */
namespace idents {
class CalLogId {

private:

	// Packed word containing log ID = (tower*8 + layer)*16 + column
        int m_packedId;
	short m_tower, m_layer, m_column;

	inline void unpackId()
	{
		m_tower = getTower();
		m_layer = getLayer();
		m_column = getColumn();
	};
	inline void packId(short tower, short layer, short column) 
		{m_packedId = (((tower<<3) + layer)<<4) + column;};

public:

      CalLogId(int packedId=0) :
	        m_packedId(packedId)
	{
		unpackId();
	};

	CalLogId(short tower, short layer, short column) 
	{
		packId(tower, layer, column);
		m_tower = tower;
		m_layer = layer;
		m_column = column;
	};
    
    ~CalLogId() {};

	void init();

	// get packed ID or unpacked tower, layer, and column
        inline int getPackedId() const {return m_packedId;};
        void getUnpackedId(short& tower, short& layer, short& column);

	// get tower, layer, and column numbers individually from packed ID
        inline short getTower() const {return (m_packedId >> 0x7) & 0xf;};
        inline short getLayer() const {return (m_packedId >> 0x4) & 0x7;};
        inline short getColumn() const {return m_packedId & 0xf;};

        operator int() const {return m_packedId;};

       	friend ostream& operator<< (ostream &stream, CalLogId logId);
    	friend istream& operator>> (istream &stream, CalLogId &logId);

};


// overload the inserter to stream unpacked tower, layer and column
ostream& operator<<(ostream &stream, CalLogId logId)
{
        stream << logId.m_tower << " ";
	stream << logId.m_layer << " ";
	stream << logId.m_column << " ";
	return stream;
}


// extract unpacked ID, and stuff packed ID from unpacked info
istream& operator>> (istream &stream, CalLogId &logId)
{
	stream >> logId.m_tower >> logId.m_layer >> logId.m_column;
	logId.packId(logId.m_tower, logId.m_layer, logId.m_column);		// ensure m_packedId is synced
	return stream;
}


/// retrieve unpacked ID, tower, layer, and column
void CalLogId::getUnpackedId(short& tower, short& layer, short& column)
{
	tower  = m_tower;
	layer  = m_layer;
	column = m_column;
}




} // namespace idents
#endif    // GlastEvent_LOGID_H
