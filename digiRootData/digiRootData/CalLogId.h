#ifndef digiRootData_CALLOGID_H
#define digiRootData_CALLOGID_H 1


// Include files
//#include <iostream>
//using namespace std;

#include "TObject.h"


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
class CalLogId : public TObject {

private:

	// Packed word containing log ID = (tower*8 + layer)*16 + column
        Int_t m_packedId;
	Short_t m_tower;
	Short_t m_layer;
	Short_t  m_column;

	inline void unpackId()
	{
		m_tower = getTower();
		m_layer = getLayer();
		m_column = getColumn();
	};
	inline void packId(Short_t tower, Short_t layer, Short_t column) 
		{m_packedId = (((tower<<3) + layer)<<4) + column;};

public:

      CalLogId(Int_t packedId=0) ;
//	{
//		unpackId();
//	};

	CalLogId(Short_t tower, Short_t layer, Short_t column); 
    
    ~CalLogId();

	void init();

	// get packed ID or unpacked tower, layer, and column
         Int_t getPackedId() const;
        void getUnpackedId(Short_t& tower, Short_t& layer, Short_t& column);

	// get tower, layer, and column numbers individually from packed ID
         Short_t getTower() const;
         Short_t getLayer() const;
         Short_t getColumn() const;

        operator int() const {return m_packedId;};

//       	friend ostream& operator<< (ostream &stream, CalLogId logId);
//    	friend istream& operator>> (istream &stream, CalLogId &logId);

//private:
ClassDef(CalLogId,1)
};

/*
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

*/

#endif    // GlastEvent_LOGID_H
