
/*!
//------------------------------------------------------------------------------
//
// ClassName:   CalXtalId        
//  
// Description: ID class for CAL Xtals                                
//		Supports both packed Xtal ID = (tower*8 + layer)*16 + column and
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
// Include files
#include "idents/CalXtalId.h"
using namespace std;


using namespace idents; 


// overload the inserter to stream unpacked tower, layer and column
void CalXtalId::write(ostream &stream) const
{
        stream << m_tower << " ";
	stream << m_layer << " ";
	stream << m_column << " ";
}


// extract unpacked ID, and stuff packed ID from unpacked info
void CalXtalId::read(istream &stream)
{
	stream >> m_tower >> m_layer >> m_column;
	packId(m_tower, m_layer, m_column);		// ensure m_packedId is synced
}


/// retrieve unpacked ID, tower, layer, and column
void CalXtalId::getUnpackedId(short& tower, short& layer, short& column)
{
	tower  = m_tower;
	layer  = m_layer;
	column = m_column;
}


