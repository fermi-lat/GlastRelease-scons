
#ifndef GlastEvent_CalDigi_H
#define GlastEvent_CalDigi_H 1


// Include files
#include <iostream>
#include <vector>
#include "idents/CalLogId.h"

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"
#include "GlastEvent/TopLevel/Definitions.h"


/*!
//------------------------------------------------------------------------------
//
// \class   CalLogReadout        
//  
// \brief   Pulse heights and energy range for both faces of log for Cal                                
//              
// Author:  J. Eric Grove, 23 Feb 2001
//
//------------------------------------------------------------------------------
 */
class CalLogReadout {  // : virtual public ContainedObject  { 

public:
        CalLogReadout(char rangeP, short adcP, char rangeM, short adcM) :
	        m_rangeP(rangeP), 
		m_adcP(adcP), 
		m_rangeM(rangeM), 
		m_adcM(adcM)
        {};

        /// Destructor
        ~CalLogReadout() {};

	/// log ends are labeled by POSitive or NEGative face
        typedef enum
	{
                POS = 0,
                NEG
        } LogFace;

	// retrieve pulse height from specified face
        inline short getAdc(LogFace face) const {return face == POS ? m_adcP : m_adcM;};

	// retrieve energy range from specified face
	inline char getRange(LogFace face) const {return face == POS ? m_rangeP : m_rangeM;};

/*        /// Serialize the object for writing
        inline StreamBuffer& serialize( StreamBuffer& s ) const
	{
		s = ContainedObject::serialize(s);
		return s << m_rangeP 
				 << m_adcP 
				 << m_rangeM 
				 << m_adcM;
	}
        /// Serialize the object for reading
        inline StreamBuffer& serialize( StreamBuffer& s )
	{
		s = ContainedObject::serialize(s);
		s >> m_rangeP
		  >> m_adcP
		  >> m_rangeM
		  >> m_adcM;
		return s;
	}
*/

private:

        short m_adcP, m_adcM;
	char  m_rangeP, m_rangeM;

};




/*!
//------------------------------------------------------------------------------
//
// \class   CalDigi        
//  
// \brief Digitizations for Cal                                
//              
// Author:  J. Eric Grove, 23 Feb 2001
//
//------------------------------------------------------------------------------
 */

extern const CLID& CLID_CalDigi;

class CalDigi : virtual public ContainedObject  { 

public:
	/// each log end can report four energy ranges
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

	/// shifts and masks for packed readout of energy range and Adc value
        enum {POS_OFFSET = 14,						// shift for POSitive face
		RANGE_OFFSET = 12, RANGE_MASK = 0x3000,		// energy range bits
                ADC_VAL_MASK = 0xfff};						// Adc value bits

        CalDigi() {};

        //CalDigi(CalTrigMode mode, idents::CalLogId CalLogId, ObjectVector<CalLogReadout> readout) : 
	    //    m_mode(mode),
        //        m_logId(CalLogId),
        //        m_readout(readout)
        //{};

        /// Destructor
        virtual ~CalDigi() { };

        //! Retrieve reference to class definition structure
        virtual const CLID& clID() const   { return CalDigi::classID(); }
        static const CLID& classID()       { return CLID_CalDigi; }

	/// Retrieve readout mode
	inline const CalTrigMode getMode() const { return m_mode; };
    inline void setMode(CalTrigMode m) { m_mode = m; };

	/// Retrieve log identifier
        inline const idents::CalLogId getPackedId() const { return m_logId; };
        inline void setPackedId(idents::CalLogId id) { m_logId = id; };

        inline void addReadout(CalLogReadout r) { m_readout.push_back(r); } ;
	
	/// Retrieve energy range for selected face and readout
	inline char getRange(short readoutIndex, CalLogReadout::LogFace face) const
	{
		return (readoutIndex < m_readout.size()) ? ((m_readout[readoutIndex])).getRange(face) : (char)-1;
	}

	/// Retrieve pulse height for selected face and readout
	inline short getAdc(short readoutIndex, CalLogReadout::LogFace face) const
	{
		return (readoutIndex < m_readout.size()) ? ((m_readout[readoutIndex])).getAdc(face) : (short)-1;
	}

	/// Retrieve ranges and pulse heights from both ends of selected readout
	inline const CalLogReadout* getLogReadout(short readoutIndex)
	{
 		//return ((readoutIndex < m_readout.size()) ? m_readout[readoutIndex] : 0);
        if ( readoutIndex < m_readout.size() )
            return &(m_readout[readoutIndex]);
        else
            return 0;

	}

	/// Retrieve pulse height from selected range
	inline short getAdcSelectedRange(char range, CalLogReadout::LogFace face) const
	{
		char nRanges = (char)m_readout.size();
		if (nRanges == 1)
			return (range == ((m_readout[0])).getRange(face)) ? ((m_readout[0])).getAdc(face) : (short)-1;
		else
			return ((m_readout[(nRanges + range - ((m_readout[0])).getRange(face)) % nRanges])).getAdc(face);
	}

        /// Serialize the object for writing
        virtual StreamBuffer& serialize( StreamBuffer& s ) const;
        /// Serialize the object for reading
        virtual StreamBuffer& serialize( StreamBuffer& s );
        /// Fill the ASCII output stream
        virtual std::ostream& fillStream( std::ostream& s ) const;

private:

        /// Cal readout mode is based on trigger type
        CalTrigMode m_mode;
        /// Cal ID
        idents::CalLogId m_logId;
        /// ranges and pulse heights
        //ObjectVector<CalLogReadout> m_readout;
        std::vector<CalLogReadout> m_readout;
       // std::vector<int> m_readout;

};

//! Definition of all container types of CalDigi
//	//m_CalRawLogs = SmartDataPtr<CalDigiVector>(eventSvc(),"/Event/Digi/CalDigis"); 
//template <class TYPE> class ObjectVector;
typedef ObjectVector<CalDigi> CalDigiVector;
typedef ObjectList<CalDigi> CalDigiList;


/// Serialize the object for writing
inline StreamBuffer& CalDigi::serialize( StreamBuffer& s ) const
{
	ContainedObject::serialize(s);
	//s << m_mode;
	//s = m_logId.serialize(s);
	//s = (*(m_readout[0])).serialize(s);
	//if (m_mode == ALLRANGE)
	//	for (int rangeIndex=1; rangeIndex<4; rangeIndex++)
//			s = (*(m_readout[rangeIndex])).serialize(s);
  
	return s;
}


/// Serialize the object for reading
inline StreamBuffer& CalDigi::serialize( StreamBuffer& s )
{
	ContainedObject::serialize(s);
//        int mode;
  //      s >> mode; m_mode = (mode==CalTrigMode::ALLRANGE) ? ALLRANGE:BESTRANGE;
// For now by default m_mode is considered to be BESTRANGE
        
    //    s = m_logId.serialize(s);
	//s = (*(m_readout[0])).serialize(s);
	//if (m_mode == ALLRANGE)
	//	for (int rangeIndex = 1; rangeIndex < 4; rangeIndex++)
	//		s = (*(m_readout[rangeIndex])).serialize(s);

	return s;
}


/// Fill the ASCII output stream
inline std::ostream& CalDigi::fillStream( std::ostream& s ) const
{
    /*
    s << "    base class CalDigi :"
    << "\n        CalTrigMode = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_mode << " )"
    << "\n        ID = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << m_logId << " )"
    << "\n        Best plus-face range and pulseheight  = ( "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << CalDigi::getRange(BESTRANGE,CalLogReadout::LogFace::POS) << ", " << CalDigi::getAdc(BESTRANGE,CalLogReadout::LogFace::POS) << " )"
    << "\n        Best minus-face range and pulseheight = "
    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
    << CalDigi::getRange(BESTRANGE,CalLogReadout::LogFace::NEG) << ", " << CalDigi::getAdc(BESTRANGE,CalLogReadout::LogFace::NEG) << " )";

	if (m_mode == ALLRANGE)
	{
		for (int rangeIndex = 1; rangeIndex < 4; rangeIndex++)
		{
			s << "\n        Next plus-face range and pulseheight  = "
		    << GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
			<< CalDigi::getRange(rangeIndex,CalLogReadout::LogFace::POS) << ", " << CalDigi::getAdc(rangeIndex,CalLogReadout::LogFace::POS) << " )"
			<< "\n        Next minus-face range and pulseheight = "
			<< GlastEventFloatFormat( GlastEvent::width, GlastEvent::precision )
			<< CalDigi::getRange(rangeIndex,CalLogReadout::LogFace::NEG) << ", " << CalDigi::getAdc(rangeIndex,CalLogReadout::LogFace::NEG) << " )";
		}
	}

	s << " )\n";
*/
	return s;
}


#endif