#ifndef RECONHEADER_H
#define RECONHEADER_H

#include "TObject.h"

//! Root object containing flags for a concerning reconstruction data 
/*! Currently contains flags denoting whether or not:
    TKR reconstruction data is available
    CAL reconstruction data is available */

class ReconHeader : public TObject
{
private:
    //! Packed word that contains all of the flags for this event
    /*!\verbatim
                    ||RESERVED  ||
                    ||For Future||   
                    || 31....16 ||

||                   CAL FLAGS                          ||
|| 15 | 14 | 13 | 12 | 11      | 10   | 09     | 08     ||
||  UNUSED           |NO ENERGY|POOR  |FAIR    |GOOD    ||
||                   |RECON    |      |        |        ||

||                   TKR FLAGS                          ||
|| 07 | 06 | 05 | 04 | 03       | 02  | 01     | 00     ||
||    UNUSED         |MANYTracks|GAMMA|NOTRACKS|TKRRECON|| 

\endverbatim */

    UInt_t m_flags;

    enum {
    //! offset to the beginning of the CAL flags
        calOffset = 8,
    //! offset to the beginning of the TKR flags
        tkrOffset = 0
    };

    //! Event Id
    UInt_t m_eventID;

public:
    //! Defines the possible values for the tkrRecon bits
    typedef enum {
        TKRRECON = UInt_t(1 << tkrOffset),     // 1 == Tkr Recon data is available
        NOTRACKS = UInt_t(1 << (tkrOffset+1)),  // 1 == No tracks were reconstructed 
        GAMMA = UInt_t(1 << (tkrOffset+2)),      // bit is on (1) when there was a reconstructed gamma
        MANYTRACKS = UInt_t(1 << (tkrOffset+3))  // bit is on (1) when there were >1 tracks found
    } tkrFlags;

    //! Defines the possible values for the calRecon bits
    typedef enum {
	//! set to 1 when the energy recon is "Good"
        GOOD = UInt_t(1 << calOffset),  
	//! set to 1 when energy recon is "Fair"
        FAIR = UInt_t(1 << (calOffset+1)),  
	//! set to 1 when the energy recon is "Poor"
        POOR = UInt_t(1 << (calOffset+2)), 
	//! set to 1 when NO energy recon is available
        NOENERGY = UInt_t(1 << (calOffset+3))
    } calFlags;

    enum {
	//! extract TKR data
        tkrTestWord = UInt_t(TKRRECON | NOTRACKS | GAMMA | MANYTRACKS), 
	//! extract CAL data
        calTestWord = UInt_t(GOOD | FAIR | POOR | NOENERGY) 
    };

    //! default constructor - initialize all flags to false
    ReconHeader();

    //! destructor
    ~ReconHeader() {Clear();};
    
    //! Reset all flags to false
    int Clear();

    //! get the word containing the flags
    UInt_t getFlags() { return m_flags; };

    //! Returns the integer value of the TKR flags
    UInt_t getTkrFlags() { return (m_flags & tkrTestWord); };

    // Returns the integer value of the CAL flags
    UInt_t getCalFlags(){ return (m_flags & calTestWord); };

    //! denotes whether there is ANY recon data available
    /*! depends on whether there is any TKR or CAL recon data */
    Bool_t HasRecon() { return (((m_flags & TKRRECON) || !(m_flags & NOENERGY)) ? kTRUE : kFALSE); };

    //! is there any TKR recon data available?
    Bool_t HasTkrRecon() { return ((m_flags & TKRRECON) ? kTRUE : kFALSE); };

    //! Is there a reconstructed gamma?
    Bool_t HasGamma() { return ((m_flags & GAMMA) ? kTRUE : kFALSE); };

    //! Did the TKR recon find more than 1 track?
    Bool_t HasManyTracks() { return ((m_flags & MANYTRACKS) ? kTRUE : kFALSE); };

    //! Is there CAL reconstruction data? Just checks that NOENERGY bit is not set
    Bool_t HasCalRecon() { return (!(m_flags & NOENERGY) ? kTRUE : kFALSE);};

    //! Was the energy reconstruction good?
    Bool_t HasCalGood() { return ((m_flags & GOOD) ? kTRUE : kFALSE);};

    //! Was energy reconstruction fair?
    Bool_t HasCalFair() { return ((m_flags & FAIR) ? kTRUE : kFALSE);};

    //! Was the energy reconstruction poor?
    Bool_t HasCalPoor() { return ((m_flags & POOR) ? kTRUE : kFALSE);};

    //! True when NOENERGY bit is set - meaning NO Energy reconstruction
    Bool_t HasCalNoEnergyRecon() { return ((m_flags & NOENERGY) ? kTRUE : kFALSE);};

    //! set the TKR flags
    /*! returns kFALSE if the routine determines the val provided is out of
        bounds...if successful kTRUE is returned */
    Bool_t setTkrFlag(UInt_t val);

    //! set the CAL flags
    /*! returns kFALSE if the routine determines the val provided is out of
        bounds...if successful kTRUE is returned */
    Bool_t setCalFlag(UInt_t val);

    //! Checks to be sure flag falls within bounds for TKR bits
    Bool_t isValidTkrFlag(UInt_t val);

    //! Checks to be sure flag falls within bounds for CAL bits
    Bool_t isValidCalFlag(UInt_t val);

    //! retrieve the event ID for this event
    UInt_t getEventID() { return m_eventID; };

    //! set the event ID for this event - should match the event ID in the raw data!!!
    void setEventID(UInt_t id) { m_eventID = id; };

    ClassDef(ReconHeader,1)
};

#endif
