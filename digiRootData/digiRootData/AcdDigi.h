#ifndef AcdDigi_H
#define AcdDigi_H

#include "TObject.h"
#include "TClonesArray.h"

#include "AcdTile.h"
#include "TileID.h"

/*! \class AcdDigi
\brief Root class to contain the ACD digitization data for
one event. 
*/
class AcdDigi: public TObject {

public:
    
    AcdDigi();

    virtual ~AcdDigi();

    void Clean(Option_t *option="");

    /// provide mechanism to update array of logs
    AcdTile* Add(TileID *id);

    /// provide access to array of logs
    TClonesArray* getTiles() { return m_tiles; };

    unsigned long eventId() { return m_eventId; };
    void eventId(unsigned int id) { m_eventId = id; };

    unsigned long timerWord() { return m_timerWord; };
    void timerWord(unsigned int time) { m_timerWord = time; };

    unsigned long status() { return m_TREQ_VETO_status; };
    void status(unsigned int stat) { m_TREQ_VETO_status = stat; };

    unsigned long deadTime() { return m_deadTime; };
    void deadTime(unsigned int time) { m_deadTime = time; };

    unsigned long deadTimeCause() { return m_deadTimeCause; };
    void deadTimeCause(unsigned int cause) { m_deadTimeCause = cause; };

private:

    TClonesArray *m_tiles;  //-> 

    static TClonesArray *m_staticAcdArray;

    unsigned int m_numTiles;

    unsigned long m_eventId;
    unsigned long m_timerWord;
    unsigned long m_TREQ_VETO_status;
    unsigned long m_deadTime;
    unsigned long m_deadTimeCause;
    
    ClassDef(AcdDigi,1)
};

#endif

