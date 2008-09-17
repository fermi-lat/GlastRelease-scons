#ifndef GCRSELECTEVENT_H
#define GCRSELECTEVENT_H

#include "TObject.h"
#include "GcrSelect.h"

/** @class GcrSelectEvent
 * @brief Primary Root Recon Object - top level object in Root Tree.
 *
 * A GcrSelectEvent object contains:
 * -# Run Number
 * -# Event Id
 * -# A pointer to a GcrSelect object
 * 	
 * $Header$
*/

/** @page gcrSelectRootData_doc
    @section GcrSelectEvent GcrSelectEvent Documentation
    \image html gcrSelectRootDataTree.gif
    

<table>
<tr>
<th>Method</th><th>Description</th><th>Return Value</th>
</tr>
<small>
<tr><td>getEventId( )</td><td>Returns the Event Id</td><td>UInt_t</td>
</tr><tr>
<td>getRunId( )</td><td>Returns the Run Id</td><td>UInt_t</td>
</tr><tr>
<td>getGcrSelect( )</td><td>Returns a pointer to the GcrSelect object</td><td>GcrSelect*</td>
</tr><tr>

</small>
</table>

*/



class GcrSelectEvent : public TObject
{

public:

    /// Define the bits in the m_eventFlag member
    typedef enum {
        GOOD = 0,
        EVTSEQ = 1
    } EventFlags;

    GcrSelectEvent();

    virtual ~GcrSelectEvent();

    void initialize(UInt_t eventId, UInt_t runId, GcrSelect *gcrSelect);

    void initEventFlags(UInt_t flags) { m_eventFlags = flags; };

    void Clear(Option_t *option="");
    void Fake( Int_t ievent, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const GcrSelectEvent &, const std::string & name = "" ) const ; // for tests
    void Print(Option_t *option="") const;

    UInt_t getEventId()  const{ return m_eventId; };

    UInt_t getRunId() const { return m_runId; };

    const GcrSelect* getGcrSelect()  const {return m_gcrSelect;}
    GcrSelect* getGcrSelect() {return m_gcrSelect;}
    


    UInt_t getEventFlags() const { return m_eventFlags; };
    Bool_t goodEvent() const { return (m_eventFlags == 0); };
    Bool_t badEvent() const { return (m_eventFlags != 0); };
    Bool_t badEventSeq() const { return (m_eventFlags && EVTSEQ); };

private:
	/// Event Id
    UInt_t m_eventId;
	/// Run number
    UInt_t m_runId;
    /// Do the easy thing and add flags here
    UInt_t m_eventFlags;
    
    GcrSelect *m_gcrSelect;

    ClassDef(GcrSelectEvent,4) 
};

#endif

