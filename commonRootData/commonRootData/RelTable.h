#ifndef ROOT_RELTABLE_H
#define ROOT_RELTABLE_H

#include "TObject.h"
#include "Relation.h"
#include "TObjArray.h"


/** @class RelTable
* @brief GLAST Relational Table class.
*
*  
*  $Header$
*/

class RelTable : public TObject  
{
public:
    
    RelTable();
    
    RelTable(const RelTable& copy);
    
    virtual ~RelTable();

    inline UInt_t getEventId() const { return m_eventId; };
    inline UInt_t getRunId() const { return m_runId; };

    void initialize(UInt_t nEvent, UInt_t nRun);
    
    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");
    
    void Print(Option_t *option="") const;
    
    void addRelation(Relation *rel);
    
    const Relation* getRelByFirst(const TObject *key) const;
    
    const TObjArray* getRelationTable() const { return m_table; };
    
    
private:
    /// unique event id for this run (for indexing)
    UInt_t m_eventId;
    
    /// Run number (for indexing)
    UInt_t m_runId;

    /// The list of relations that we really wwant
    TObjArray *m_table;
    
    ClassDef(RelTable,1) // Relational Table Class
};

#endif
