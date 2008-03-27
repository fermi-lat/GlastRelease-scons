#ifndef DigiObjectManager_H
#define DigiObjectManager_H

#include <list>

#include "digiRootData/AcdDigi.h"
#include "digiRootData/CalDigi.h"
#include "digiRootData/TkrDigi.h"
#include "digiRootData/TkrDiagnosticData.h"
#include "digiRootData/CalDiagnosticData.h"
#include "digiRootData/Tem.h"

/** 
 *  @class DigiObjectManager
 *
 *  @brief A singleton manager for handling the DIGI root objects
 *         This uses stl lists to manage a pool of objects for each root object in
 *         order to avoid using "new" and "delete". The pool can grow at a sensible 
 *         rate, it will automatically contract to a sensible size as well to try to
 *         prevent memory issues. 
 * 
 *         The objects are retrieved from the pool by overriding the "new" and "delete"
 *         operators in the specific classes. See, for example, TkrDigi.h in this package
 *
 *  ****   Currently only active for TkrDigi objects (March 27, 2008)
 *
 *  @author Iama Crazy
 *
 * $Header$
 */
class DigiObjectManager {
public:
    virtual ~DigiObjectManager() {m_pointer = 0;};

    /// The static pointer retrival method of the singleton
    static DigiObjectManager* getPointer();

    /// This method returns a "new" AcdDigi object allocated from the local pool
    AcdDigi*           getNewAcdDigi();

    /// This method returns a "new" CalDigi object allocated from the local pool
    CalDigi*           getNewCalDigi();

    /// This method returns a "new" TkrDigi object allocated from the local pool
    TkrDigi*           getNewTkrDigi();

    /// This method returns a "new" TkrDiagnosticData object allocated from the local pool
    TkrDiagnosticData* getNewTkrDiagnosticData();

    /// This method returns a "new" CalDiagnosticData object allocated from the local pool
    CalDiagnosticData* getNewCalDiagnosticData();

    /// This method returns a "new" TEM object allocated from the local pool
    Tem*               getNewTem();

    /// "Delete" all objects
    void Delete();

private:
    /// The constructor is private since this is a singleton
    DigiObjectManager();   

    /// The static pointer of the singleton
    static DigiObjectManager*              m_pointer;

    /// Define a "pool" for AcdDigi and an iterator for accessing them
    std::list<AcdDigi>                     m_acdDigiPool;
    std::list<AcdDigi>::iterator           m_acdDigiPoolIdx;

    /// Define a "pool" for CalDigi and an iterator for accessing them
    std::list<CalDigi>                     m_calDigiPool;
    std::list<CalDigi>::iterator           m_calDigiPoolIdx;

    /// Define a "pool" for TkrDigi and an iterator for accessing them
    std::list<TkrDigi>                     m_tkrDigiPool;
    std::list<TkrDigi>::iterator           m_tkrDigiPoolIdx;

    /// Define a "pool" for TkrDiagnosticData and an iterator for accessing them
    std::list<TkrDiagnosticData>           m_tkrDiagnosticDataPool;
    std::list<TkrDiagnosticData>::iterator m_tkrDiagnosticDataPoolIdx;

    /// Define a "pool" for CalDiagnosticData and an iterator for accessing them
    std::list<CalDiagnosticData>           m_calDiagnosticDataPool;
    std::list<CalDiagnosticData>::iterator m_calDiagnosticDataPoolIdx;

    /// Define a "pool" for Tem and an iterator for accessing them
    std::list<Tem>                         m_temPool;
    std::list<Tem>::iterator               m_temPoolIdx;
};
#endif //DigiObjectManager_H
