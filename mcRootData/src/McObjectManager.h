#ifndef McObjectManager_H
#define McObjectManager_H

#include <list>

#include "mcRootData/McParticle.h"
#include "mcRootData/McPositionHit.h"
#include "mcRootData/McIntegratingHit.h"
#include "mcRootData/McTrajectory.h"

/** 
 *  @class McObjectManager
 *
 *  @brief A singleton manager for handling the MC root objects
 *
 *
 *  @author Am. Crazy
 *
 * $Header$
 */
class McObjectManager {
public:
    virtual ~McObjectManager(){ m_pointer=0; };

    /// The static pointer retrival method of the singleton
    static McObjectManager* getPointer();

    /// This method returns a "new" McParticle allocated from the local pool
    McParticle* getNewMcParticle();

    /// This method returns a "new" McPositionHit allocated from the local pool
    McPositionHit* getNewMcPositionHit();

    /// This method returns a "new" McIntegratingHit allocated from the local pool
    McIntegratingHit* getNewMcIntegratingHit();

    /// This method returns a "new" McTrajectory allocated from the local pool
    McTrajectory* getNewMcTrajectory();

    /// This method returns a "new" McTrajectoryPoint allocated from the local pool
    McTrajectoryPoint* getNewMcTrajectoryPoint();

    /// "Delete" all objects
    void Delete();

private:
    /// The constructor is private since this is a singleton
    McObjectManager();   

    /// The static pointer of the singleton
    static McObjectManager*                m_pointer;

    /// Define a "pool" for McParticles and an iterator for accessing them
    std::list<McParticle>                  m_mcPartPool;
    std::list<McParticle>::iterator        m_mcPartPoolIdx;

    /// Define a "pool" for McParticles and an iterator for accessing them
    std::list<McPositionHit>               m_mcPosHitPool;
    std::list<McPositionHit>::iterator     m_mcPosHitPoolIdx;

    /// Define a "pool" for McParticles and an iterator for accessing them
    std::list<McIntegratingHit>            m_mcIntHitPool;
    std::list<McIntegratingHit>::iterator  m_mcIntHitPoolIdx;

    /// Define a "pool" for McTrajectory's and an iterator for accessing them
    std::list<McTrajectory>                m_mcTrajectoryPool;
    std::list<McTrajectory>::iterator      m_mcTrajectoryPoolIdx;

    /// Define a "pool" for McTrajectoryPoint's and an iterator for accessing them
    std::list<McTrajectoryPoint>           m_mcTrajPointPool;
    std::list<McTrajectoryPoint>::iterator m_mcTrajPointPoolIdx;
};
#endif //McObjectManager_H
