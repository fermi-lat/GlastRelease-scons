// $Header$

#include <iostream>
#include <math.h>
#include "Event/MonteCarlo/McTrajectory.h"
#include "Event/Utilities/CLHEPStreams.h"


namespace Event{

/// Retrieve pointer to McParticle (const or non-const)
//const McParticle* McTrajectory::getMcParticle() const
//{
//    return m_mcParticle; 
//}
//McParticle* McTrajectory::getMcParticle()
//{
//    return m_mcParticle; 
//}

McTrajectory::~McTrajectory()
{
    for(std::vector<McTrajectoryPoint*>::iterator trajIter = m_points.begin(); trajIter != m_points.end(); trajIter++)
    {
        delete *trajIter;
    }
}

/// Update pointer to McParticle (by a C++ pointer or a smart reference)
//void McTrajectory::setMcParticle( McParticle* value )
//{
//    m_mcParticle = value; 
//}
//void McTrajectory::setMcParticle( SmartRef<McParticle> value )
//{ 
//    m_mcParticle = value; 
//}


void McTrajectory::addPoints(std::vector<Event::McTrajectoryPoint*>& points)
{
  m_points = points;
}

void McTrajectory::addPoint(Event::McTrajectoryPoint* point)
{
    m_points.push_back(point);
}

}
