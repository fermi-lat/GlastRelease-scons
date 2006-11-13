#include "Event/Recon/CalRecon/GcrSelectClasses.h"

//-----------------------------------------------------------------------------------------------------------------
void Event::GcrSelectedXtal::initialize(idents::CalXtalId xtalId, float rawEnergy, float pathLength, float corrEnergy, int selectGrade, double closestFaceDist, int crossedFaces, Point entryPoint, Point exitPoint)
{
    GcrXtal::initialize(xtalId, pathLength, closestFaceDist, crossedFaces, entryPoint, exitPoint);
    
    m_rawEnergy   = rawEnergy;
    m_corrEnergy = corrEnergy;
    m_selectGrade = selectGrade;
     
}

//-----------------------------------------------------------------------------------------------------------------
void Event::GcrSelectedXtal::writeOut(MsgStream& log) const
{
    log << MSG::DEBUG;
    if (log.isActive() ) 
    {
	int it=getXtalId().getTower();
	int il=getXtalId().getLayer();
	int ic=getXtalId().getColumn();
        log << "---> writeOut GcrSelectedXtal: " << it << " " << il << " " << ic << endreq;
        log << "--->rawEnergy = " << m_rawEnergy << endreq;
        log << "---> corrEnergy = " << m_corrEnergy << endreq;
       
	
    }
}

//-----------------------------------------------------------------------------------------------------------------
 std::ostream& Event::GcrSelectedXtal::fillStream( std::ostream& s ) const 
{ 
    s << "rawEnergy ="   << m_rawEnergy   << " corrEnergy = " << m_corrEnergy << "\n";
  return s; 
}


//-----------------------------------------------------------------------------------------------------------------
void Event::GcrSelectVals::initialize(int inferedZ, int acdZ, int interactionParams)
{
    m_inferedZ   = inferedZ;
    m_acdZ = acdZ;
    m_interactionParams = interactionParams;
}
	   




