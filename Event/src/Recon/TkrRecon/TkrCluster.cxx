#include "Event/Recon/TkrRecon/TkrCluster.h"
#include "idents/TowerId.h"

//---------------------------------------------------
//       TkrCluster
//---------------------------------------------------

using namespace Event;

void TkrCluster::writeOut(MsgStream& log) const

{
    // Purpose: writes out debug info
	// Inputs:  message stream
	// Outputs: data written to message stream
	
    log << MSG::DEBUG;
    if (log.isActive()) {
      log << " plane " << layer() << " XY " << (int) v()
            << " xpos  " << m_position.x()  << " ypos   " << m_position.y()
            << " zpos  " << m_position.z()
            << " i0-if " << m_strip0 <<"-"<< m_stripf;
    }
    log << endreq;
}
//---------  Private --------------------------------


int TkrCluster::chip()     
{ 
	// Purpose: returns chip number of this cluster
	// Method: calculated from the strip number
	return m_strip0/64;
}

double TkrCluster::strip() 
{ 
	// Purpose: return the center strip number
	// Method: average of first and last strip
	// Caveat: may be half-integer
	return 0.5*(m_strip0+m_stripf);
}

double TkrCluster::size()  
{ 
	// Purpose: returns number of strips in the cluster
    return std::abs(m_stripf-m_strip0) + 1.;
}

TkrCluster::view TkrCluster::intToView(int iv)

{
	// Purpose: converts a view number to an enum
	
	TkrCluster::view v = XY;
	if (iv == 0) v = X;
	else if (iv == 1) v =Y;
	return v;
}


int TkrCluster::viewToInt(TkrCluster::view v)

{
	// Purpose: converts an enum to an integer
	
	if (v == TkrCluster::XY) return 2;
	return (v == TkrCluster::X? 0:1);
}

int TkrCluster::tower() const
{
  idents::TowerId tId(m_tkrId.getTowerX(), m_tkrId.getTowerY());
  return tId.id();
}
