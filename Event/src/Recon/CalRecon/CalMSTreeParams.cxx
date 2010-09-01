// File and Version information:
// $Header:
//
//  Implementation file of CalMSTreeParams container to be filled with
//  quantities from minimum spanning tree MSTTree
//  
// Authors:
//
//        Luca Baldini, Johan Bregeon, Carmelo Sgro'
//
//

#include "Event/Recon/CalRecon/CalMSTreeParams.h"


// // Constructor from the MSTTree object
// // Forget about it for now, as it would require
// // to link the package with CalRecon
// Event::CalMSTreeParams::CalMSTreeParams(const MSTTree& mstree):
// {
//     setNumberOfEdges(mstree.size())
// }

// Constructor by giving out all the arguments
Event::CalMSTreeParams::CalMSTreeParams(double totalEnergy,
                                        double maxXtalEnergy,  int    numEdges,
                                        double minEdgeLength,  double maxEdgeLength,
                                        double meanEdgeLength, double meanEdgeLengthTrunc,
                                        double rmsEdgeLength,  double rmsEdgeLengthTrunc):
       m_totalEnergy(totalEnergy),
       m_maxXtalEnergy(maxXtalEnergy),   m_numberOfEdges(numEdges),
       m_minEdgeLength(minEdgeLength),   m_maxEdgeLength(maxEdgeLength),
       m_meanEdgeLength(meanEdgeLength), m_meanEdgeLengthTrunc(meanEdgeLengthTrunc),
       m_rmsEdgeLength(rmsEdgeLength),   m_rmsEdgeLengthTrunc(rmsEdgeLengthTrunc)
{
return;
}

// Clear up the container from any sensible value
void Event::CalMSTreeParams::clear()
{
    m_totalEnergy    = 0.;
    m_maxXtalEnergy  = 0.;
    m_numberOfEdges  = 0;
    m_minEdgeLength  = 0.;
    m_maxEdgeLength  = 0.;
    m_meanEdgeLength = 0.;
    m_meanEdgeLengthTrunc = 0.;
    m_rmsEdgeLength  = 0.;
    m_rmsEdgeLengthTrunc  = 0.;
}

// Dump the content of the container: the MST tree parameters
std::ostream& Event::CalMSTreeParams::fillStream( std::ostream& s ) const
{
    s << m_totalEnergy	       << " " << m_maxXtalEnergy << " " << m_numberOfEdges << " "
      << m_minEdgeLength       << " " << m_maxEdgeLength << " " << m_meanEdgeLength<< " "
      << m_meanEdgeLengthTrunc << " " << m_rmsEdgeLength << " " << m_rmsEdgeLengthTrunc;
      
    return s; 
}

