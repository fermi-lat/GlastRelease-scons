#ifndef CalMSTreeParams_H
#define CalMSTreeParams_H

/** 
* @class CalMSTreeParams
*
* @brief Gaudi TDS class to store the parameters of the cluster Minimum Spanning Tree
* 
* @authors Luca Baldini, Johan Bregeon, Carmelo Sgro'
*
* $Header:
*/

#include <iostream>

namespace Event { //Namespace Event


class CalMSTreeParams
{
public:
    /// Default constructor
    CalMSTreeParams() { clear() ; }

    ///// Constructor given a minimum spanning tree
    //    Forget about it for now, as it would require
    //    to connect the package with CalRecon
    //CalMSTreeParams(const MSTTree& mstree);

    /// Direct construction from all the elements
    CalMSTreeParams(double totalEnergy,
                    double maxXtalEnergy,  int    numberOfEdges,
                    double minEdgeLength,  double maxEdgeLength,
                    double meanEdgeLength, double meanEdgeLengthTrunc,
                    double rmsEdgeLength,  double rmsEdgeLengthTrunc);

    /// Destructor
    ~CalMSTreeParams() {}

    /// Reset method
    void clear() ;
  
    ///--------------------------------------------------- 
    ///--- Get methods
    /// Retrieve the total enegy
    inline double  getTotalEnergy()      const {return m_totalEnergy;}
    /// Retrieve the energy in the crystal with the maximum enegy
    inline double  getMaxXtalEnergy()    const {return m_maxXtalEnergy;}
    /// Retrieve the number of edges
    inline int     getNumberOfEdges()    const {return m_numberOfEdges;}
    /// Retrieve the minimum edge length
    inline double  getMinEdgeLength()         const {return m_minEdgeLength;}
    /// Retrieve the minimum edge length
    inline double  getMaxEdgeLength()         const {return m_maxEdgeLength;}
    /// Retrieve the maximum edge length
    inline double  getMeanEdgeLength()        const {return m_meanEdgeLength;}
    /// Retrieve the average edge length
    inline double  getMeanEdgeLengthTrunc()   const {return m_meanEdgeLengthTrunc;}
    /// Retrieve the RMS of edges length
    inline double  getRmsEdgeLength()         const {return m_rmsEdgeLength;}
    /// Retrieve the RMS of edges length  after truncation
    inline double  getRmsEdgeLengthTrunc()    const {return m_rmsEdgeLengthTrunc;}

    /// Retrieve the number of cristals from the number of edges
    inline int getNumXtals()          const {return m_numberOfEdges + 1;}


    ///--------------------------------------------------- 
    ///--- Set methods
    /// Set the total energy
    inline void    setTotalEnergy(   double totEnergy     )    {m_totalEnergy = totEnergy;}
    /// Set the energy in the crystal with the maximum enegy
    inline void    setMaxXtalEnergy( double maxXtalEnergy )    {m_maxXtalEnergy = maxXtalEnergy;}
    /// Set the number of edges
    inline void    setNumberOfEdges( int numEdges      )       {m_numberOfEdges = numEdges;}
    /// Set the minimum edge length
    inline void  setMinEdgeLength(   double length )           {m_minEdgeLength = length;}
    /// Set the minimum edge length
    inline void  setMaxEdgeLength(   double length )           {m_maxEdgeLength = length;}
    /// Set the maximum edge length
    inline void  setMeanEdgeLength(  double length )           {m_meanEdgeLength = length;}
    /// Set the average edge length
    inline void  setMeanEdgeLengthTrunc( double length )       {m_meanEdgeLengthTrunc = length;}
    /// Set the RMS of edges length
    inline void  setRmsEdgeLength(   double length )           {m_rmsEdgeLength = length;}
    /// Set the RMS of edges length  after truncation
    inline void  setRmsEdgeLengthTrunc( double length )        {m_rmsEdgeLengthTrunc = length;}

    ///--------------------------------------------------- 
    ///--- Other methods
    /// Dump the content of the container: the MST tree parameters
    /// How does this work, it's also implemented in the .cxx ? -- Johan
    std::ostream& fillStream( std::ostream& s ) const;
    friend std::ostream& operator<< ( std::ostream& s, const CalMSTreeParams& obj ) 
    {
	    return obj.fillStream(s);
    }

private:

    /// Sum of the energy of the crystal within the Tree
    /// Total energy is also present in the CalParams container: m_energy
    /// Pan is to keep this as a cross check value here
    double m_totalEnergy;
    /// Energy in the crystal with the maximum energy
    double m_maxXtalEnergy;
    /// Number of edges
    int m_numberOfEdges;

    /// More quantities
    double m_minEdgeLength;
    double m_maxEdgeLength;
    double m_meanEdgeLength;
    double m_meanEdgeLengthTrunc;
    double m_rmsEdgeLength;
    double m_rmsEdgeLengthTrunc;
    
};


}; //Namespace Event


#endif

