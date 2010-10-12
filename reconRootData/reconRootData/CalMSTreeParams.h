
#ifndef reconRootData_CalMSTreeParams_H
#define reconRootData_CalMSTreeParams_H

#include <TObject.h>

/** 
* @class CalMSTreeParams
* @brief Root clone for the TDS CalMSTreeParams
* @authors Luca Baldini, Johan Bregeon, Carmelo Sgro'
*
* $Header: 
*/

class CalMSTreeParams
{
public:

    /// Default constructor
    CalMSTreeParams() { Clear() ; }

    /// Direct construction from all the elements
    CalMSTreeParams
    ( Double32_t totalEnergy,
      Double32_t maxXtalEnergy,  Int_t    numberOfEdges,		      
      Double32_t minEdgeLength,  Double32_t maxEdgeLength,		      
      Double32_t meanEdgeLength, Double32_t meanEdgeLengthTrunc,	      
      Double32_t rmsEdgeLength,  Double32_t rmsEdgeLengthTrunc);	      

    /// init
    void init
    ( Double32_t totalEnergy,
      Double32_t maxXtalEnergy,  Int_t    numberOfEdges,		      
      Double32_t minEdgeLength,  Double32_t maxEdgeLength,		      
      Double32_t meanEdgeLength, Double32_t meanEdgeLengthTrunc,	      
      Double32_t rmsEdgeLength,  Double32_t rmsEdgeLengthTrunc);

    /// Destructor
    virtual ~CalMSTreeParams() {}

    ///--------------------------------------------------- 
    ///--- Get methods
    /// Retrieve the total enegy
    inline Double32_t  getTotalEnergy()      const {return m_totalEnergy;}
    /// Retrieve the energy in the crystal with the maximum enegy
    inline Double32_t  getMaxXtalEnergy()    const {return m_maxXtalEnergy;}
    /// Retrieve the number of edges
    inline Int_t     getNumberOfEdges()    const {return m_numberOfEdges;}
    /// Retrieve the minimum edge length
    inline Double32_t  getMinEdgeLength()         const {return m_minEdgeLength;}
    /// Retrieve the minimum edge length
    inline Double32_t  getMaxEdgeLength()         const {return m_maxEdgeLength;}
    /// Retrieve the maximum edge length
    inline Double32_t  getMeanEdgeLength()        const {return m_meanEdgeLength;}
    /// Retrieve the average edge length
    inline Double32_t  getMeanEdgeLengthTrunc()   const {return m_meanEdgeLengthTrunc;}
    /// Retrieve the RMS of edges length
    inline Double32_t  getRmsEdgeLength()         const {return m_rmsEdgeLength;}
    /// Retrieve the RMS of edges length  after truncation
    inline Double32_t  getRmsEdgeLengthTrunc()    const {return m_rmsEdgeLengthTrunc;}

    /// Retrieve the number of cristals from the number of edges
    inline Int_t getNumXtals()          const {return m_numberOfEdges + 1;}

    /// Other methods
    void Print( Option_t * option ="" ) const ;
    void Clear( Option_t * option ="" ) ;
    /// for tests
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; 
    Bool_t CompareInRange( const CalMSTreeParams &, const std::string & name = "" ) const ;

private:

    /// Sum of the energy of the crystal within the Tree
    /// Total energy is also present in the CalParams container: m_energy
    /// Pan is to keep this as a cross check value here
    Double32_t m_totalEnergy;
    /// Energy in the crystal with the maximum energy
    Double32_t m_maxXtalEnergy;
    /// Number of edges
    Int_t m_numberOfEdges;

    /// More quantities from test_MST.py
    /// Not implemented yet in MSTTree
    Double32_t m_minEdgeLength;
    Double32_t m_maxEdgeLength;
    Double32_t m_meanEdgeLength;
    Double32_t m_meanEdgeLengthTrunc;
    Double32_t m_rmsEdgeLength;
    Double32_t m_rmsEdgeLengthTrunc;

    ClassDef(CalMSTreeParams,1)
} ;




#endif

