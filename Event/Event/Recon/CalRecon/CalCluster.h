#ifndef CalCluster_H
#define CalCluster_H

#include <vector>
#include <string>
#include "geometry/Point.h"
#include "geometry/Vector.h"
#include "GaudiKernel/ObjectVector.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/MsgStream.h"
#include "Event/RelTable/RelTable.h"
#include "Event/Recon/CalRecon/CalFitParams.h"
#include "Event/Recon/CalRecon/CalParams.h"
#include "Event/Recon/CalRecon/CalXtalRecData.h"

#define NUMCALLAYERS 8

static const CLID& CLID_CalClusterCol = InterfaceID("CalClusterCol", 1, 0);

/**
*  @class CalCluster
*
*
*  @brief Defines a CalRecon TDS object which will contain the results of
*         clustering - the association of CalXtalRecData objects. This includes
*         defining a RelTable relating clusters with the CalXtalRecData
*         objects comprising them. 
* 
*         Main changes from old CalCluster.h:
*         1) Remove all "final" energy recon information (leakage, profile fitting, etc.)
*            This will be moved to the higher level object - CalEventEnergy
*         2) Incorporate layer data into a subclass and store in a vector. The idea is 
*            that the vector will always have 8 "empty" elements, presuming that external
*            usage is to loop over Cal "layers".
*         3) Provide individal "set" and "get" routines (get was there, now include set)
*  
*  @author The CalRecon Rewrite Group
*
* $Header$
*/

namespace Event 
{

/**
*   Define a subclass which will hold data by later
*
*/
class CalClusterLayerData
{
public:
    CalClusterLayerData() : m_energy(0.), m_position(0.,0.,0.), m_rmsSpread(0.,0.,0.) {}
    CalClusterLayerData(double energy, Point position, Vector rmsSpread) :
                  m_energy(energy), m_position(position), m_rmsSpread(rmsSpread) {}

    ~CalClusterLayerData() {}

    double         getEnergy()    const {return m_energy;}
    const  Point&  getPosition()  const {return m_position;}
    const  Vector& getRmsSpread() const {return m_rmsSpread;}
private:
    double m_energy;         // Energy deposition in crystals in this cluster and layer
    Point  m_position;       // Average position in this layer
    Vector m_rmsSpread;      // Quadratic position spread for this layer
};

typedef std::vector<CalClusterLayerData> CalClusterLayerDataVec;
    
class CalCluster : public CalClusterLayerDataVec, virtual public ContainedObject 
{     
public:
  
    // Define a "null" constructor
    CalCluster(int numCalLayers=NUMCALLAYERS) : CalClusterLayerDataVec(numCalLayers) {iniCluster();}
        
    virtual ~CalCluster() {}

	/// Status word bits organized like:
    /// low:   |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         < Not used   >    <Analysis Status> <   Found By   >  < Cluster Type >
    /// high:  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0  |  0   0   0   0   |
    ///         <  Not yet used                                                        >
	enum StatusBits {
    /// @param ZERO
                     ZERO            = 0x00000000, 
    /// @param ALLXTALS Bit set means all xtals in the event are here - super set
                     ALLXTALS        = 0x00000001, 
// DC: redundant with ALLXTALS
//    /// @param ISOLATEDCLUSTER Bit set means this is a subset of all Xtals
//                     ISOLATEDCLUSTER = 0x00000002,  
// DC: not used ?
//	/// @param MIPTRACK Bit set means this is a subset of all Xtals associate with MIP
//                     MIPTRACK        = 0x00000004,  
// DC: to be replaced with a string
//	/// @param SIMPLECLUSTER Bit set means cluster found using SimpleClusterTool
//                     SIMPLECLUSTER   = 0x00000010,  
//	/// @param FUZZYCLUSTER Bit set means cluster found using FuzzyClusterTool
//                     FUZZYCLUSTER    = 0x00000020,
//	/// @param MIPCLUSTER Bit set means cluster found using CalMIPFinderTool
//                     MIPCLUSTER      = 0x00000040,
// DC: CENTROID SHOULD BE REPLACED WITH AXIS
    /// @param CENTROID Bit set means centroid found
                     CENTROID        = 0x00000100, 
	/// @param MOMENTS Bit set means moments analysis run
                     MOMENTS         = 0x00000200,
// DC: not used ?
//	/// @param ENERGYCORR Bit set means one or more energy corrections were run
//	                 ENERGYCORR      = 0x00000400,
// DC: not used ?
//	/// @param MIPFIT Bit set means a MIP like track was fitted
//	                 MIPFIT          = 0x00000800
	};

    /**
    *   initializing a subset of CalCluster data members (this is
    *   primarily used by reconRootReaderAlg in translating from
    *   PDS to TDS).
    *  
    *   @param params Calorimeter Parameters for this cluster
    *   @param layerData vector of CalReconLayerData objects
    *   @param rms_long RMS of longitudinal position measurements 
    *   @param rms_trans RMS of transversal position measurements 
    */
    void initialize(const CalFitParams& fitParams,
                    const CalParams&    params,
                    double              rms_long,
                    double              rms_trans,
					double              long_asym,
                    int                 numSaturatedXtals,
					int                 numTruncXtals)
    {
        m_fitParams         = fitParams;
        m_params            = params;
        m_rmslong           = rms_long;
        m_rmstrans          = rms_trans;
		m_rmslongAsym       = long_asym;
        m_numSaturatedXtals = numSaturatedXtals;
		m_numTruncXtals     = numTruncXtals;
    }

    /*
     *   Define individual set methods here for all variables
     */
    void setFitParams(const CalFitParams& params) {m_fitParams         = params;}
    void setCalParams(const CalParams& params)    {m_params            = params;  }
    void setRmsLong(double rmsLong)               {m_rmslong           = rmsLong; }
	void setRmsLongAsym(double rmsLongAsym)       {m_rmslongAsym       = rmsLongAsym;}
    void setRmsTrans(double rmsTrans)             {m_rmstrans          = rmsTrans;}
    void setNumSaturatedXtals(int nSat)           {m_numSaturatedXtals = nSat;}
    void setNumXtals(int numXtals)                {m_numTruncXtals     = numXtals;}

    /*
     * Provide access to the data
     */
    /// Direct access to the CalFitParams
    const CalFitParams& getFitParams() const {return m_fitParams;}
    /// Direct access to the CalParams
    const CalParams& getCalParams()    const {return m_params;}
    /// get RMS of longitudinal position measurements
    double getRmsLong()		           const {return m_rmslong;} 
    /// get RMS of transverse position measurements
	double getRmsLongAsym()            const {return m_rmslongAsym;} 
    /// get RMS of transverse position measurements
    double getRmsTrans()	           const {return m_rmstrans;}
    /// get number of saturated Xtals in Cluster
    int getNumSaturatedXtals()         const {return m_numSaturatedXtals;}
	/// get Number of Truncated Xtals in Cluster
    int getNumTruncXtals()	           const {return m_numTruncXtals;}
    /// get reconstructed position
    const Point & getPosition()        const {return m_params.getCentroid();}
    /// get reconstructed direction
    const Vector & getDirection()      const {return m_params.getAxis();}

	/// Access individual status bits
    inline void setStatusBit( StatusBits bitToSet ) { m_statusBits |=  bitToSet ; }
    inline void clearStatusBit( StatusBits bitToClear ) { m_statusBits &= ~bitToClear ; }
    inline bool checkStatusBit( StatusBits bitToCheck ) const { return ((m_statusBits&bitToCheck)!=ZERO) ; }

    /// Access the status bits globally
    inline unsigned int getStatusBits() const { return m_statusBits ; }
    inline void setStatusBits( unsigned int statusBits ) { m_statusBits = statusBits ; }

    /// Access the algo name
    inline const std::string & getProducerName() const { return m_producerName ; }
    inline void setProducerName( const std::string & producerName ) { m_producerName = producerName ; }

    /// write some of CalCluster data to the ASCII output file
    /// for debugging purposes
    void writeOut(MsgStream& stream) const;
        
private:

    ///reset CalCluster data
    inline void iniCluster();
        
    //! name of the producer
    std::string m_producerName;
    //! Results of the "fit" to the cluster parameters
    CalFitParams m_fitParams;
    //! Cal Parameters
    CalParams m_params;
    //! RMS of longitudinal position measurement
    double m_rmslong;
	//! difference in RMS of longitudinal position measurement moments
    double m_rmslongAsym;
    //! RMS of transverse position measurement
    double m_rmstrans;
	//! number of "saturated" Xtals
	int m_numSaturatedXtals;
	//! number of Xtals with > 1% Total Cluster Energy
	int m_numTruncXtals;
	//! Status Bits
	unsigned int m_statusBits;

};

inline void CalCluster::iniCluster()
{
    m_fitParams         = CalFitParams();
    m_params            = CalParams();
    m_rmslong           = 0.;
	m_rmslongAsym       = 0.;
    m_rmstrans          = 0.;
	m_statusBits        = 0;
    m_numSaturatedXtals = 0;
	m_numTruncXtals     = 0; 
}

inline void CalCluster::writeOut(MsgStream& stream) const

// Purpose: provide ascii output of some data members for
//          debugging purposes
// Input:
//        stream - Gaudi message stream
{
    stream << "Energy " << m_params.getEnergy();
	stream << " No.Trunc Xtals " << m_numTruncXtals;
    stream << " " << getPosition().x() 
           << " " << getPosition().y() 
           << " " << getPosition().z();
    stream << " " << getDirection().x() 
           << " " << getDirection().y() 
           << " " << getDirection().z();
    stream << endreq;
}

//typedef for the Gaudi TDS Container
typedef ObjectVector<CalCluster>      CalClusterCol;
typedef CalClusterCol::iterator       CalClusterColItr;
typedef CalClusterCol::const_iterator CalClusterColConItr;

// Define the relational table taking us back to CalXtalRecData objects
typedef Event::RelTable<Event::CalXtalRecData, Event::CalCluster> CalClusterHitTab;
typedef Event::Relation<Event::CalXtalRecData, Event::CalCluster> CalClusterHitRel;
typedef RelationList<CalXtalRecData, CalCluster>                  CalClusterHitTabList;

}

#endif	






