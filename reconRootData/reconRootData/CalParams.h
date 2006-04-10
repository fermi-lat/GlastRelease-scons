
#ifndef reconRootData_CalParams_H
#define reconRootData_CalParams_H

#include <TObject.h>
#include <TVector3.h>

/** 
* @class CalParams
*
* @brief Root clone for the TDS CalParams
* @author Bill Atwood, Tracy Usher, David Chamont
*
* $Header$
*/

class CalParams
{
public:

    /// Default constructor
    CalParams() { Clear() ; }

    /// Direct construction from all the elements (the old fashioned way)
    CalParams
    ( Double_t energy, Double_t eneError,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz ) ;

    /// init
    void init
    ( Double_t energy, Double_t eneError,
      TVector3 centroid,
      Double_t cntdxx, Double_t cntdxy, Double_t cntdxz,
      Double_t cntdyy, Double_t cntdyz, Double_t cntdzz,
      TVector3 axis,
      Double_t axsdxx, Double_t axsdxy, Double_t axsdxz,
      Double_t axsdyy, Double_t axsdyz, Double_t axsdzz ) ;

    virtual ~CalParams() {}

    /// Retrieve the energy
    inline const Double_t getEnergy() const {return m_energy;}
    inline const Double_t getEnergyErr() const {return m_eneError;}

    /// Retrieve the centroid
    inline const TVector3  getCentroid()     const {return m_clusterCentroid;}

    /// Retrieve the centroid errors
    inline const Double_t getxPosxPos()     const {return m_cenxx; }
    inline const Double_t getxPosyPos()     const {return m_cenxy; }
    inline const Double_t getxPoszPos()     const {return m_cenxz; }
    inline const Double_t getyPosyPos()     const {return m_cenyy; }
    inline const Double_t getyPoszPos()     const {return m_cenyz; }
    inline const Double_t getzPoszPos()     const {return m_cenzz; }

    /// Retrieve the axis
    inline const TVector3 getAxis()         const { return m_clusterAxis ; }

    /// Retrieve the axis errors
    inline const Double_t getxDirxDir()     const {return m_axisxx; }
    inline const Double_t getxDiryDir()     const {return m_axisxy; }
    inline const Double_t getxDirzDir()     const {return m_axisxz; }
    inline const Double_t getyDiryDir()     const {return m_axisyy; }
    inline const Double_t getyDirzDir()     const {return m_axisyz; }
    inline const Double_t getzDirzDir()     const {return m_axiszz; }

    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalParams &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:

    /// Energy and associated error
    Double32_t m_energy;
    Double32_t m_eneError;

    /// Cluster centroid
    TVector3  m_clusterCentroid;

    /// Error matrix for centroid stored in upper diagonal form
    Double32_t m_cenxx;     // Cov(1,1) = dx * dx 
    Double32_t m_cenxy;     // Cov(1,2) = Cov(2,1) = dx * dy
    Double32_t m_cenxz;     // Cov(1,3) = Cov(3,1) = dx * dz
    Double32_t m_cenyy;     // Cov(2,2) = dy * dy
    Double32_t m_cenyz;     // Cov(2,3) = Cov (3,2)= dy * dz
    Double32_t m_cenzz;     // Cov(3,3) = dz * dz

    /// Cluster axis, assumed to pass through centroid
    TVector3 m_clusterAxis;

    /// Error matrix for centroid stored in upper diagonal form
    Double32_t m_axisxx;    // Cov(1,1) = dx * dx 
    Double32_t m_axisxy;    // Cov(1,2) = Cov(2,1) = dx * dy
    Double32_t m_axisxz;    // Cov(1,3) = Cov(3,1) = dx * dz
    Double32_t m_axisyy;    // Cov(2,2) = dy * dy
    Double32_t m_axisyz;    // Cov(2,3) = Cov (3,2)= dy * dz
    Double32_t m_axiszz;    // Cov(3,3) = dz * dz
    
    ClassDef(CalParams,2)
} ;




#endif

