#ifndef reconRootData_CalClusterLayerData_H
#define reconRootData_CalClusterLayerData_H

#include <TObject.h>
#include <TVector3.h>

class CalClusterLayerData
{
public:

    CalClusterLayerData() ;
    CalClusterLayerData( Double_t energy, TVector3 position, TVector3 rmsSpread ) ;
    void init( Double_t energy, TVector3 position, TVector3 rmsSpread ) ;
    virtual ~CalClusterLayerData() ;

    const Double_t getEnergy() const { return m_energy ; }
    const TVector3 getPosition() const { return m_position ; }
    const TVector3 getRmsSpread() const { return m_rmsSpread ; }

    void Clear( Option_t * option ="" ) ;
    void Fake( Int_t ievent, UInt_t rank, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( const CalClusterLayerData &, const std::string & name = "" ) const ; // for tests
    void Print( Option_t * option ="" ) const ;

private:


    Double32_t m_energy ; // Energy deposition in crystals in this cluster and layer
    TVector3   m_position ; // Average position in this layer
    TVector3   m_rmsSpread ; // Quadratic position spread for this layer

    ClassDef(CalClusterLayerData,3)

} ;
    
#endif        






