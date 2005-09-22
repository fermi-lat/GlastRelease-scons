#ifndef MCTRAJECTORY_H
#define MCTRAJECTORY_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"

#include "McParticle.h"

/** @class McTrajectory
* @brief GLAST Monte Carlo Particle Trajectory class.
*
* The McTrajectory class contains the step by step position information for
* particles tracked through GLAST during the simulation stage. The class contains:
* - A reference to the McParticle to which it belongs
* - A TObjArray of positions as the particle traversed the GLAST universe
* - The charge of the particle 
* The MC data stored to ROOT files mirrors the data available within the Gaudi 
* TDS.
*  
*  $Header$
*/

class McTrajectory : public TObject  
{
public: 
    McTrajectory();
    
    virtual ~McTrajectory();
    
    /// clear lists, free pointers, etc., after read from / write to file
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;
    
    void initialize(McParticle* mcPart, UInt_t charge);

    /// get/set McParticle associated with this trajectory
    const McParticle* getMcParticle() const {return (McParticle*)m_mcParticle.GetObject();}
    void setMcParticle(McParticle* mcPart) {m_mcParticle = mcPart;}

    /// get/set charge of trajectory
    UInt_t getMcTrajectoryCharge() const {return m_charge;}
    void setMcTrajectoryCharge(UInt_t charge) {m_charge = charge;}

    /// return the full TObjArray containing the points on the trajectory
    TObjArray* getMcPointCol() const { return m_pointCol; };

    /// retrieve number of points stored in the collection
    inline UInt_t getMcPointCount() const 
    { 
        return ((m_pointCol) ? m_pointCol->GetEntries() : 0);
    };
    
    /// add a TVector3 to points collection
    void addMcPoint(TVector3* point) {m_pointCol->Add(point);}

    /// clear of the array (necessary for converters)
    void clearMcPointCol()  { m_pointCol->Clear(); };

private:
    /// Pointer to the McParticle to which this trajectory is associated
    TRef       m_mcParticle;

    /// An array of points defining the trajectory
    TObjArray* m_pointCol;

    /// Charge of the trajectory (mirrors TDS version)
    UInt_t     m_charge;
    
    ClassDef(McTrajectory,1) // Monte Carlo Trajectory Class
};

#endif
