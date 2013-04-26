#ifndef ROOT_CALRECON_H
#define ROOT_CALRECON_H

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TMap.h"
#include "TString.h"

#include "CalEventEnergy.h"
#include "CalEventEnergyMap.h"
#include "CalXtalRecData.h"
#include "CalCluster.h"
#include "CalMipTrack.h"
#include "GcrXtal.h"
#include "GcrTrack.h"

/** @class CalRecon
 * @brief Primary Root object containing CAL reconstruction data.
 *
 * Portal through which one gains access to the CAL recon data. 
 * Currently there are 2 lists available:
 * - A collection of reconstructed CAL crystals (CalXtalRecData)
 * - A collection of CAL clusters (CalCluster)
 *
 * $Header$
*/

class CalRecon : public TObject
{

public:

    CalRecon();

    virtual ~CalRecon();

    void initialize();

    void Clear(Option_t *option="");
    void Fake( Int_t ievent, Float_t randNum ) ; // for tests
    Bool_t CompareInRange( CalRecon &, const std::string & name = "" ) ; // for tests
    void Print(Option_t *option="") const;

    /// provide access to list of event energies
    TObjArray* getCalEventEnergyCol() { return m_eventEnergyCol; };

    /// provide access to list of CalClusters
    TObjArray* getCalClusterCol()  { return m_clusterCol; };

    /// provide access to the list of CalXtalRecData
    TObjArray* getCalXtalRecCol()  { return m_xtalRecCol; };

    /// provide access to the list of CalMipTracks
    TObjArray* getCalMipTrackCol() { return m_mipTrackCol; };

    /// C.L. 08/22/06 provide access to the list of GcrXtals and to GcrTrack
    const TClonesArray* getGcrXtalCol() const {  return m_gcrXtalCol; };
    TObject* getGcrTrack(){return m_gcrTrack;};

    void addCalEventEnergy(CalEventEnergy* energy) { m_eventEnergyCol->Add(energy); };
    void addCalCluster(CalCluster* cluster)     { m_clusterCol->Add(cluster); };
    void addXtalRecData(CalXtalRecData *xtal)   { m_xtalRecCol->Add(xtal); };
    void addCalMipTrack(CalMipTrack* track)     { m_mipTrackCol->Add(track); };

    /// C.L. 08/22/06:
    GcrXtal* addGcrXtal(const CalXtalId& xtalId, Double_t pathLength,
                        Double_t closestFaceDist, Int_t crossedFaces, 
                        const TVector3& entry, const TVector3& exit);
    GcrXtal* addGcrXtal();
    //void addGcrXtal(GcrXtal* gcrXtal)     { m_gcrXtalCol->Add(gcrXtal); };
    void addGcrTrack(GcrTrack* gcrTrack) {m_gcrTrack = gcrTrack;};

    /// Provide access to the CalClusterMap
    TMap* getCalClusterMap() { return m_calClusterMap;}

    /// Provide access to the CalEventEnergyMap
    TMap* getCalEventEnergyMap() {return m_calEventEnergyMap;}
 
private:
    /// list of CalXtalRecData
    TObjArray* m_xtalRecCol;

    /// list of CalClusters
    TObjArray* m_clusterCol;

    /// list of CalMipTracks
    TObjArray* m_mipTrackCol;
  
    /// list of CalEventEnergy
    TObjArray* m_eventEnergyCol;

    ///CL: 08/22/06 list of gcrXtals
    TClonesArray* m_gcrXtalCol;
    TObject* m_gcrTrack;
    Int_t m_indGcrXtal;

    /// This will contain the CalClusterMap
    TMap* m_calClusterMap;

    /// This will contain the CalEventEnergyMap
    TMap* m_calEventEnergyMap;

    ClassDef(CalRecon,8)
};

#endif
