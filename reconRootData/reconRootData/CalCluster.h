#ifndef CALCLUSTER_H
#define CALCLUSTER_H

// Root Includes
#include "TObject.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TArrayF.h"
//! Root object that stores a cluster of CsI detections
/*! Provides access to the positions and energy depositions for a
    cluster of CsI logs */
class CalCluster : public TObject
{

private:
    //! sum of logs
    Float_t  m_Esum;

    //! corrected energy
    Float_t  m_Ecorr;

    //! energy corrected for leakage ( valid at high energies)
    Float_t m_Eleak;

    // energy and position for each layer
    //! energy per layer
    TArrayF m_eneLayer;
    //! x pos for each layer
    TArrayF m_xposLayer;
    //! y position for each layer
    TArrayF m_yposLayer;

    //  profile fitting parameters
    //! alpha parameter used for the fit
    Float_t  m_CsiAlpha;
    //! lambda parameter used for the fit
    Float_t  m_CsiLambda;
    
    //! Fitted energy
    Float_t  m_fitEnergy;
    
    //! Chi square of profile fitting
    Float_t  m_ProfChisq;

    //! starting point of the shower
    Float_t  m_start;

    Int_t m_id;                 // index in Object array?

    //! energy-weighted center of gravity of the cluster
    TVector3  *m_position;   
    //! direction of the cluster  (Px/P, Py/P, Pz/P)
    TVector3  *m_direction;     
  
    //! list of CalLogEne used in this cluster
    TObjArray *m_calLogs;       
	
public:

    //! class constructor
        CalCluster (Int_t index=0);
    
    //! destructor
    virtual ~CalCluster();

    //! frees dynamically allocated memory
    void Clean();

    //! allocates memory as necessary
    void Create();

    //! Access energy sum
    void setEsum (Float_t ene) { m_Esum = ene; };
    Float_t getEsum () { return m_Esum; };
  
    //! Access corrected energy
    void setEcorr (Float_t ene) { m_Ecorr = ene; };
    Float_t getEcorr () { return m_Ecorr; };

    //! Access leaking energy
    void setEleak (Float_t ene) { m_Eleak = ene; };
    Float_t getEleak () { return m_Eleak; };



    //! Access energy-weighted center of the cluster
    TVector3* getPosition() { return m_position; };
    inline void setPosition(TVector3 *p) { m_position = p; };
    inline void setPosition(Double_t x, Double_t y, Double_t z) { m_position->SetXYZ(x,y,z); };

    //! Access the direction
    TVector3* getDirection() { return m_direction; };
    inline void setDirection(TVector3 *d) { m_direction = d; };
    inline void setDirection(Double_t x, Double_t y, Double_t z) { m_direction->SetXYZ(x,y,z); };

	//!  Access to profile fitting parameters
	Float_t getCsiAlpha() { return m_CsiAlpha;}
	Float_t getCsiLambda() { return m_CsiLambda;}
	Float_t getCsiStart(){ return m_start;}
	Float_t getProfChisq(){ return m_ProfChisq;}
	Float_t getFitEnergy(){ return m_fitEnergy;}

	void setCsiAlpha(Float_t a){m_CsiAlpha = a;}
	void setCsiLambda(Float_t l){m_CsiLambda = l;}
	void setCsiStart(Float_t s){m_start = s;}
	void setFitEnergy(Float_t e){m_fitEnergy = e;}
	void setProfChisq(Float_t c){m_ProfChisq = c;}

	// Access to energy and position in each layer

        //! Access all X positions within a layer
	const TArrayF& getXLayer() const {return m_xposLayer;}
        //! Access a specific X position within a layer
	Float_t getXLayer(int i) {return m_xposLayer[i];}
        //! Access all Y positions within a layer
	const TArrayF& getYLayer() const {return m_xposLayer;}
        //! Access a specific Y position within a layer
	Float_t getYLayer(int i) {return m_yposLayer[i];}

        //! Access all energy depositions within a layer
	const TArrayF& getEneLayer() const {return m_eneLayer;}
        //! Access a specific energy within a layer
	Float_t getEneLayer(int i) {return m_eneLayer[i];}

	void setEneLayer(int i,Float_t e) {m_eneLayer[i]=e;}
	void setXLayer (int i,Float_t x) {m_xposLayer[i]=x;}
	void setYLayer (int i,Float_t y) {m_yposLayer[i]=y;}

    //! Access the list of logs associated with this cluster
    TObjArray* getCalLog () { return m_calLogs; };
    
    //! Access the ID of this cluster
    Int_t getId() { return m_id; };
    //! Set the ID of this cluster
    void setId(Int_t id) { m_id = id; };

    ClassDef(CalCluster,1)
};


#endif 


