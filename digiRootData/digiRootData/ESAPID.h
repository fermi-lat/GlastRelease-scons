// ESAPID.h
// Jan 1999 Daniel Flath - ROOT HTML comments added
// Jan 1999 Daniel Flath - new glast_pid info added
// Dec 05,1999 Richard Dubois Clone from CalHit

#ifndef ESAPID_H
#define ESAPID_H

#include "TObject.h"

class ESAPID : public TObject{
private:
  Bool_t	m_Clean;      // 1 if clean beam (no halo), else 0
  Float_t	m_XBeam;      // x of primary beam in mm
  Float_t	m_YBeam;      // y of primary beam in mm
  Int_t	        m_NPart;      // Number particles found by tagger
  UChar_t	m_ID[10];     // id 1=e, 2=gamma, 3=proton, 4=pion, 5=kaon
  Float_t       m_Mom[10];    // momentum each particle
  Float_t       m_BeamPb[5];  // Energy in GeV for lead glass
  Float_t       m_CADC;       
  Float_t       m_TOF[4];
  Float_t       m_TOFC;

  Bool_t        m_OneGoodP;
  Bool_t        m_OneGoodG;
  Bool_t        m_OneGoodE;
  Bool_t        m_OneGoodK;
  Bool_t        m_OneGoodPi;
  Bool_t        m_OneGoodH;
  
    void Clean();

public:
    ESAPID();
    ~ESAPID();
    inline Bool_t getClean() { return m_Clean; };
    inline Float_t getXbeam() { return m_XBeam; };
    inline Float_t getYbeam() { return m_YBeam; };
    inline Int_t getNpart() { return m_NPart; };
    inline Int_t getIDN(Int_t index) { return m_ID[index]; };
    inline Float_t getMomN(Int_t index) { return m_Mom[index]; };
    inline Float_t getBeamPbN(Int_t index) { return m_BeamPb[index]; };
    inline Float_t getCADC() { return m_CADC; };
    inline Float_t getTOFN(Int_t index) { return m_TOF[index]; };
    inline Float_t getTOFC() { return m_TOFC; };
    inline Bool_t getOneGoodP() { return m_OneGoodP; };
    inline Bool_t getOneGoodG() { return m_OneGoodG; };
    inline Bool_t getOneGoodE() { return m_OneGoodE; };
    inline Bool_t getOneGoodK() { return m_OneGoodK; };
    inline Bool_t getOneGoodPi() { return m_OneGoodPi; };
    inline Bool_t getOneGoodH() { return m_OneGoodH; };

    inline void setClean(Bool_t newVal) { m_Clean = newVal; };
    inline void setXBeam(Float_t newVal) { m_XBeam = newVal; };
    inline void setYBeam(Float_t newVal) { m_YBeam = newVal; };
    inline void setNPart(Int_t newVal) { m_NPart = newVal; };
    inline void setID(UChar_t newVal, UChar_t index = 0) { m_ID[index] = newVal; };
    inline void setMom(Float_t newVal, UChar_t index = 0) { m_Mom[index] = newVal; };
    inline void setBeamPbN(Float_t newVal, Int_t index = 0) { m_BeamPb[index] = newVal; };
    inline void setCADC(Float_t newVal) { m_CADC = newVal; };
    inline void setTOFN(Float_t newVal, Int_t index = 0) { m_TOF[index] = newVal; };
    inline void setTOFC(Float_t newVal) { m_TOFC = newVal; };
    inline void setOneGoodP(Bool_t newVal) { m_OneGoodP = newVal; };
    inline void setOneGoodG(Bool_t newVal) { m_OneGoodG = newVal; };
    inline void setOneGoodE(Bool_t newVal) { m_OneGoodE = newVal; };
    inline void setOneGoodK(Bool_t newVal) { m_OneGoodK = newVal; };
    inline void setOneGoodPi(Bool_t newVal) { m_OneGoodPi = newVal; };
    inline void setOneGoodH(Bool_t newVal) { m_OneGoodH = newVal; };
    

    ClassDef(ESAPID,2)      // ESA Particle Identification information
};

#endif
 
