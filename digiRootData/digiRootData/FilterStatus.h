#ifndef ROOT_FILTERSTATUS_H
#define ROOT_FILTERSTATUS_H

#include "TObject.h"
#include "TObjArray.h"
#include <vector>
#ifndef R__GLOBALSTL
#ifndef WIN32
using std::vector;
#else
using namespace std;
#endif
#endif

#include "enums/EfcGammaStatus.h"


/** $Header$*/

/** @class CalLogInfo
@brief Available in v3 of OBF  
*/
class CalLogInfo : public TObject {
public:
    CalLogInfo() { Clear(""); };
    virtual ~CalLogInfo() { Clear(""); };
    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

private:
    Int_t m_tower;
    Int_t m_layer;
    Int_t m_column;

    Int_t m_valN;
    Int_t m_rangeN;
    Int_t m_eN;
    Int_t m_pedN;
    Int_t m_gainN;
    Int_t m_shiftN;

    Int_t m_valP;
    Int_t m_rangeP;
    Int_t m_eP;
    Int_t m_pedP;
    Int_t m_gainP;
    Int_t m_shiftP;

    ClassDef(CalLogInfo, 1) // CalLogInfo
};


/** @class TfcProjection
@brief   
*/
class TfcProjection : public TObject {
public:
    TfcProjection() { };

    TfcProjection(Int_t intercept, Int_t slope, Int_t acdTopMask, 
                   Int_t acdXMask, Int_t acdYMask, UInt_t layers, 
                   Short_t* hits, UChar_t skirtMask, UChar_t min, UChar_t max, 
                   UChar_t nhits);

    TfcProjection(const TfcProjection& copy);

    /// define assignment operator
    TfcProjection& operator=(const TfcProjection& copy); 

    virtual ~TfcProjection() {};
    void initialize(Int_t intercept, Int_t slope, Int_t acdTopMask, 
                    Int_t acdXMask, Int_t acdYMask, UInt_t layers,
                    const Short_t* hits, UChar_t skirtMask, UChar_t min, 
                    UChar_t max, UChar_t nhits);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;


    inline Int_t getIntercept() const { return m_intercept; }
    inline Int_t getSlope() const { return m_slope; }
    inline Int_t getAcdTopMask() const { return m_acdTopMask; }
    inline Int_t getAcdXMask() const { return m_acdXMask; }
    inline Int_t getAcdYMask() const { return m_acdYMask; }
    inline UInt_t getLayers() const { return m_layers; }
    inline const Short_t* getHits() const { return m_hits; }
    inline UChar_t getSkirtMask() const { return m_skirtMask; }
    inline UChar_t getMin() const { return m_min; }
    inline UChar_t getMax() const { return m_max; }
    inline UChar_t getNhits() const { return m_nhits; }

private:
    /// intercept at beginning layer
    Int_t  m_intercept;
    Int_t  m_slope;
    /// ACD top tile candidates
    Int_t  m_acdTopMask; 
    /// ACD x facing candidates
    Int_t  m_acdXMask; 
    /// ACD y facing candidates
    Int_t  m_acdYMask; 
    /// Bit mask representing the struck layers
    UInt_t m_layers;
    /// Hits assigned to projection
    Short_t m_hits[18]; 
    /// Mask of which skirt region projection strikes
    UChar_t m_skirtMask; 
    /// Beginning layer number of the projection
    UChar_t m_min; 
    /// Ending layer number of the projection
    UChar_t m_max; 
    /// Number of hits assigned
    UChar_t m_nhits;

    ClassDef(TfcProjection,1) // TfcProjection
};

/** @class TfcProjectionDir
@brief   
*/
class TfcProjectionDir: public TObject {
public:
    TfcProjectionDir() { };
    TfcProjectionDir(UShort_t idx, UChar_t x, UChar_t y);
    TfcProjectionDir(const TfcProjectionDir& copy);
    /// define assignment operator
    TfcProjectionDir& operator=(const TfcProjectionDir& copy); 

    virtual ~TfcProjectionDir(){ };

    void initialize(UShort_t idx, UChar_t x, UChar_t y);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    inline UShort_t getIndex() const { return m_idx; }
    inline UChar_t getXCount() const { return m_xCnt; }
    inline UChar_t getYCount() const { return m_yCnt; }

private:
    /// Index to the first projection
    UShort_t  m_idx; 
    /// Count of the X projections
    UChar_t   m_xCnt; 
    /// Count of Y projections
    UChar_t   m_yCnt; 
    ClassDef(TfcProjectionDir,1) // TfcProjectionDir
};


/** @class TfcProjectionCol
@brief   
*/
class TfcProjectionCol : public TObject {
public:
    TfcProjectionCol():TObject() { m_maxCnt = 0; m_prjs = 0; Clear(""); };
    
    TfcProjectionCol(Int_t maxCnt, UShort_t curCnt, UShort_t twrMas, 
        const TfcProjectionDir *dir, const TfcProjection* projectionCol);

    TfcProjectionCol(const TfcProjectionCol& other);

    /// define assignment operator
    TfcProjectionCol& operator=(const TfcProjectionCol& copy); 

    virtual ~TfcProjectionCol();

    void initialize(Int_t maxCnt, UShort_t curCnt, UShort_t twrMas, 
        const TfcProjectionDir *dir, const TfcProjection* projectionCol);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    void Fake( Int_t ievent, Float_t randNum ); 

    inline Int_t getMaxCount() const { return m_maxCnt; };
    inline UShort_t getCurrentCount() const { return m_curCnt; };
    inline UShort_t getTwrMsk() const { return m_twrMsk; };
    inline const TfcProjectionDir* getProjectionDirCol() const { return m_dir; }
    inline const TfcProjection* getProjectionCol() const { return m_prjs; }

private:
    ///Maximum # of projections available
    Int_t    m_maxCnt; 
    //Current # of projections in use   
    UShort_t    m_curCnt;
    ///Mask of the tower with projections  
    UShort_t    m_twrMsk; 
    /// Directory of the projections by twr 
    TfcProjectionDir    m_dir[16];
    TfcProjection* m_prjs; //[m_curCnt]

    ClassDef(TfcProjectionCol,2) // TfcProjectionCol
};

class ObfTrack : public TObject {
public:
    ObfTrack(){ Clear(""); };
    ObfTrack(Int_t numHits, Double_t phi, Double_t theta, Double_t len,
        const vector<Double_t>& lowCoord, const vector<Double_t> highCoord,
        const vector<Double_t>& exLowCoord, const vector<Double_t>& exHighCoord);

    ObfTrack(const ObfTrack &copy);
    virtual ~ObfTrack();

    void initialize(Int_t numHits, Double_t phi, Double_t theta, Double_t len,
        const vector<Double_t>& lowCoord, const vector<Double_t>& highCoord,
        const vector<Double_t>& exLowCoord, const vector<Double_t>& exHighCoord);

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    inline Int_t getNumHits() const { return m_numHits; }
    inline Double_t getPhi() const { return m_phi_rad; }
    inline Double_t getTheta() const { return m_theta_rad; }
    inline const vector<Double_t>& getLowCoord() const { return m_lowCoord; }
    inline const vector<Double_t>& getHighCoord() const { return m_highCoord; }
    inline const vector<Double_t>& getExLowCoord() const { return m_exLowCoord; }
    inline const vector<Double_t>& getExHighCoord() const { return m_exHighCoord; }
    inline Double_t getLength() const { return m_length;};

private:
    Int_t m_numHits;
    Double_t m_phi_rad;
    Double_t m_theta_rad;
    Double_t m_length;
    vector<Double_t> m_lowCoord;
    vector<Double_t> m_highCoord;
    vector<Double_t> m_exLowCoord;
    vector<Double_t> m_exHighCoord;

    ClassDef(ObfTrack,1) // ObfTrack
};

/** @class FilterStatus
* @brief Store the Onboard Filter Status
*
*/

class FilterStatus : public TObject {
public:

    FilterStatus(){ Clear(""); };
    FilterStatus(const FilterStatus& copy);
    virtual ~FilterStatus();


    /// define assignment operator
    FilterStatus& operator=(const FilterStatus& copy); 

    /// For tests
    void Fake( Int_t ievent, Float_t randNum ); 
    Bool_t CompareInRange( FilterStatus &, const std::string & name = "" );


    void init(UInt_t status, Int_t stageEnergy, Int_t tcids, Int_t ebfSize);

    void initGem(Int_t thrTkr, Int_t calHiLo, Int_t condsumCno, 
                 Int_t acd_vetoes_XZ, Int_t acd_vetoes_YZ, Int_t acd_vetoes_XY,
                 Int_t acd_vetoes_RU, Int_t livetime, Int_t trgtime, 
                 Int_t ppstime, Int_t discarded, Int_t prescaled, Int_t sent);

    void initAcd(Int_t xz, Int_t yz, Int_t xy, Int_t vetoword, 
                 const Int_t *status);

    void initTkr(const Int_t *xcapture, const Int_t *ycapture, 
                 const Int_t *xy00, const Int_t *xy11, const Int_t *xy22, 
                 const Int_t *xy33, Int_t tmsk);

    void initCal(Int_t numLogsHit, const Float_t *layerEnergy);

    void initLayers(const Int_t *layers);

    void initSeparation(Double_t separation);

    ///Projections for the towers
    void initProjectionCol(const TfcProjectionCol& projectionCol);

    void initTracks(const TObjArray &tracks);

    void initBestTrack( Int_t xHits, Int_t yHits,
                              Double_t slopeXZ,Double_t slopeYZ,
                              Double_t intXZ,Double_t intYZ);
       

    //void initCalLogInfo(const CalLogInfo *logData); 

    void Clear(Option_t *option ="");

    void Print(Option_t *option="") const;

    ///Return the statuscode of the filter
    inline UInt_t getStatus() const { return m_status; };
    ///Returns the 16 most significant bits of the filter code
    unsigned int getHigh() const { return m_status>>15; };
    ///Returns the 16 least signifcat bits of the filter code
    unsigned int getLow() const { return m_status & 0x7FFF; };

    ///Returns the value stored in CalEnergy
    inline Float_t getCalEnergy() const { 
        return (Float_t)((m_stageEnergy & enums::EFC_GAMMA_STAGE_M_ENERGY)/4.0);  }
    inline Int_t getStageEnergy() const{ return m_stageEnergy;  }
    ///Return the Code specifying the towers with triggers or possible triggers
    inline Int_t getTcids()const{  return m_tcids; }
    inline Int_t getGemThrTkr()const { return m_gem_thrTkr; };
    inline Int_t getGemCalHiLo() const { return m_gem_calHiLo; };
    inline Int_t getGemCondsumCno() const { return m_gem_condsumCno;} ;

    inline Int_t getGemAcd_vetoes_XZ()const{ return m_gem_acd_vetoes_XZ;  }
    inline Int_t getGemAcd_vetoes_YZ()const{ return m_gem_acd_vetoes_YZ;  }
    inline Int_t getGemAcd_vetoes_XY()const{   return m_gem_acd_vetoes_XY;  }
    inline Int_t getGemAcd_vetoes_RU()const{   return m_gem_acd_vetoes_RU;  }
    inline Int_t getGemLivetime()const{    return m_gem_livetime;  }
    inline Int_t getGemTrgtime()const{ return m_gem_trgtime;  }  
    inline Int_t getGemPpstime()const{    return m_gem_ppstime;  }
    inline Int_t getGemDiscarded()const{ return m_gem_discarded;  } 
    inline Int_t getGemPrescaled()const{    return m_gem_prescaled;  } 
    inline Int_t getGemSent()const{    return m_gem_sent;  }

    ///Return the ACD hit map results
    inline void getAcdMap(int &xz, int &yz, int &xy) const {
        xz=m_acd_xz;
        yz=m_acd_yz;
        xy=m_acd_xy;
    }
    ///Return the ACD faces intersected by projections
    inline const Int_t* getAcdStatus() const {
        return m_acdStatus;
    }

    ///Return pointer to array of layers that were hit in each tower
    inline const Int_t* getLayerCol()const { return m_layers;  }

    inline const TfcProjectionCol* getTfcProjectionCol(){  return &m_projectionCol;  }

    ///Return all available tracks
    inline const TObjArray& getTrackCol()const { return m_tracks;   }

    ///Return angular separation between best track and incomming particle
    inline Double_t getSeparation() const{    return m_separation;  }

    inline Bool_t tracksExist() const{
        if(m_tracks.GetEntriesFast()!=0)
            return true;
        return false;
    }

    inline const Int_t* getXcapture() const{    return m_xcapture;  }

    inline const Int_t* getYcapture() const{    return m_ycapture;  }

    ///return the Cal Layer Energy
    inline const Float_t* getLayerEnergy() const{    return m_layerEnergy;  }
    ///Return xy00
    inline const Int_t* getXY00() const{    return m_xy00;  }
    ///Return xy11
    inline const Int_t* getXY11() const{    return m_xy11;  }
    ///Return xy22
    inline const Int_t* getXY22() const{    return m_xy22;  }
    ///Return xy33
    inline const Int_t* getXY33() const{    return m_xy33;  }


    ///Return tmsk
    inline int getTmsk() const{    return m_tmsk;  }

    inline int getVetoWord() const{    return m_vetoword;  }

    inline void getBestTrack(Int_t &xHits,Int_t &yHits,
        Double_t &slopeXZ,Double_t &slopeYZ,
        Double_t &intXZ,Double_t &intYZ) const {
            xHits    = m_xHits;
            yHits    = m_yHits;
            slopeXZ  = m_slopeXZ;
            slopeYZ  = m_slopeYZ;
            intXZ    = m_intXZ;
            intYZ    = m_intYZ;
        }

        inline Int_t getXhitsBestTrack() const { return m_xHits; }
        inline Int_t getYhitsBestTrack() const { return m_yHits; }
        inline Double_t getSlopeXzBestTrack() const { return m_slopeXZ; }
        inline Double_t getSlopeYzBestTrack() const { return m_slopeYZ; }
        inline Double_t getIntXzBestTrack() const { return m_intXZ; }
        inline Double_t getIntYzBestTrack() const { return m_intYZ; }

private:

    ///Filter status code
    UInt_t m_status;
    ///Energy in CAL
    Int_t m_stageEnergy;
    ///Towers with triggers
    Int_t m_tcids;

    ///Gem info
    Int_t m_gem_thrTkr;  
    Int_t m_gem_calHiLo;  
    Int_t m_gem_condsumCno;
    Int_t m_gem_acd_vetoes_XZ;
    Int_t m_gem_acd_vetoes_YZ;
    Int_t m_gem_acd_vetoes_XY;
    Int_t m_gem_acd_vetoes_RU;
    Int_t m_gem_livetime;
    Int_t m_gem_trgtime;   
    Int_t m_gem_ppstime;
    Int_t m_gem_discarded; 
    Int_t m_gem_prescaled; 
    Int_t m_gem_sent;

    ///ACD hit map
    Int_t m_acd_xz;
    Int_t m_acd_yz;
    Int_t m_acd_xy;
    Int_t m_vetoword;
    ///ACD faces
    Int_t m_acdStatus[16];

    Int_t m_xcapture[16];
    Int_t m_ycapture[16];
    Int_t m_xy00[16];
    Int_t m_xy11[16];
    Int_t m_xy22[16];
    Int_t m_xy33[16];
    Int_t m_tmsk;

    Int_t m_ebfSize;  /// In new OBF
    Int_t m_numLogsHit; /// In new OBF

    ///Layers hit in each tower
    Int_t m_layers[16];
    
    /// Best Track data
    Int_t m_xHits;   
    Int_t m_yHits;   
    Double_t m_slopeXZ;
    Double_t m_slopeYZ;
    Double_t m_intXZ;
    Double_t m_intYZ;

    ///Angular separation between best track and incomming particle
    Double_t m_separation;

    Float_t m_layerEnergy[8];

    ///Projections for the towers
    TfcProjectionCol m_projectionCol;  //Disable for right now

    ///Tracks found for this event
    TObjArray m_tracks;

    //CalLogInfo m_calLogData[16*8*12];    // 16 towers * 8 layers * 12 logs  in newer version of OBF

    ClassDef(FilterStatus,2) // OBF Filter Status
};

#endif
