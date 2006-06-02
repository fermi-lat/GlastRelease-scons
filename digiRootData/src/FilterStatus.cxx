

#include "digiRootData/FilterStatus.h"
#include "Riostream.h"

ClassImp(LogInfo)
ClassImp(TFC_projection)
ClassImp(TFC_projectionDir)
ClassImp(TFC_projectionCol)
ClassImp(ObfTrack)
ClassImp(FilterStatus)

void LogInfo::Clear(Option_t *option) {
    m_tower = 0;
    m_layer = 0;
    m_column = 0;

    m_valN = 0;
    m_rangeN = 0;
    m_eN = 0;
    m_pedN = 0;
    m_gainN = 0;
    m_shiftN = 0;

    m_valP = 0;
    m_rangeP = 0;
    m_eP = 0;
    m_pedP = 0;
    m_gainP = 0;
    m_shiftP = 0;

}

void LogInfo::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "LogInfo: " << endl;
}



TFC_projection::TFC_projection(Int_t intercept, Int_t slope, Int_t acdTopMask, Int_t acdXMask, Int_t acdYMask, 
                               UInt_t layers, Short_t* hits, UChar_t skirtMask, UChar_t min, UChar_t max, UChar_t nhits)
{
    initialize(intercept, slope, acdTopMask, acdXMask, acdYMask, layers, hits, skirtMask, min, max, nhits);
}

void TFC_projection::initialize(Int_t intercept, Int_t slope, Int_t acdTopMask, Int_t acdXMask, Int_t acdYMask, UInt_t layers,
                                Short_t* hits, UChar_t skirtMask, UChar_t min, UChar_t max, UChar_t nhits) 
{
    m_intercept = intercept;
    m_slope = slope;
    m_acdTopMask = acdTopMask;
    m_acdXMask = acdXMask;
    m_acdYMask = acdYMask;
    m_layers = layers;
    for (iHit=0; iHit<18; iHit++) { 
        m_hits[iHit] = hits[iHit];
    }
    m_skirtMask = skirtMask;
    m_min = min;
    m_max = max;
    m_nhits = nhits;
}

void TFC_projection::Clear(Option_t *option) {
    m_intercept = 0;
    m_slope = 0;
    m_acdTopMask = '0';
    m_acdXMask = '0';
    m_acdYMask = '0';
    m_layers = 0;
    UInt_t iHit;
    for (iHit=0; iHit<18; iHit++) { 
        m_hits[iHit] = 0];
    }
    m_skirtMask = '0';
    m_min = 0;
    m_max = 0;
    m_nhits = 0;
}

void TFC_projection::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "Intercept: " << m_intercept << " Slope: " << m_slope << " layers: " << m_layers << endl;
    cout << "Min: " << m_min << " Max: " << m_max << " nHits: " << m_nhits << endl;
    cout << "acdTopMask " << m_acdTopMask << " AcdXMask " << m_acdXMask << " AcdYMask: " << m_acdYMask 
        << " SkirtMask: " << m_skirtMask << endl;
    cout << "Hits: ";
    UInt_t iHit;
    for (iHit = 0; iHit<18; iHit++) cout << m_hits[iHit] << ":";
    cout << endl;
}

TFC_projectionDir::TFC_projectionDir(UShort_t idx, UChar_t x, UChar_t y) {
    initialize(idx, x, y);
}

TFC_projectionDir::TFC_projectionDir(const TFC_projectionDir& other):TObject(other) {
    m_idx = other.m_idx;
    m_xCnt = other.m_xCnt;
    m_yCnt = other.m_yCnt;
}




void TFC_projectionDir::initialize(UShort_t idx, UChar_t x, UChar_t y) {
    m_idx = idx;
    m_xCnt = x;
    m_yCnt = y;
}

void TFC_projectionDir::Clear(Option_t *option) {
    m_idx = 0;
    m_xCnt = '0';
    m_yCnt = '0';
}

void TFC_projectionDir::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "Idx: " << m_idx << " x: " << m_x << " y: " << m_y << endl;
}

TFC_projectionCol::TFC_projectionCol(UShort_t maxCnt, UShort_t curCnt, UShort_t twrMsk, 
                                     const TFC_projectionDir *dir, const TFC_projection* projectionCol) 
{
    initialize(maxCnt, curCnt, twrMsk, dir, projectionCol);
}

TFC_projectionCol::TFC_projectionCol(const TFC_projectionCol& copy):TObject(copy) {

    initialize(copy.m_maxCnt, copy.m_curCnt, copy.m_twrMsk, copy.m_dir, copy.m_prjs);
}

void TFC_projectionCol::initialize(UShort_t maxCnt, UShort_t curCnt, UShort_t twrMsk, 
                                   const TFC_projectionDir *dir, 
                                   const TFC_projection* projectionCol)  
{
    m_maxCnt = maxCnt;
    m_curCnt = curCnt;
    m_twrMsk = twrMsk;
    UInt_t iDir;
    for (iDir = 0; iDir < 16; iDir++) {
        m_dir[iDir] = dir[iDir];
    }
    UInt_t iProj;
    for (iProj=0; iProj<m_maxCnt; iProj++) {
        m_prjs[iProj] = projectionCol[iProj];
    }
}
virtual TFC_projectionCol::~TFC_projectionCol() {

}
void TFC_projectionCol::Clear(Option_t *option) {

}

void TFC_projectionCol::Print(Option_t *option="") const {
    using namespace std;
    TObject::Print(option);
    cout << "maxCount: " << m_maxCnt << " curCnt: " << m_curCnt << " twrMsk: " << m_twrMask << endl;
}

ObjTrack::ObjTrack(Int_t numHits, Double_t phi, Double_t theta, Double_t len,
                   const vector<Double_t>& lowCoord, const vector<Double_t> highCoord,
                   const vector<Double_t>& exLowCoord, const vector<Double_t>& exHighCoord)
{
    Clear("");
    initialize(numHits, phi, theta, lowCoord, highCoord, exLowCoord, exHighCoord, len);
}

ObfTrack::ObfTrack(const ObfTrack& copy) : TObject(copy) {
    initialize(copy.m_numHits, copy.m_phi_rad, copy.m_theta_rad, copy.m_length, copy.m_lowCoord,
        copy.m_highCoord, copy.m_exLowCoord, copy.m_exHighCoord);
}

ObfTrack::~ObfTrack() {
    Clear("");
}

void ObfTrack::initialize(Int_t numHits, Double_t phi, Double_t theta, Double_t len,
                          const vector<Double_t>& lowCoord, const vector<Double_t> highCoord,
                          const vector<Double_t>& exLowCoord, const vector<Double_t>& exHighCoord) 
{
    m_numhits = numHits;
    m_phi_rad = phi;
    m_theta_rad = theta;
    m_length = len;
    m_lowCoord = lowCoord;
    m_highCoord = highCoord;
    m_exLowCoord = exLowCoord;
    m_exHighCoord = exHighCoord;

}

void ObfTrack::Clear(Option_t *option) {
    m_numhits =0;
    m_phi_rad = 0.0;
    m_theta_rad = 0.0;
    m_length = 0.0;
    m_lowCoord.clear();
    m_highCoord.clear();
    m_exLowCoord.clear();
    m_exHighCoord.clear();
}

void ObfTrack::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "numHits: " << m_numhits << " Phi: " << m_phi_rad << " Theta: " << m_theta_rad << " Len: " << m_length << endl;
}

FilterStatus::FilterStatus(const FilterStatus& copy) : TObject(copy)
{
init(copy.m_status, copy.m_stageEnergy, copy.m_tcids, copy.m_ebfSize);

initGem(copy.m_gem_thrTkr, copy.m_gem_calHiLo, copy.m_gem_condsumCno, copy.m_gem_acd_vetoes_XZ, copy.m_gem_acd_vetoes_YZ,
        copy.m_gem_acd_vetoes_XY, copy.m_gem_acd_vetoes_RU, copy.m_gem_livetime, copy.m_gem_trgtime, copy.m_gem_ppstime,
        copy.m_gem_discarded, copy.m_gem_prescaled, copy.m_gem_sent);

initAcd(copy.m_acd_xz, copy.m_acd_yz, copy.m_acd_xy, copy.m_vetoword, copy.m_acdStatus;

initTkr(copy.m_xcapture, copy.m_ycapture, copy.m_xy00, copy.m_xy11, copy.m_xy22, copy.m_xy33, copy.m_tmsk);


initSeparation(copy.m_separation);

}

virtual FilterStatus;:~FilterStatus() {

}

void FilterStatus::init(UInt_t status, Int_t stageEnergy, Int_t tcids, Int_t ebfSize) {
    m_status = status;
    m_stageEnergy = stageEnergy;
    m_tcids = tcids;
    m_ebfsize = ebfSize;
}

void FilterStatus::initGem(Int_t thrTkr, Int_t calHiLo, Int_t condsumCno, Int_t acd_vetoes_XZ, Int_t acd_vetoes_YZ, 
                           Int_t acd_vetoes_XY, Int_t acd_vetoes_RU, Int_t livetime, Int_t trgtime, Int_t ppstime,
                           Int_t discarded, Int_t prescaled, Int_t sent)
{
    m_gem_thrTkr = thrTkr;  
    m_gem_calHiLo = caHiLo;  
    m_gem_condsumCno = condsumCno;
    m_gem_acd_vetoes_XZ = acd_vetoes_XZ;
    m_gem_acd_vetoes_YZ = acd_vetoes_YZ;
    m_gem_acd_vetoes_XY = acd_veteos_XY;
    m_gem_acd_vetoes_RU = acd_vetoes_RU;
    m_gem_livetime = livetime;
    m_gem_trgtime = trgtime;   
    m_gem_ppstime = ppstime;
    m_gem_discarded = discarded; 
    m_gem_prescaled = prescaled; 
    m_gem_sent = sent;
}

void FilterStatus::initAcd(Int_t xz, Int_t yz, Int_t xy, Int_t vetoword, Int_t *status) {
    m_acd_xz = xz;
    m_acd_yz = yz;
    m_acd_xy = xy;
    m_vetoword = vetoword;
    UInt_t iStatus;
    for (iStatus = 0; iStatus < 16; iStatus++) {
        m_acdStatus[iStatus] = status[iStatus];
    }
}

void FilterStaus::initTkr(Int_t *xcapture, Int_t *ycapture, Int_t *xy00, Int_t *xy11, 
                          Int_t *xy22, Int_t *xy33, Int_t tmsk) 
{
    UInt_t iTower;
    for (iTower=0; iTower<16; iTower++) {
        m_xcapture[iTower] = xcapture[iTower];
        m_ycapture[iTower] = ycapture[iTower];
        m_xy00[iTower] = xy00[iTower];
        m_xy11[iTower] = xy11[iTower];
        m_xy22[iTower] = xy22[iTower];
        m_xy33[iTower] = xy33[iTower];
    }
    m_tmsk = tmsk;

}

void FilterStatus::initCal(Int_t numLogsHit, Float_t *layerEnergy) 
{
   m_numLogsHit = numLogsHit;
   UInt_t iLayer;
   for (iLayer=0; iLayer<8; iLayer++){
       m_layerEnergy[iLayer] = layerEnergy[iLayer];
   }


}

void FilterStatus::initLayers(Int_t *layers) 
{
    UInt_i iLayer;
    for (iLayer = 0; iLayer<16; iLayer++) {
        m_layers[iLayer] = layers[iLayer];
    }

}

void FilterStatus::initSeparation(Double_t separation) 
{
m_separation = separation;
}

///Projections for the towers
void FilterStatus::initProjectionCol(const TFC_projectionCol& projectionCol) 
{
    m_projectionCol = projectionCol;
}

void FilterStatus::initTracks(const TObjArray *tracks) 
{
//m_tracks;
}

//void FilterStatus::initLogInfo(const LogInfo *logData); 

void FilterStatus::Clear(Option_t *option) 
{
    m_status = 0;
    m_stageEnergy = 0;
    m_tcids = 0;

    m_gem_thrTkr = 0;  
    m_gem_calHiLo = 0;  
    m_gem_condsumCno = 0;
    m_gem_acd_vetoes_XZ = 0;
    m_gem_acd_vetoes_YZ = 0;
    m_gem_acd_vetoes_XY = 0;
    m_gem_acd_vetoes_RU = 0;
    m_gem_livetime = 0;
    m_gem_trgtime = 0;   
    m_gem_ppstime = 0;
    m_gem_discarded = 0; 
    m_gem_prescaled = 0; 
    m_gem_sent = 0;

    m_acd_xz = 0;
    m_acd_yz = 0;
    m_acd_xy = 0;
    m_vetoword = 0;
    UInt_t i;
    for(i = 0; i<16; i++) {
       m_acdStatus[i] = 0;

    m_xcapture[i] = 0;
    m_ycapture[i] = 0;
    m_xy00[i] = 0;
    m_xy11[i] = 0;
    m_xy22[i] = 0;
    m_xy33[i] = 0;
    m)layers[i] = 0;
    }
    
    m_tmsk = 0;

    m_ebfsize = 0;
    m_numLogsHit = 0;

    m_separation = 0.0;

    for (i=0; i<8; i++){
        m_layerEnergy[i] = 0.0;
    }

    m_projectionCol.Clear("");

    m_tracks.Delete();

}

void FilterStatus::Print(Option_t *option) const 
{
    using namespace std;
    TObject::Print(option);
}



