

#include "digiRootData/FilterStatus.h"
#include "commonRootData/RootDataUtil.h"
#include "Riostream.h"

ClassImp(LogInfo)
ClassImp(TfcProjection)
ClassImp(TfcProjectionDir)
ClassImp(TfcProjectionCol)
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



TfcProjection::TfcProjection(Int_t intercept, Int_t slope, Int_t acdTopMask, 
                               Int_t acdXMask, Int_t acdYMask, UInt_t layers, 
                               Short_t* hits, UChar_t skirtMask, UChar_t min, 
                               UChar_t max, UChar_t nhits)
{
    initialize(intercept, slope, acdTopMask, acdXMask, acdYMask, layers, hits, 
               skirtMask, min, max, nhits);
}

TfcProjection::TfcProjection(const TfcProjection& other):TObject(other) {

    initialize(other.m_intercept, other.m_slope, other.m_acdTopMask, 
               other.m_acdXMask, other.m_acdYMask, other.m_layers, 
               other.m_hits, other.m_skirtMask, other.m_min, other.m_max,  
               other.m_nhits);
}

TfcProjection& TfcProjection::operator=(const TfcProjection& other) {
    initialize(other.m_intercept, other.m_slope, other.m_acdTopMask, 
               other.m_acdXMask, other.m_acdYMask, other.m_layers, 
               other.m_hits, other.m_skirtMask, other.m_min, other.m_max,  
               other.m_nhits);
    return *this;
}

void TfcProjection::initialize(Int_t intercept, Int_t slope, Int_t acdTopMask,
                                Int_t acdXMask, Int_t acdYMask, UInt_t layers,
                                const Short_t* hits, UChar_t skirtMask, 
                                UChar_t min, UChar_t max, UChar_t nhits) 
{
    m_intercept = intercept;
    m_slope = slope;
    m_acdTopMask = acdTopMask;
    m_acdXMask = acdXMask;
    m_acdYMask = acdYMask;
    m_layers = layers;
    UInt_t iHit;
    for (iHit=0; iHit<18; iHit++) { 
        m_hits[iHit] = hits[iHit];
    }
    m_skirtMask = skirtMask;
    m_min = min;
    m_max = max;
    m_nhits = nhits;
}

void TfcProjection::Clear(Option_t *option) {
    m_intercept = 0;
    m_slope = 0;
    m_acdTopMask = '0';
    m_acdXMask = '0';
    m_acdYMask = '0';
    m_layers = 0;
    UInt_t iHit;
    for (iHit=0; iHit<18; iHit++) { 
        m_hits[iHit] = 0;
    }
    m_skirtMask = '0';
    m_min = 0;
    m_max = 0;
    m_nhits = 0;
}

void TfcProjection::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "Intercept: " << m_intercept << " Slope: " << m_slope 
         << " layers: " << m_layers << endl;
    cout << "Min: " << m_min << " Max: " << m_max << " nHits: " << m_nhits 
         << endl;
    cout << "acdTopMask " << m_acdTopMask << " AcdXMask " << m_acdXMask 
         << " AcdYMask: " << m_acdYMask 
         << " SkirtMask: " << m_skirtMask << endl;
    cout << "Hits: ";
    UInt_t iHit;
    for (iHit = 0; iHit<18; iHit++) cout << m_hits[iHit] << ":";
    cout << endl;
}

TfcProjectionDir::TfcProjectionDir(UShort_t idx, UChar_t x, UChar_t y) {
    initialize(idx, x, y);
}

TfcProjectionDir::TfcProjectionDir(const TfcProjectionDir& other):TObject(other) {
    initialize(other.m_idx, other.m_xCnt, other.m_yCnt);
}


TfcProjectionDir& TfcProjectionDir::operator=(const TfcProjectionDir& other) {
    initialize(other.m_idx, other.m_xCnt, other.m_yCnt);
    return *this;
}


void TfcProjectionDir::initialize(UShort_t idx, UChar_t x, UChar_t y) {
    m_idx = idx;
    m_xCnt = x;
    m_yCnt = y;
}

void TfcProjectionDir::Clear(Option_t *option) {
    m_idx = 0;
    m_xCnt = '0';
    m_yCnt = '0';
}

void TfcProjectionDir::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "Idx: " << m_idx << " x: " << m_xCnt << " y: " << m_yCnt << endl;
}

TfcProjectionCol::TfcProjectionCol(Int_t maxCnt, UShort_t curCnt, 
                                     UShort_t twrMsk, 
                                     const TfcProjectionDir *dir, 
                                     const TfcProjection* projectionCol) 
{
    initialize(maxCnt, curCnt, twrMsk, dir, projectionCol);
}

TfcProjectionCol::TfcProjectionCol(const TfcProjectionCol& copy):TObject(copy) {

    initialize(copy.m_maxCnt, copy.m_curCnt, copy.m_twrMsk, copy.m_dir, 
               copy.m_prjs);
}


TfcProjectionCol& TfcProjectionCol::operator=(const TfcProjectionCol& copy) {
    initialize(copy.m_maxCnt, copy.m_curCnt, copy.m_twrMsk, copy.m_dir, 
               copy.m_prjs);
    return *this;
}

void TfcProjectionCol::initialize(Int_t maxCnt, UShort_t curCnt, 
                                   UShort_t twrMsk, 
                                   const TfcProjectionDir *dir, 
                                   const TfcProjection* projectionCol)  
{
    m_maxCnt = maxCnt;
    m_curCnt = curCnt;
    m_twrMsk = twrMsk;
    UInt_t iDir;
    for (iDir = 0; iDir < 16; iDir++) {
        m_dir[iDir] = dir[iDir];
    }
    if (projectionCol) {
    m_prjs = new TfcProjection[m_maxCnt];
    UInt_t iProj;
    for (iProj=0; iProj<m_maxCnt; iProj++) {
        std::cout << "creating TfcProjection  " << iProj << std::endl;
        if (!projectionCol) std::cout << "input col is null" << std::endl;
        m_prjs[iProj] = projectionCol[iProj];
    }
    } else {
        m_prjs = 0;
        m_maxCnt = 0;
    }
}

TfcProjectionCol::~TfcProjectionCol() {
    Clear("");

}
void TfcProjectionCol::Clear(Option_t *option) {

    if (m_maxCnt && m_prjs) {
        delete [] m_prjs;
    }
    m_maxCnt = 0;
    m_prjs = 0;
    m_curCnt = 0;
    m_twrMsk = 0;
    UInt_t iDir;
    for (iDir=0; iDir < 16; iDir++) {
        m_dir[iDir].Clear("");
    }

}

void TfcProjectionCol::Print(Option_t *option) const {
    using namespace std;
    TObject::Print(option);
    cout << "maxCount: " << m_maxCnt << " curCnt: " << m_curCnt << " twrMsk: " 
         << m_twrMsk << endl;
    Int_t i;
    for (i=0; i<m_maxCnt; i++) {
        cout << "TfcProjection: " << i << std::endl;
        m_prjs[i].Print();
    }
}

void TfcProjectionCol::Fake( Int_t ievent, Float_t randNum ) {

    std::cout << "Called TfcProjectionColFake" << std::endl;
    m_curCnt = 5;
    m_twrMsk = 4;
    UInt_t iDir;
    for (iDir=0; iDir < 16; iDir++) {
        m_dir[iDir].initialize(6, '1', '2');
    }
    m_maxCnt = 2;
    m_prjs = new TfcProjection[m_maxCnt];
    UInt_t iCnt;
    const Short_t hits[18] = {10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27};
    for (iCnt = 0; iCnt < m_maxCnt; iCnt++) {

        m_prjs[iCnt].initialize(1, 2, 3, 4, 5, 6, hits, '8', '9', '0', '1');

    }
}

ObfTrack::ObfTrack(Int_t numHits, Double_t phi, Double_t theta, Double_t len,
                   const vector<Double_t>& lowCoord, 
                   const vector<Double_t> highCoord,
                   const vector<Double_t>& exLowCoord,  
                   const vector<Double_t>& exHighCoord)
{
    Clear("");
    initialize(numHits, phi, theta, len, lowCoord, highCoord, exLowCoord, 
               exHighCoord);
}

ObfTrack::ObfTrack(const ObfTrack& copy) : TObject(copy) {
    initialize(copy.m_numHits, copy.m_phi_rad, copy.m_theta_rad, copy.m_length,
               copy.m_lowCoord, copy.m_highCoord, copy.m_exLowCoord, 
               copy.m_exHighCoord);
}

ObfTrack::~ObfTrack() {
    Clear("");
}

void ObfTrack::initialize(Int_t numHits, Double_t phi, Double_t theta, 
                          Double_t len, const vector<Double_t>& lowCoord, 
                          const vector<Double_t>& highCoord,
                          const vector<Double_t>& exLowCoord, 
                          const vector<Double_t>& exHighCoord) 
{
    m_numHits = numHits;
    m_phi_rad = phi;
    m_theta_rad = theta;
    m_length = len;
    m_lowCoord = lowCoord;
    m_highCoord = highCoord;
    m_exLowCoord = exLowCoord;
    m_exHighCoord = exHighCoord;

}

void ObfTrack::Clear(Option_t *option) {
    m_numHits =0;
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
    cout << "numHits: " << m_numHits << " Phi: " << m_phi_rad << " Theta: " 
         << m_theta_rad << " Len: " << m_length << endl;
}

FilterStatus::FilterStatus(const FilterStatus& copy) : TObject(copy)
{
    init(copy.m_status, copy.m_stageEnergy, copy.m_tcids, copy.m_ebfSize);

    initTracks(copy.m_tracks);
    initLayers(copy.m_layers);

    initGem(copy.m_gem_thrTkr, copy.m_gem_calHiLo, copy.m_gem_condsumCno, 
            copy.m_gem_acd_vetoes_XZ, copy.m_gem_acd_vetoes_YZ,
            copy.m_gem_acd_vetoes_XY, copy.m_gem_acd_vetoes_RU, 
            copy.m_gem_livetime, copy.m_gem_trgtime, copy.m_gem_ppstime,
            copy.m_gem_discarded, copy.m_gem_prescaled, copy.m_gem_sent);

    initAcd(copy.m_acd_xz, copy.m_acd_yz, copy.m_acd_xy, copy.m_vetoword, 
            copy.m_acdStatus);

    initTkr(copy.m_xcapture, copy.m_ycapture, copy.m_xy00, copy.m_xy11, 
            copy.m_xy22, copy.m_xy33, copy.m_tmsk);

    initCal(copy.m_numLogsHit, copy.m_layerEnergy);

    initSeparation(copy.m_separation);

    initProjectionCol(copy.m_projectionCol);
}

FilterStatus::~FilterStatus() {

}



FilterStatus& FilterStatus::operator=(const FilterStatus& copy) {

    init(copy.m_status, copy.m_stageEnergy, copy.m_tcids, copy.m_ebfSize);

    initTracks(copy.m_tracks);

    initLayers(copy.m_layers);

    initGem(copy.m_gem_thrTkr, copy.m_gem_calHiLo, copy.m_gem_condsumCno, 
            copy.m_gem_acd_vetoes_XZ, copy.m_gem_acd_vetoes_YZ,
            copy.m_gem_acd_vetoes_XY, copy.m_gem_acd_vetoes_RU, 
            copy.m_gem_livetime, copy.m_gem_trgtime, copy.m_gem_ppstime,
            copy.m_gem_discarded, copy.m_gem_prescaled, copy.m_gem_sent);

    initAcd(copy.m_acd_xz, copy.m_acd_yz, copy.m_acd_xy, copy.m_vetoword, 
            copy.m_acdStatus);

    initTkr(copy.m_xcapture, copy.m_ycapture, copy.m_xy00, copy.m_xy11, 
            copy.m_xy22, copy.m_xy33, copy.m_tmsk);


    initCal(copy.m_numLogsHit, copy.m_layerEnergy);

    initSeparation(copy.m_separation);

    initProjectionCol(copy.m_projectionCol);

    return *this;
  }

void FilterStatus::init(UInt_t status, Int_t stageEnergy, Int_t tcids, 
                        Int_t ebfSize) {
    m_status = status;
    m_stageEnergy = stageEnergy;
    m_tcids = tcids;
    m_ebfSize = ebfSize;
}

void FilterStatus::initGem(Int_t thrTkr, Int_t calHiLo, Int_t condsumCno, 
                           Int_t acd_vetoes_XZ, Int_t acd_vetoes_YZ, 
                           Int_t acd_vetoes_XY, Int_t acd_vetoes_RU, 
                           Int_t livetime, Int_t trgtime, Int_t ppstime,
                           Int_t discarded, Int_t prescaled, Int_t sent)
{
    m_gem_thrTkr = thrTkr;  
    m_gem_calHiLo = calHiLo;  
    m_gem_condsumCno = condsumCno;
    m_gem_acd_vetoes_XZ = acd_vetoes_XZ;
    m_gem_acd_vetoes_YZ = acd_vetoes_YZ;
    m_gem_acd_vetoes_XY = acd_vetoes_XY;
    m_gem_acd_vetoes_RU = acd_vetoes_RU;
    m_gem_livetime = livetime;
    m_gem_trgtime = trgtime;   
    m_gem_ppstime = ppstime;
    m_gem_discarded = discarded; 
    m_gem_prescaled = prescaled; 
    m_gem_sent = sent;
}

void FilterStatus::initAcd(Int_t xz, Int_t yz, Int_t xy, Int_t vetoword, 
                           const Int_t *status) {
    m_acd_xz = xz;
    m_acd_yz = yz;
    m_acd_xy = xy;
    m_vetoword = vetoword;
    UInt_t iStatus;
    for (iStatus = 0; iStatus < 16; iStatus++) {
        m_acdStatus[iStatus] = status[iStatus];
    }
}

void FilterStatus::initTkr(const Int_t *xcapture, const Int_t *ycapture, 
                          const Int_t *xy00, const Int_t *xy11, 
                          const Int_t *xy22, const Int_t *xy33, Int_t tmsk) 
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

void FilterStatus::initCal(Int_t numLogsHit, const Float_t *layerEnergy) 
{
   m_numLogsHit = numLogsHit;
   UInt_t iLayer;
   for (iLayer=0; iLayer<8; iLayer++){
       m_layerEnergy[iLayer] = layerEnergy[iLayer];
   }


}

void FilterStatus::initLayers(const Int_t *layers) 
{
    UInt_t iLayer;
    for (iLayer = 0; iLayer<16; iLayer++) {
        m_layers[iLayer] = layers[iLayer];
    }

}

void FilterStatus::initSeparation(Double_t separation) 
{
m_separation = separation;
}

///Projections for the towers
void FilterStatus::initProjectionCol(const TfcProjectionCol& projectionCol) 
{
    m_projectionCol = projectionCol;
} 

void FilterStatus::initTracks(const TObjArray &tracks) 
{
    m_tracks.Delete();
    Int_t i;
    for (i = 0; i<tracks.GetEntriesFast(); i++) {
        m_tracks.Add(tracks.At(i));
    }
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
       m_layers[i] = 0;
    }
    
    m_tmsk = 0;

    m_ebfSize = 0;
    m_numLogsHit = 0;

    m_separation = 0.0;

    for (i=0; i<8; i++){
        m_layerEnergy[i] = 0.0;
    }

    //m_projectionCol.Clear("");

    m_tracks.Delete();

}

void FilterStatus::Print(Option_t *option) const 
{
    using namespace std;
    TObject::Print(option);

    cout << "Status: " << m_status << " stageEnergy: " << m_stageEnergy
         << " Tower Triggers: " << m_tcids << endl;

    cout << "GEM: " << endl;
    cout << "ThrTkr: " << m_gem_thrTkr << " calHiLo: " << m_gem_calHiLo
         << " CondsumCno: " << m_gem_condsumCno << endl;
    cout << "ACD Vetos:  XZ: " << m_gem_acd_vetoes_XZ << " YZ: "
         << m_gem_acd_vetoes_YZ << " XY: " << m_gem_acd_vetoes_XY 
         << " RU: " << m_gem_acd_vetoes_RU << endl;
    cout << "LiveTime: " << m_gem_livetime << " TrgTime: " << m_gem_trgtime
         << " PPSTime: " << m_gem_ppstime << endl;
    cout << "Discarded: " << m_gem_discarded << " Prescaled: " 
         << m_gem_prescaled << " Sent: " << m_gem_sent << endl;

    cout << "ACD Hit Map: " << endl;
    cout << "XZ: " << m_acd_xz << " YZ: " << m_acd_yz << " XY: " << m_acd_xy
         << " vetoword: " << m_vetoword << endl;
    UInt_t i;
    for (i=0; i<16; i++) 
        cout << m_acdStatus[i] << ":";
    cout << endl;

    cout << "Num\t" << "xCapture \t" << "yCapture\t" << "xy00\t" << "xy11\t" 
         << "xy22\t" << "xy33\t" << endl;
    for (i=0; i<16;i++) {
        cout << i << "\t" << m_xcapture[i] << "\t" << m_ycapture[i] << "\t"
             << m_xy00[i] << "\t" << m_xy11[i] << "\t" << m_xy22[i] << "\t"
             << m_xy33[i] << endl;
    }
    cout << "TMSK: " << m_tmsk << endl;

    cout << "ebfSize: " << m_ebfSize << " numLogsHit: " << m_numLogsHit << endl;

    cout << "Layers: " << endl;
    for (i=0; i<16; i++)
        cout << m_layers[i] << ":";
    cout << endl;

    cout << "separation: " << m_separation << endl;

    cout << "LayerEnergy:" <<endl;
    for (i = 0; i < 8; i++) 
        cout << m_layerEnergy[i] << ":";
    cout << endl;

    TIter trackIt(&m_tracks);
    ObfTrack *curTrack;
    while ((curTrack = ((ObfTrack*)trackIt.Next()))) {
        curTrack->Print(option);
    }

    m_projectionCol.Print(option);
}

//======================================================
// For Unit Tests
//======================================================

void FilterStatus::Fake( Int_t ievent, Float_t randNum ) {

  init(1,2,3,0);
  initGem(4,5,6,7,8,9,10,11,12,13,14,15,16);
  const Int_t stat[16] = {30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};
  initAcd(17,18,19,20,stat);
  Int_t xCapture[16] = {70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85};
  Int_t yCapture[16] = {90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105};
  Int_t xy00[16]={40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
  Int_t xy11[16]={40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55};
  Int_t xy22[16]={20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,55};
  Int_t xy33[16]={10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
  initTkr(xCapture,yCapture,xy00,xy11,xy22,xy33,20);
  Int_t layers[16] = {50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65};
  initLayers(layers);
  Float_t layerEnergy[8] = {1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8};
  initCal(0,layerEnergy);
  initSeparation(10.0);

  TObjArray trackObf;
  initTracks(trackObf);
  TfcProjectionCol proj;
  proj.Fake(ievent, randNum);
  initProjectionCol(proj);
    
}

#define COMPARE_TOBJ_ARRAY_IN_RANGE(T,m) rootdatautil::TObjArrayCompareInRange<T>(m,ref.m)

Bool_t FilterStatus::CompareInRange( FilterStatus & ref, const std::string & name ) {

    bool result = true ;
    result = rootdatautil::CompareInRange(getStatus(),ref.getStatus(),"Status") && result;
    result = rootdatautil::CompareInRange(getHigh(),ref.getHigh(),"High") && result;
    result = rootdatautil::CompareInRange(getLow(),ref.getLow(),"Low") && result;
    result = rootdatautil::CompareInRange(getTcids(),ref.getTcids(),"Tcids") && result;
    result = rootdatautil::CompareInRange(getCalEnergy(),ref.getCalEnergy(),"CalEnergy") && result;
    result = rootdatautil::CompareInRange(getGemThrTkr(),ref.getGemThrTkr(),"GemThrTkr") && result;
    result = rootdatautil::CompareInRange(getGemCalHiLo(),ref.getGemCalHiLo(),"GemCalHiLo") && result;
    result = rootdatautil::CompareInRange(getGemCondsumCno(),ref.getGemCondsumCno(),"GemCondsumCno") && result;
    result = rootdatautil::CompareInRange(getGemAcd_vetoes_XZ(),ref.getGemAcd_vetoes_XZ(),"GemAcd_vetoes_XZ") && result;
    result = rootdatautil::CompareInRange(getGemAcd_vetoes_YZ(),ref.getGemAcd_vetoes_YZ(),"GemAcd_vetoes_YZ") && result;
    result = rootdatautil::CompareInRange(getGemAcd_vetoes_XY(),ref.getGemAcd_vetoes_XY(),"GemAcd_vetoes_XY") && result;
    result = rootdatautil::CompareInRange(getGemAcd_vetoes_RU(),ref.getGemAcd_vetoes_RU(),"GemAcd_vetoes_RU") && result;
    result = rootdatautil::CompareInRange(getGemLivetime(),ref.getGemLivetime(),"GemLivetime") && result;
    result = rootdatautil::CompareInRange(getGemTrgtime(),ref.getGemTrgtime(),"GemTrgtime") && result;
    result = rootdatautil::CompareInRange(getGemPpstime(),ref.getGemPpstime(),"GemPpstime") && result;
    result = rootdatautil::CompareInRange(getGemDiscarded(),ref.getGemDiscarded(),"GemDiscarded") && result;
    result = rootdatautil::CompareInRange(getGemPrescaled(),ref.getGemPrescaled(),"GemPrescaled") && result;
    result = rootdatautil::CompareInRange(getGemSent(),ref.getGemSent(),"GemSent") && result;

    if (!result) {
        if (name == "") {
            std::cout<<"Comparison ERROR for " << ClassName() << std::endl;
        } else {
            std::cout << "Comparison ERROR for " << name<<std::endl;
        }
    }

    return result ;

}



