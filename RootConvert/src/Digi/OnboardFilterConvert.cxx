#include "RootConvert/Digi/OnboardFilterConvert.h"

#include <RootConvert/Utilities/Toolkit.h>


namespace RootPersistence {

void convert(const TFC_prj &tdsObj, TfcPrj &rootObj) 
{
    TfcPrjPrms rootTopPrms,rootBotPrms;
    convert(tdsObj.top, rootTopPrms);
    convert(tdsObj.bot, rootBotPrms);

    TfcHit rootHit[18];

    for(int idx=0; idx<tdsObj.nhits; idx++)
    {
        convert(tdsObj.hits[idx], rootHit[idx]);
    }

    rootObj.initialize(rootTopPrms, rootBotPrms, tdsObj.acdTopMask,
                       tdsObj.acdXMask, tdsObj.acdYMask, tdsObj.skirtMask, tdsObj.min, tdsObj.max, 
                       tdsObj.nhits, tdsObj.layers, rootHit);  
}

void convert(const TfcPrj &rootObj, TFC_prj &tdsObj) 
{
    convert(rootObj.getTop(), tdsObj.top);
    convert(rootObj.getBot(), tdsObj.bot);

    tdsObj.acdTopMask = rootObj.getAcdTopMask();
    tdsObj.acdXMask   = rootObj.getAcdXMask();
    tdsObj.acdYMask   = rootObj.getAcdYMask();
    tdsObj.skirtMask  = rootObj.getSkirtMask();
    tdsObj.min        = rootObj.getMin();
    tdsObj.max        = rootObj.getMax();
    tdsObj.nhits      = rootObj.getNhits();
    tdsObj.layers     = rootObj.getLayers(); 

    const TfcHit* hitsRoot = rootObj.getHits();
    unsigned int i;
    for (i=0;i<18;i++) 
        convert(hitsRoot[i], tdsObj.hits[i]); 
}

void convert(const TFC_hit &tdsObj, TfcHit& rootObj) 
{
    rootObj.initialize(tdsObj.strip, tdsObj.tower, tdsObj.idx);
}

void convert(const TfcHit &rootObj, TFC_hit& tdsObj) 
{
    tdsObj.strip = rootObj.getStrip();
    tdsObj.tower = rootObj.getTower();
    tdsObj.idx   = rootObj.getIdx();
}

void convert(const TFC_prjPrms &tdsObj, TfcPrjPrms& rootObj) 
{
    rootObj.initialize(tdsObj.xy0,tdsObj.dxy,tdsObj.dz,tdsObj.tower,tdsObj.layer);
}

void convert(const TfcPrjPrms &rootObj, TFC_prjPrms& tdsObj) 
{
    tdsObj.xy0   = rootObj.getxy0();
    tdsObj.dxy   = rootObj.getxy0();
    tdsObj.dz    = rootObj.getxy0();
    tdsObj.tower = rootObj.getxy0();
    tdsObj.layer = rootObj.getxy0();
}

void convert(const TFC_prjList &tdsObj, TfcPrjList& rootObj) 
{
    rootObj.initialize(tdsObj.node.flnk, tdsObj.node.blnk);
}

void convert(const TfcPrjList &rootObj, TFC_prjList& tdsObj) 
{
    tdsObj.node.flnk = (TFC_prjNode*)rootObj.getFlnk();
    tdsObj.node.blnk = (TFC_prjNode*)rootObj.getBlnk();
}

void convert(const TFC_prjDir &tdsObj, TfcPrjDir& rootObj) 
{
    rootObj.initialize(tdsObj.idx, tdsObj.xCnt, tdsObj.yCnt);
}

void convert(const TfcPrjDir &rootObj, TFC_prjDir& tdsObj) 
{
    tdsObj.idx = rootObj.getIndex();
    tdsObj.xCnt = rootObj.getXCount();
    tdsObj.yCnt = rootObj.getYCount();
}

void convert( const TFC_prjs* tdsObj, TfcProjectionCol& rootObj) 
{
    int i;
    TfcPrjDir  dirRoot[16];
    TfcPrjList listTopRoot[32];
    TfcPrjList listBotRoot[32];

    for (int idx=0; idx<16; idx++)
    {
        convert(tdsObj->dir[idx], dirRoot[idx]);
        convert(tdsObj->top[0][idx], listTopRoot[idx]);
        convert(tdsObj->top[1][idx], listTopRoot[idx+16]);
        convert(tdsObj->bot[0][idx], listBotRoot[idx]);
        convert(tdsObj->bot[1][idx], listBotRoot[idx+16]);
    }

    TfcPrj *prjColRoot = 0;
    prjColRoot = new TfcPrj[tdsObj->curCnt];

    for (i=0; i<tdsObj->curCnt; i++) 
        convert(tdsObj->prjs[i], prjColRoot[i]);

    rootObj.initialize(tdsObj->maxCnt, tdsObj->curCnt, tdsObj->twrMsk, 
        dirRoot, listTopRoot, listBotRoot, prjColRoot);

    delete [] prjColRoot;
}

void convert ( const TfcProjectionCol& rootObj, TFC_prjs*& tdsObj ) 
{
    tdsObj->maxCnt = rootObj.getMaxCount();
    tdsObj->curCnt = rootObj.getCurrentCount();
    tdsObj->twrMsk = rootObj.getTwrMsk();

    const TfcPrjDir*  dirRoot     = rootObj.getProjectionDirCol();
    const TfcPrjList* topListRoot = rootObj.getTopList();
    const TfcPrjList* botListRoot = rootObj.getBotList();
    UInt_t i;
    for (i=0;i<16;i++)
    {
        convert(dirRoot[i], tdsObj->dir[i]);
        convert(topListRoot[i],    tdsObj->top[0][i]);
        convert(topListRoot[i+16], tdsObj->top[1][i]);
        convert(botListRoot[i],    tdsObj->bot[0][i]);
        convert(botListRoot[i+16], tdsObj->bot[1][i]);
    }
 
    const TfcPrj* projectionArrRoot = rootObj.getProjectionCol();
    // Check to be sure the number of elements in the array doesn't exceed 
    // 1000, that's the setting for this fixed array on the TDS
    unsigned int numProj = rootObj.getCurrentCount();
    if (rootObj.getCurrentCount() >= 1000) 
        numProj=1000;
    for (i=0; i<numProj; i++) 
       convert(projectionArrRoot[i], tdsObj->prjs[i]); 

}

void convert( const std::vector<OnboardFilterTds::track>& tdsObj, 
                TObjArray& rootObj) 
{
    std::vector<OnboardFilterTds::track>::const_iterator trackIt;
    for(trackIt = tdsObj.begin(); trackIt != tdsObj.end(); trackIt++) {
        ObfTrack *rootTrack = new ObfTrack();
        rootTrack->initialize(trackIt->numhits, trackIt->phi_rad,
                   trackIt->theta_rad, trackIt->length, trackIt->lowCoord, 
                   trackIt->highCoord, trackIt->exLowCoord, 
                   trackIt->exHighCoord);
        rootObj.Add(rootTrack);
    }
}

void convert ( const TObjArray* rootObj, 
                 std::vector<OnboardFilterTds::track>& tdsObj) 
{
    TIter rootIt(rootObj);
    ObfTrack* rootTrack;
    while((rootTrack = (ObfTrack*)rootIt.Next())) 
    {
        OnboardFilterTds::track curTrack;
        curTrack.numhits = rootTrack->getNumHits();
        curTrack.phi_rad = rootTrack->getPhi();
        curTrack.theta_rad = rootTrack->getTheta();
        curTrack.lowCoord = rootTrack->getLowCoord();
        curTrack.highCoord = rootTrack->getHighCoord();
        curTrack.exLowCoord = rootTrack->getExLowCoord();
        curTrack.exHighCoord = rootTrack->getExHighCoord();
        curTrack.length = rootTrack->getLength();
        tdsObj.push_back(curTrack);
    }
}


void convert( OnboardFilterTds::FilterStatus& tdsObj, 
                FilterStatus& rootObj) 
{
    Int_t ebfSize = 0;
    rootObj.init(tdsObj.get(),tdsObj.getStageEnergy(),tdsObj.getTcids(),
                 ebfSize);
    rootObj.initGem(tdsObj.getGemThrTkr(), tdsObj.getGemCalHiLo(),
            tdsObj.getGemCondsumCno(),
            tdsObj.getGemAcd_vetoes_XZ(),tdsObj.getGemAcd_vetoes_YZ(),
            tdsObj.getGemAcd_vetoes_XY(), tdsObj.getGemAcd_vetoes_RU(),
            tdsObj.getGemLivetime(), tdsObj.getGemTrgtime(), 
            tdsObj.getGemPpstime(), tdsObj.getGemDiscarded(),
            tdsObj.getGemPrescaled(), tdsObj.getGemSent());

    int xz, yz, xy;
    tdsObj.getAcdMap(xz, yz, xy);
    const Int_t *acdStatus = tdsObj.getAcdStatus();

    rootObj.initAcd(xz, yz, xy, tdsObj.getVetoWord(),acdStatus);

    rootObj.initTkr(tdsObj.getXcapture(), tdsObj.getYcapture(),
                    tdsObj.getXY00(), tdsObj.getXY11(), tdsObj.getXY22(),
                    tdsObj.getXY33(), tdsObj.getTmsk());

    rootObj.initLayers(tdsObj.getLayers());

    Int_t numLogsHit = 0; //  can't find this right now
    rootObj.initCal(numLogsHit, tdsObj.getLayerEnergy());

    rootObj.initSeparation(tdsObj.getSeparation());
    
    int xHits, yHits;
    double slopeXZ,slopeYZ,intXZ,intYZ;
    int nGrbHitsX,nGrbHitsY;
    double grbSlpX,grbSlpY;
    tdsObj.getBestTrack(xHits,yHits,slopeXZ,slopeYZ,intXZ,intYZ,nGrbHitsX,nGrbHitsY,grbSlpX,grbSlpY);
    rootObj.initBestTrack(xHits,yHits,slopeXZ,slopeYZ,intXZ,intYZ);

    TObjArray trackRootCol;
    if (tdsObj.tracksExist()) {
        const std::vector<OnboardFilterTds::track> trackTdsCol = tdsObj.getTracks();
        convert(trackTdsCol, trackRootCol);
    }
    rootObj.initTracks(trackRootCol);

    TfcProjectionCol projRootCol;
    const TFC_prjs* projTds = tdsObj.getProjections();
    convert(projTds, projRootCol);
    rootObj.initProjectionCol(projRootCol);
};

void convert( const FilterStatus& rootObj, 
              OnboardFilterTds::FilterStatus& tdsObj) 
{
    tdsObj.set(rootObj.getStatus());
    tdsObj.setStageEnergy(rootObj.getStageEnergy());
    tdsObj.setTcids(rootObj.getTcids());

    tdsObj.setGemThrTkr(rootObj.getGemThrTkr());
    tdsObj.setGemCalHiLo(rootObj.getGemCalHiLo());
    tdsObj.setGemCondsumCno(rootObj.getGemCondsumCno());
    tdsObj.setGemAcd_vetoes_XZ(rootObj.getGemAcd_vetoes_XZ());
    tdsObj.setGemAcd_vetoes_YZ(rootObj.getGemAcd_vetoes_YZ());
    tdsObj.setGemAcd_vetoes_XY(rootObj.getGemAcd_vetoes_XY());
    tdsObj.setGemAcd_vetoes_RU(rootObj.getGemAcd_vetoes_RU());
    tdsObj.setGemLivetime(rootObj.getGemLivetime());
    tdsObj.setGemTrgtime(rootObj.getGemTrgtime());
    tdsObj.setGemPpstime(rootObj.getGemPpstime());
    tdsObj.setGemDiscarded(rootObj.getGemDiscarded());
    tdsObj.setGemPrescaled(rootObj.getGemPrescaled());
    tdsObj.setGemSent(rootObj.getGemSent());
    tdsObj.setBestTrack(rootObj.getXhitsBestTrack(),rootObj.getYhitsBestTrack(),
        rootObj.getSlopeXzBestTrack(), rootObj.getSlopeYzBestTrack(),rootObj.getIntXzBestTrack(),
        rootObj.getIntYzBestTrack(),0,0,0.,0.);

    Int_t xz, yz, xy;
    rootObj.getAcdMap(xz,yz,xy);
    tdsObj.setAcdMap(xz,yz,xy);
      
    Int_t iAcd;
    const Int_t* statusRoot = rootObj.getAcdStatus();
    for (iAcd=0; iAcd < 16; iAcd++) 
    {
        tdsObj.setAcdStatus(iAcd, statusRoot[iAcd]);
    }

    tdsObj.setVetoWord(rootObj.getVetoWord());
    
    tdsObj.setLayers(rootObj.getLayerCol());

    tdsObj.setSeparation(rootObj.getSeparation());

    //Ask Richard to provide suitable set method using floats rather than int
    //tdsObj->setLayerEnergy(rootObj.getLayerEnergy());

    tdsObj.setCapture(rootObj.getXcapture(), rootObj.getYcapture());

    tdsObj.setXY(rootObj.getXY00(), rootObj.getXY11(), rootObj.getXY22(),
                 rootObj.getXY33());

    tdsObj.setTmsk(rootObj.getTmsk());
};

/// Fill persitent verison of ObfFilterStatus 
void convert( const OnboardFilterTds::ObfFilterStatus& tdsObj, ObfFilterStatus& rootObj)
{
    // Initialize
    rootObj.initialize();

    // Now step through and initialize the status objects one by one
    const OnboardFilterTds::IObfStatus* tdsStatus = 0;

    // We do this one by one explicitly for now. Start with the results of the gamma filter
    ObfGammaStatus* gammaStatus = new ObfGammaStatus();
    if ((tdsStatus = tdsObj.getFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter)))
    {
        convert(*tdsStatus, *gammaStatus);
    }
    rootObj.addFilterStatus(ObfFilterStatus::GammaFilter, gammaStatus);

    // HFC Filter is next
    ObfHipStatus* hfcStatus = new ObfHipStatus();
    if ((tdsStatus = tdsObj.getFilterStatus(OnboardFilterTds::ObfFilterStatus::HIPFilter)))
    {
        convert(*tdsStatus, *hfcStatus);
    }
    rootObj.addFilterStatus(ObfFilterStatus::HipFilter, hfcStatus);

    // MIP Filter is next
    ObfMipStatus* mipStatus = new ObfMipStatus();
    if ((tdsStatus = tdsObj.getFilterStatus(OnboardFilterTds::ObfFilterStatus::MIPFilter)))
    {
        convert(*tdsStatus, *mipStatus);
    }
    rootObj.addFilterStatus(ObfFilterStatus::MipFilter, mipStatus);

    // DFC Filter is next
    ObfDgnStatus* dfcStatus = new ObfDgnStatus();
    if ((tdsStatus = tdsObj.getFilterStatus(OnboardFilterTds::ObfFilterStatus::DGNFilter)))
    {
        convert(*tdsStatus, *dfcStatus);
    }
    rootObj.addFilterStatus(ObfFilterStatus::DgnFilter, dfcStatus);

    return;
}

/// Fill transient verison of ObfFilterStatus 
void convert( const ObfFilterStatus& rootObj, OnboardFilterTds::ObfFilterStatus& tdsObj) 
{
    // Now step through and initialize the status objects one by one
    const IObfStatus* rootStatus = 0;

    if (rootStatus = rootObj.getFilterStatus(ObfFilterStatus::GammaFilter))
    {
        OnboardFilterTds::ObfGammaStatus* tdsStatus = new OnboardFilterTds::ObfGammaStatus(0,0,0,0);
        convert(*rootStatus, *tdsStatus);
        tdsObj.addFilterStatus(OnboardFilterTds::ObfFilterStatus::GammaFilter, tdsStatus);
    }

    if (rootStatus = rootObj.getFilterStatus(ObfFilterStatus::HipFilter))
    {
        OnboardFilterTds::ObfHipStatus* tdsStatus = new OnboardFilterTds::ObfHipStatus(0,0,0);
        convert(*rootStatus, *tdsStatus);
        tdsObj.addFilterStatus(OnboardFilterTds::ObfFilterStatus::HIPFilter, tdsStatus);
    }

    if (rootStatus = rootObj.getFilterStatus(ObfFilterStatus::MipFilter))
    {
        OnboardFilterTds::ObfMipStatus* tdsStatus = new OnboardFilterTds::ObfMipStatus(0,0,0);
        convert(*rootStatus, *tdsStatus);
        tdsObj.addFilterStatus(OnboardFilterTds::ObfFilterStatus::MIPFilter, tdsStatus);
    }

    if (rootStatus = rootObj.getFilterStatus(ObfFilterStatus::DgnFilter))
    {
        OnboardFilterTds::ObfDgnStatus* tdsStatus = new OnboardFilterTds::ObfDgnStatus(0,0,0);
        convert(*rootStatus, *tdsStatus);
        tdsObj.addFilterStatus(OnboardFilterTds::ObfFilterStatus::DGNFilter, tdsStatus);
    }

    return;
}

/// Fill persistent verison of ObfGammaStatus 
void convert( const OnboardFilterTds::IObfStatus& tdsObj, ObfGammaStatus& rootObj )
{
    const OnboardFilterTds::ObfGammaStatus& gamObj = dynamic_cast<const OnboardFilterTds::ObfGammaStatus&>(tdsObj);

    rootObj = ObfGammaStatus(gamObj.getFilterId(), gamObj.getStatusWord(), gamObj.getFiltersb(), gamObj.getEnergy());

    return;
}

/// Fill transient verison of ObfGammaStatus 
void convert( const IObfStatus& rootObj, OnboardFilterTds::ObfGammaStatus& tdsObj) 
{
    const ObfGammaStatus& gamObj = dynamic_cast<const ObfGammaStatus&>(rootObj);

    tdsObj = OnboardFilterTds::ObfGammaStatus(gamObj.getFilterId(), gamObj.getStatusWord(), gamObj.getFiltersb(), gamObj.getEnergy());

    return;
}

/// Fill persistent verison of ObfHipStatus 
void convert( const OnboardFilterTds::IObfStatus& tdsObj, ObfHipStatus& rootObj)
{
    rootObj = ObfHipStatus(tdsObj.getFilterId(), tdsObj.getStatusWord(), tdsObj.getFiltersb());

    return;
}

/// Fill transient verison of ObfHipStatus 
void convert( const IObfStatus& rootObj, OnboardFilterTds::ObfHipStatus& tdsObj) 
{
    tdsObj = OnboardFilterTds::ObfHipStatus(rootObj.getFilterId(), rootObj.getStatusWord(), rootObj.getFiltersb());

    return;
}

/// Fill persistent verison of ObfMipStatus 
void convert( const OnboardFilterTds::IObfStatus& tdsObj, ObfMipStatus& rootObj)
{
    rootObj = ObfMipStatus(tdsObj.getFilterId(), tdsObj.getStatusWord(), tdsObj.getFiltersb());

    return;
}

/// Fill transient verison of ObfMipStatus 
void convert( const IObfStatus& rootObj, OnboardFilterTds::ObfMipStatus& tdsObj) 
{
    tdsObj = OnboardFilterTds::ObfMipStatus(rootObj.getFilterId(), rootObj.getStatusWord(), rootObj.getFiltersb());

    return;
}

/// Fill persistent verison of ObfDFCStatus 
void convert( const OnboardFilterTds::IObfStatus& tdsObj, ObfDgnStatus& rootObj)
{
    rootObj = ObfDgnStatus(tdsObj.getFilterId(), tdsObj.getStatusWord(), tdsObj.getFiltersb());

    return;
}

/// Fill transient verison of ObfDFCStatus 
void convert( const IObfStatus& rootObj, OnboardFilterTds::ObfDgnStatus& tdsObj) 
{
    tdsObj = OnboardFilterTds::ObfDgnStatus(rootObj.getFilterId(), rootObj.getStatusWord(), rootObj.getFiltersb());

    return;
}

  /// Fill persistent verison of ObfFilterTrack 
void convert( const OnboardFilterTds::ObfFilterTrack& tdsObj, ObfFilterTrack& rootObj)
{
    rootObj.initialize(tdsObj.get_nXhits(), 
                       tdsObj.get_nYhits(), 
                       tdsObj.get_xInt(), 
                       tdsObj.get_yInt(), 
                       tdsObj.get_z(), 
                       tdsObj.get_slpXZ(), 
                       tdsObj.get_slpYZ() );
}
  
/// Fill transient verison of ObfFilterStatus 
void convert( const ObfFilterTrack& rootObj, OnboardFilterTds::ObfFilterTrack& tdsObj)
{
    tdsObj.initialize(rootObj.get_nXhits(), 
                      rootObj.get_nYhits(), 
                      rootObj.get_xInt(), 
                      rootObj.get_yInt(), 
                      rootObj.get_z(), 
                      rootObj.get_slpXZ(), 
                      rootObj.get_slpYZ() );
}


}
