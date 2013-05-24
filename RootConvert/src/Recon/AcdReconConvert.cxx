
#include <RootConvert/Recon/AcdReconConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include <CLHEP/Matrix/Matrix.h>
#include <TMatrixD.h>
#include <TMatrixF.h>

namespace RootPersistence {

using namespace CLHEP;

void convert( const Event::AcdRecon & tdsAcdRec, AcdRecon & rootAcdRec ) {
    
  std::vector<AcdId> rootIdCol ;
  std::vector<idents::AcdId>::const_iterator tdsIdIt ;
  for ( tdsIdIt = tdsAcdRec.getIdCol().begin() ;
        tdsIdIt != tdsAcdRec.getIdCol().end() ; tdsIdIt++ ) {    
    AcdId acdId = convert( *tdsIdIt );
    rootIdCol.push_back( acdId );
  }

  idents::AcdId tdsDocaId = tdsAcdRec.getMinDocaId();
  AcdId rootDocaId(tdsDocaId.layer(), tdsDocaId.face(),
    tdsDocaId.row(), tdsDocaId.column());
  idents::AcdId tdsActDistId = tdsAcdRec.getMaxActDist3DId();
  AcdId rootActDistId(tdsActDistId.layer(),tdsActDistId.face(),
    tdsActDistId.row(), tdsActDistId.column());
  idents::AcdId tdsRibActDistId = tdsAcdRec.getRibbonActiveDistId();
  AcdId rootRibActDistId(tdsRibActDistId.layer(), tdsRibActDistId.face(),
    tdsRibActDistId.row(), tdsRibActDistId.column());

  rootAcdRec.initialize(
    tdsAcdRec.getEnergy(), tdsAcdRec.getRibbonEnergy(),
    tdsAcdRec.getTileCount(), tdsAcdRec.getRibbonCount(),
        tdsAcdRec.getGammaDoca(), tdsAcdRec.getDoca(), rootDocaId, 
    tdsAcdRec.getActiveDist3D(), rootActDistId, 
    tdsAcdRec.getRibbonActiveDist(), rootRibActDistId,
    tdsAcdRec.getRowDocaCol(), tdsAcdRec.getRowActDist3DCol(),
        rootIdCol, tdsAcdRec.getEnergyCol(), tdsAcdRec.getCornerDoca());

  idents::AcdId tdsActDistId_Down = tdsAcdRec.getMaxActDist3DId_Down();
  AcdId rootActDistId_Down = convert(tdsActDistId_Down);
  rootAcdRec.initActDist_Down(tdsAcdRec.getActiveDist3D_Down(), 
                              rootActDistId_Down,
                              tdsAcdRec.getRowActDist3DCol_Down());
    
  const Event::AcdTkrIntersectionCol& tdsAcdTkrIntersects = tdsAcdRec.getAcdTkrIntersectionCol() ;
  AcdTkrIntersection rootAcdInter;
  for ( Event::AcdTkrIntersectionCol::const_iterator itrInter = tdsAcdTkrIntersects.begin() ;
        itrInter != tdsAcdTkrIntersects.end(); itrInter++ ) {
    const Event::AcdTkrIntersection* tdsAcdInter = *itrInter;
    convert(*tdsAcdInter,rootAcdInter);
    rootAcdRec.addAcdTkrIntersection(rootAcdInter);
  }

  const Event::AcdTkrPocaCol& tdsAcdTkrPocas = tdsAcdRec.getAcdTkrPocaCol() ;
  AcdTkrPoca rootAcdPoca;
  for ( Event::AcdTkrPocaCol::const_iterator itrPoca = tdsAcdTkrPocas.begin();
        itrPoca != tdsAcdTkrPocas.end(); itrPoca++ ) {
    const Event::AcdTkrPoca* acdTkrPocaTDS = *itrPoca;
    convert(*acdTkrPocaTDS,rootAcdPoca);
    rootAcdRec.addAcdTkrPoca(rootAcdPoca);
  }  
  
  const Event::AcdHitCol& tdsAcdHits = tdsAcdRec.getAcdHitCol() ;    
  AcdHit rootAcdHit;
  for ( Event::AcdHitCol::const_iterator itrHit = tdsAcdHits.begin();
        itrHit != tdsAcdHits.end(); itrHit++ ) {
    const Event::AcdHit* acdHitTDS = *itrHit;
    convert(*acdHitTDS,rootAcdHit);
    rootAcdRec.addAcdHit(rootAcdHit);
  } 

  const Event::AcdTkrHitPocaCol& tdsAcdTkrHitPocas = tdsAcdRec.getAcdTkrHitPocaCol() ;
  AcdTkrHitPoca rootAcdHitPoca;
  for ( Event::AcdTkrHitPocaCol::const_iterator itrHitPoca = tdsAcdTkrHitPocas.begin();
        itrHitPoca != tdsAcdTkrHitPocas.end(); itrHitPoca++ ) {
    const Event::AcdTkrHitPoca* acdTkrHitPocaTDS = *itrHitPoca;
    convert(*acdTkrHitPocaTDS,rootAcdHitPoca);
    rootAcdRec.addAcdTkrHitPoca(rootAcdHitPoca);
  }  
  
  const Event::AcdTkrGapPocaCol& tdsAcdTkrGapPocas = tdsAcdRec.getAcdTkrGapPocaCol() ;
  AcdTkrGapPoca rootAcdGapPoca;
  for ( Event::AcdTkrGapPocaCol::const_iterator itrGapPoca = tdsAcdTkrGapPocas.begin();
        itrGapPoca != tdsAcdTkrGapPocas.end(); itrGapPoca++ ) {
    const Event::AcdTkrGapPoca* acdTkrGapPocaTDS = *itrGapPoca;
    convert(*acdTkrGapPocaTDS,rootAcdGapPoca);
    rootAcdRec.addAcdTkrGapPoca(rootAcdGapPoca);
  }  

  const Event::AcdTkrPointCol& tdsAcdTkrPoints = tdsAcdRec.getAcdTkrPointCol() ;
  AcdTkrPoint rootAcdPoint;
  for ( Event::AcdTkrPointCol::const_iterator itrPoint = tdsAcdTkrPoints.begin();
        itrPoint != tdsAcdTkrPoints.end(); itrPoint++ ) {
    const Event::AcdTkrPoint* acdTkrPointTDS = *itrPoint;
    convert(*acdTkrPointTDS,rootAcdPoint);
    rootAcdRec.addAcdTkrPoint(rootAcdPoint);
  }  

  const Event::AcdSplashVarsCol& tdsAcdSplashVars = tdsAcdRec.getAcdSplashVarsCol() ;
  AcdSplashVars rootAcdSplash;
  for ( Event::AcdSplashVarsCol::const_iterator itrSplash = tdsAcdSplashVars.begin();
        itrSplash != tdsAcdSplashVars.end(); itrSplash++ ) {
    const Event::AcdSplashVars* acdSplashVarsTDS = *itrSplash;
    convert(*acdSplashVarsTDS,rootAcdSplash);
    rootAcdRec.addAcdSplashVar(rootAcdSplash);
  }  
}
 
void convert( const AcdRecon & rootAcdRec, Event::AcdRecon & tdsAcdRec )
{
    std::vector<Event::AcdTkrIntersection*> acdTkrIntersections;
    int nInter = rootAcdRec.nAcdIntersections() ;    
    for ( int iInter(0); iInter < nInter; iInter++ ) {
      const AcdTkrIntersection * rootAcdInter = rootAcdRec.getAcdTkrIntersection(iInter);      
      Event::AcdTkrIntersection* tdsAcdInter = new Event::AcdTkrIntersection();
      convert(*rootAcdInter,*tdsAcdInter);
      acdTkrIntersections.push_back(tdsAcdInter) ;
    }

    std::vector<Event::AcdTkrPoca*> acdTkrPocas;
    int nPoca = rootAcdRec.nAcdTkrPoca();
    for ( int iPoca(0); iPoca < nPoca; iPoca++ ) {
      const AcdTkrPoca* acdPocaRoot = rootAcdRec.getAcdTkrPoca(iPoca);
      Event::AcdTkrPoca* acdPocaTds = new Event::AcdTkrPoca();
      convert(*acdPocaRoot,*acdPocaTds);
      acdTkrPocas.push_back(acdPocaTds);
    }
    
    // do the same w/ the hits   
    std::vector<Event::AcdHit*> acdHits;
    int nHit = rootAcdRec.nAcdHit();
    for ( int iHit(0); iHit < nHit; iHit++ ) {
      const AcdHit* acdHitRoot = rootAcdRec.getAcdHit(iHit);    
      Event::AcdHit* acdHitTds = new Event::AcdHit();
      convert(*acdHitRoot,*acdHitTds);
      acdHits.push_back(acdHitTds);
    }

    std::vector<Event::AcdTkrHitPoca*> acdTkrHitPocas;
    int nHitPoca = rootAcdRec.nAcdTkrHitPoca();
    for ( int iHitPoca(0); iHitPoca < nHitPoca; iHitPoca++ ) {
      const AcdTkrHitPoca* acdPocaRoot = rootAcdRec.getAcdTkrHitPoca(iHitPoca);
      Event::AcdTkrHitPoca* acdPocaTds = new Event::AcdTkrHitPoca();
      convert(*acdPocaRoot,*acdPocaTds);
      acdTkrHitPocas.push_back(acdPocaTds);
    }

    std::vector<Event::AcdTkrGapPoca*> acdTkrGapPocas;
    int nGapPoca = rootAcdRec.nAcdTkrGapPoca();
    for ( int iGapPoca(0); iGapPoca < nGapPoca; iGapPoca++ ) {
      const AcdTkrGapPoca* acdPocaRoot = rootAcdRec.getAcdTkrGapPoca(iGapPoca);
      Event::AcdTkrGapPoca* acdPocaTds = new Event::AcdTkrGapPoca();
      convert(*acdPocaRoot,*acdPocaTds);
      acdTkrGapPocas.push_back(acdPocaTds);
    }

    std::vector<Event::AcdTkrPoint*> acdTkrPoints;
    int nPoint = rootAcdRec.nAcdTkrPoint();
    for ( int iPoint(0); iPoint < nPoint; iPoint++ ) {
      const AcdTkrPoint* acdPocaRoot = rootAcdRec.getAcdTkrPoint(iPoint);
      Event::AcdTkrPoint* acdPocaTds = new Event::AcdTkrPoint();
      convert(*acdPocaRoot,*acdPocaTds);
      acdTkrPoints.push_back(acdPocaTds);
    }

    std::vector<Event::AcdSplashVars*> acdSplashVars;
    int nSplash = rootAcdRec.nAcdSplash();
    for ( int iSplash(0); iSplash < nSplash; iSplash++ ) {
      const AcdSplashVars* acdSplashRoot = rootAcdRec.getAcdSplashVars(iSplash);
      Event::AcdSplashVars* acdSplashTds = new Event::AcdSplashVars();
      convert(*acdSplashRoot,*acdSplashTds);
      acdSplashVars.push_back(acdSplashTds);
    }
  
    const AcdId rootDocaId = rootAcdRec.getMinDocaId();
    const idents::AcdId tdsDocaId(rootDocaId.getLayer(), rootDocaId.getFace(), 
        rootDocaId.getRow(), rootDocaId.getColumn());

    const AcdId rootActDistId = rootAcdRec.getMaxActDistId();
    const idents::AcdId tdsActDistId(rootActDistId.getLayer(), 
        rootActDistId.getFace(), rootActDistId.getRow(), 
        rootActDistId.getColumn());

    const AcdId rootRibActDistId = rootAcdRec.getRibbonActDistId();
    const idents::AcdId tdsRibActDistId(rootRibActDistId.getLayer(),
        rootRibActDistId.getFace(), rootRibActDistId.getRow(),
        rootRibActDistId.getColumn());

    std::vector<idents::AcdId> tdsIdCol ;
    std::vector<AcdId>::const_iterator rootIdIt ;
    for ( rootIdIt = rootAcdRec.getIdCol().begin() ;
          rootIdIt != rootAcdRec.getIdCol().end() ;
          rootIdIt++) {
        tdsIdCol.push_back(
            idents::AcdId(
                rootIdIt->getLayer(), 
                rootIdIt->getFace(), rootIdIt->getRow(), 
                rootIdIt->getColumn()
                )
            ) ;
    }
    
    std::vector<double> tdsEnergyCol = rootAcdRec.getEnergyCol();
    std::vector<double> tdsRowActDist3D = rootAcdRec.getRowActDistCol();

    // Note that ActDist stored in ROOT is only the new 3D Active Dist.
    tdsAcdRec.init(
        rootAcdRec.getEnergy(), 
        rootAcdRec.getRibbonEnergy(),
        rootAcdRec.getTileCount(),
        rootAcdRec.getRibbonCount(),
        rootAcdRec.getGammaDoca(),
        rootAcdRec.getDoca(),
        tdsDocaId,
        rootAcdRec.getActiveDist(),
        tdsActDistId, 
        rootAcdRec.getRowDocaCol(),
        rootAcdRec.getRowActDistCol(),
        tdsIdCol, 
        tdsEnergyCol,
        rootAcdRec.getActiveDist(), 
        tdsActDistId,
        rootAcdRec.getRowActDistCol(),
        acdTkrIntersections,
        acdTkrPocas,
        acdHits,
        acdTkrHitPocas,
        acdTkrGapPocas,
        acdTkrPoints,
        acdSplashVars,
        rootAcdRec.getRibbonActiveDist(),
        tdsRibActDistId, 
        rootAcdRec.getCornerDoca()
        ) ;

     AcdId rootActDistId_down = rootAcdRec.getMaxActDistId_Down();
     idents::AcdId tdsActDistId_down = convert(rootActDistId_down);
     tdsAcdRec.initActDist3D_down(rootAcdRec.getActiveDist_Down(),
                                  tdsActDistId_down,
                                  rootAcdRec.getRowActDistCol_Down());
}


  void convert( const Event::AcdReconV2& tds, AcdReconV2& root) {
    convert( tds.getEventTopology(), root.getEventTopology() );

    const Event::AcdHitCol& tdsAcdHits = tds.getHitCol() ;    
    AcdHit rootAcdHit;
    for ( Event::AcdHitCol::const_iterator itrHit = tdsAcdHits.begin();
          itrHit != tdsAcdHits.end(); itrHit++ ) {
      const Event::AcdHit* acdHitTDS = *itrHit;
      AcdHit* addedHit = root.addAcdHit(rootAcdHit);
      convert(*acdHitTDS,*addedHit);      
    } 

    const Event::AcdTkrAssocCol& tdsAcdTkrAssocs = tds.getTkrAssocCol() ;    
    AcdAssoc rootAcdTkrAssoc;
    for ( Event::AcdTkrAssocCol::const_iterator itrAssoc = tdsAcdTkrAssocs.begin();
          itrAssoc != tdsAcdTkrAssocs.end(); itrAssoc++ ) {
      const Event::AcdAssoc* acdTkrAssocTDS = *itrAssoc;
      AcdAssoc* addedAssoc = root.addAcdTkrAssoc(rootAcdTkrAssoc);
      convert(*acdTkrAssocTDS,*addedAssoc);      
    } 

    const Event::AcdCalAssocCol& tdsAcdCalAssocs = tds.getCalAssocCol() ;    
    AcdAssoc rootAcdCalAssoc;
    for ( Event::AcdCalAssocCol::const_iterator itrAssoc = tdsAcdCalAssocs.begin();
          itrAssoc != tdsAcdCalAssocs.end(); itrAssoc++ ) {
      const Event::AcdAssoc* acdCalAssocTDS = *itrAssoc;
      AcdAssoc* addedAssoc = root.addAcdCalAssoc(rootAcdCalAssoc);
      convert(*acdCalAssocTDS,*addedAssoc);      
    } 

  }
  
  void convert( const AcdReconV2 & root, Event::AcdReconV2& tds) {
    convert( root.getEventTopology(), tds.getEventTopology() );

    // do the same w/ the hits   
    const TClonesArray& rootHits = root.getHitCol();
    for ( Int_t iHit(0); iHit < rootHits.GetEntriesFast(); iHit++ ) {
      const AcdHit* acdHitRoot = (const AcdHit*)(rootHits.UncheckedAt(iHit));
      Event::AcdHit* acdHitTds = new Event::AcdHit();
      convert(*acdHitRoot,*acdHitTds);
      tds.getHitCol().push_back(acdHitTds);
    }

    const TClonesArray& rootTkrAssocs = root.getTkrAssocCol();
    for ( Int_t iAssoc(0); iAssoc < rootTkrAssocs.GetEntriesFast(); iAssoc++ ) {
      const AcdAssoc* acdTkrAssocRoot = (const AcdAssoc*)(rootTkrAssocs.UncheckedAt(iAssoc));
      Event::AcdAssoc* acdTkrAssocTds = new Event::AcdAssoc();
      convert(*acdTkrAssocRoot,*acdTkrAssocTds);
      tds.getTkrAssocCol().push_back(acdTkrAssocTds);
    }

    const TClonesArray& rootCalAssocs = root.getCalAssocCol();
    for ( Int_t iAssoc(0); iAssoc < rootCalAssocs.GetEntriesFast(); iAssoc++ ) {
      const AcdAssoc* acdCalAssocRoot = (const AcdAssoc*)(rootCalAssocs.UncheckedAt(iAssoc));
      Event::AcdAssoc* acdCalAssocTds = new Event::AcdAssoc();
      convert(*acdCalAssocRoot,*acdCalAssocTds);
      tds.getCalAssocCol().push_back(acdCalAssocTds);
    }

  }   

  void convert( const Event::AcdHit& tds, AcdHit& root) {
    AcdId theId = convert( tds.getAcdId());
    root.set( theId,
              tds.getFlags(Event::AcdHit::A), tds.getFlags(Event::AcdHit::B),
              tds.getPha(Event::AcdHit::A), tds.getPha(Event::AcdHit::B),
              tds.mips(Event::AcdHit::A), tds.mips(Event::AcdHit::B) ); 
  } 
  
  void convert( const AcdHit& root, Event::AcdHit& tds) {
    idents::AcdId theId = convert( root.getId());
    tds.set( theId, 
             root.getFlags(AcdHit::A), root.getFlags(AcdHit::B),
             root.getPha(AcdHit::A), root.getPha(AcdHit::B),
             root.getMips(AcdHit::A), root.getMips(AcdHit::B) ); 
  }

  void convert( const Event::AcdEventTopology& tds, AcdEventTopology& root) {
    UInt_t tileCountSideRow[4];
    UInt_t tileCountSideFace[4];
    UInt_t vetoCountSideRow[4];
    UInt_t vetoCountSideFace[4];
    Float_t tileEnergySideRow[4];
    Float_t tileEnergySideFace[4];
    for ( UInt_t i(0); i < 4; i++ ) {
      tileCountSideRow[i] = tds.getTileCountSideRow(i);
      tileCountSideFace[i] = tds.getTileCountSideFace(i+1);
      vetoCountSideRow[i] = tds.getVetoCountSideRow(i);
      vetoCountSideFace[i] = tds.getVetoCountSideFace(i+1);
      tileEnergySideRow[i] = tds.getTileEnergySideRow(i);
      tileEnergySideFace[i] = tds.getTileEnergySideFace(i+1);
    }

    root.set( tds.getTileCount(), tds.getRibbonCount(), tds.getVetoCount(), tds.getTileVeto(),
              tds.getTotalTileEnergy(), tds.getTotalRibbonEnergy(), tds.getTileEnergy(), tds.getRibbonEnergy(), 
              tds.getGhostTileEnergy(), tds.getGhostRibbonEnergy(), tds.getTriggerTileEnergy(), tds.getTriggerRibbonEnergy(), 
              tds.getTileCountTop(), tileCountSideRow, tileCountSideFace,
              tds.getVetoCountTop(), vetoCountSideRow, vetoCountSideFace,
              tds.getTileEnergyTop(), tileEnergySideRow, tileEnergySideFace,
              tds.getNSidesHit(), tds.getNSidesVeto());
  }

  void convert( const AcdEventTopology& root, Event::AcdEventTopology& tds){
    UInt_t tileCountSideRow[4];
    UInt_t tileCountSideFace[4];
    UInt_t vetoCountSideRow[4];
    UInt_t vetoCountSideFace[4];
    Float_t tileEnergySideRow[4];
    Float_t tileEnergySideFace[4];
    for ( UInt_t i(0); i < 4; i++ ) {
      tileCountSideRow[i] = root.getTileCountSideRow(i);
      tileCountSideFace[i] = root.getTileCountSideFace(i+1);
      vetoCountSideRow[i] = root.getVetoCountSideRow(i);
      vetoCountSideFace[i] = root.getVetoCountSideFace(i+1);
      tileEnergySideRow[i] = root.getTileEnergySideRow(i);
      tileEnergySideFace[i] = root.getTileEnergySideFace(i+1);
    }
    tds.set( root.getTileCount(), root.getRibbonCount(), root.getVetoCount(), root.getTileVeto(),
             root.getTotalTileEnergy(), root.getTotalRibbonEnergy(), root.getTileEnergy(), root.getRibbonEnergy(), 
             root.getGhostTileEnergy(), root.getGhostRibbonEnergy(), root.getTriggerTileEnergy(), root.getTriggerRibbonEnergy(), 
             root.getTileCountTop(), tileCountSideRow, tileCountSideFace,
             root.getVetoCountTop(), vetoCountSideRow, vetoCountSideFace,
             root.getTileEnergyTop(), tileEnergySideRow, tileEnergySideFace,
             root.getNSidesHit(), root.getNSidesVeto());
  }

  void convert( const Event::AcdAssoc& tds, AcdAssoc& root) {
    static TMatrixFSym covStart(5);
    static TMatrixFSym covEnd(5);

    TVector3 start = convert( tds.getStart() );
    TVector3 dir = convert3vector( tds.getDir() );
    convert( tds.getCovStart(), covStart );
    convert( tds.getCovEnd(), covEnd );
    root.set( tds.getTrackIndex(), tds.getUpward(), tds.getEnergy(),
              start, dir, 
              tds.getArcLength(), 
              covStart, covEnd,
              tds.getTkrSSDVeto(), tds.getCornerDoca(),
	      tds.GetEnergy15(),tds.GetEnergy30(),tds.GetEnergy45(),
	      tds.GetTriggerEnergy15(),tds.GetTriggerEnergy30(),tds.GetTriggerEnergy45());
	      

    AcdTkrHitPocaV2 rootAcdHitPoca;
    int nHitPoca = tds.nHitPoca();
    for ( int iHitPoca(0); iHitPoca < nHitPoca; iHitPoca++ ) {
      const Event::AcdTkrHitPoca* acdTkrHitPocaTDS = tds.getHitPoca(iHitPoca);
      AcdTkrHitPocaV2* addedHitPoca = root.addHitPoca(rootAcdHitPoca);;
      convert(*acdTkrHitPocaTDS,*addedHitPoca);
      
    }  
    
    AcdTkrGapPocaV2 rootAcdGapPoca;
    int nGapPoca = tds.nGapPoca();
    for ( int iGapPoca(0); iGapPoca < nGapPoca; iGapPoca++ ) {
      const Event::AcdTkrGapPoca* acdTkrGapPocaTDS = tds.getGapPoca(iGapPoca);
      AcdTkrGapPocaV2* addedGapPoca = root.addGapPoca(rootAcdGapPoca);
      convert(*acdTkrGapPocaTDS,*addedGapPoca);      
    }  
    
    AcdTkrPointV2 rootAcdPoint;
    const Event::AcdTkrPoint* acdTkrPointTDS = tds.getPoint();
    convert(*acdTkrPointTDS,rootAcdPoint);
    root.setPoint(rootAcdPoint);
  }
  
  void convert( const AcdAssoc& root, Event::AcdAssoc& tds) { 
    static HepSymMatrix covStart(5);
    static HepSymMatrix covEnd(5);
    Point start = convertPoint( root.getStart() );
    Vector dir = convertVector( root.getDir() );
    convert( root.getCovStart(), covStart );
    convert( root.getCovEnd(), covEnd );
    tds.set( root.getTrackIndex(), root.getUpward(), root.getEnergy(),
             start, dir, 
             root.getArcLength(), 
             covStart, covEnd,
             root.getTkrSSDVeto(), root.getCornerDoca(),
	     root.GetEnergy15(),root.GetEnergy30(),root.GetEnergy45(),
	     root.GetTriggerEnergy15(),root.GetTriggerEnergy30(),root.GetTriggerEnergy45());
	     

    int nHitPoca = root.nAcdHitPoca();
    for ( int iHitPoca(0); iHitPoca < nHitPoca; iHitPoca++ ) {
      const AcdTkrHitPocaV2* acdPocaRoot = root.getHitPoca(iHitPoca);
      Event::AcdTkrHitPoca* acdPocaTds = new Event::AcdTkrHitPoca();
      convert(*acdPocaRoot,*acdPocaTds);
      tds.addHitPoca(*acdPocaTds);
    }

    int nGapPoca = root.nAcdGapPoca();
    for ( int iGapPoca(0); iGapPoca < nGapPoca; iGapPoca++ ) {
      const AcdTkrGapPocaV2* acdPocaRoot = root.getGapPoca(iGapPoca);
      Event::AcdTkrGapPoca* acdPocaTds = new Event::AcdTkrGapPoca();
      convert(*acdPocaRoot,*acdPocaTds);
      tds.addGapPoca(*acdPocaTds);
    }

    AcdTkrPointV2& acdPocaRoot = const_cast<AcdTkrPointV2&>(root.getPoint());
    Event::AcdTkrPoint* acdPocaTds = new Event::AcdTkrPoint();
    convert(acdPocaRoot,*acdPocaTds);
    tds.setPoint(*acdPocaTds);
  }

  void convert( const Event::AcdTkrIntersection& tds, AcdTkrIntersection& root) {
    static TVector3 globalPosition;
    static Double_t localPosition[2];
    static TMatrixD covMatrix(2,2);    
    const Point& glbPos = tds.getGlobalPosition();
    globalPosition.SetXYZ(glbPos.x(),glbPos.y(),glbPos.z());
    localPosition[0] = tds.getLocalX();
    localPosition[1] = tds.getLocalY();
    covMatrix[0][0] = tds.getLocalXXCov();
    covMatrix[1][1] = tds.getLocalYYCov();
    covMatrix[0][1] = tds.getLocalXYCov();
    covMatrix[1][0] = tds.getLocalXYCov();
    AcdId acdId = convert( tds.getTileId() );
    root.initialize( acdId, tds.getTrackIndex(),
                     globalPosition, 
                     localPosition, covMatrix, 
                     tds.getArcLengthToIntersection(),
                     tds.getPathLengthInTile(),
                     tds.tileHit(),tds.getCosTheta());
  }

  
  void convert( const AcdTkrIntersection& root, Event::AcdTkrIntersection& tds) { 
    static Point globalPosition;
    static HepMatrix covMatrix(2,2);
    static double localPosition[2];
    const TVector3& glbPos = root.getGlobalPosition();
    globalPosition.set(glbPos.X(),glbPos.Y(),glbPos.Z());
    localPosition[0] = root.getLocalX();
    localPosition[1] = root.getLocalY();
    covMatrix[0][0] = root.getLocalXXCov();
    covMatrix[1][1] = root.getLocalYYCov();
    covMatrix[0][1] = covMatrix[1][0] = root.getLocalXYCov();
    idents::AcdId acdId = convert(root.getTileId());
    tds.set(acdId, root.getTrackIndex(),
            globalPosition,localPosition, covMatrix,
            root.getArcLengthToIntersection(),
            root.getPathLengthInTile(),
            root.tileHit(),root.getCosTheta()) ;
  }

  void convert( const Event::AcdTkrPoca& tds, AcdTkrPoca& root) {
    static TVector3 poca;
    const Point& pocaTDS = tds.poca();
    poca.SetXYZ(pocaTDS.x(),pocaTDS.y(),pocaTDS.z());
    const Event::TkrTrackParams& tkrParTds = tds.paramsAtPoca();
    TkrTrackParams paramsAtPoca(tkrParTds.getxPosition(), tkrParTds.getxSlope(), 
                                tkrParTds.getyPosition(), tkrParTds.getySlope(),
                                tkrParTds.getxPosxPos(),  tkrParTds.getxPosxSlp(), 
                                tkrParTds.getxPosyPos(),  tkrParTds.getxPosySlp(),
                                tkrParTds.getxSlpxSlp(),  tkrParTds.getxSlpyPos(), 
                                tkrParTds.getxSlpySlp(),  tkrParTds.getyPosyPos(), 
                                tkrParTds.getyPosySlp(),  tkrParTds.getySlpySlp() );
    AcdId acdId = convert(tds.acdId());
    root.initialize( acdId, tds.trackIndex(),
                     tds.doca(),tds.docaErr(),tds.docaRegion(),
                     poca,paramsAtPoca);
  }

  void convert( const AcdTkrPoca& root, Event::AcdTkrPoca& tds) {
    static Point poca;
    const TVector3& rootPoca = root.getPoca();
    poca.set(rootPoca.X(),rootPoca.Y(),rootPoca.Z());        
    idents::AcdId acdIdTds = convert(root.getId());
    const TkrTrackParams& paramsRoot = root.getParamsAtPoca();
    Event::TkrTrackParams paramsAtPoca(paramsRoot.getxPosition(), paramsRoot.getxSlope(), 
                                       paramsRoot.getyPosition(), paramsRoot.getySlope(),
                                       paramsRoot.getxPosxPos(),  paramsRoot.getxPosxSlp(), 
                                       paramsRoot.getxPosyPos(),  paramsRoot.getxPosySlp(),
                                       paramsRoot.getxSlpxSlp(),  paramsRoot.getxSlpyPos(), 
                                       paramsRoot.getxSlpySlp(),  paramsRoot.getyPosyPos(), 
                                       paramsRoot.getyPosySlp(),  paramsRoot.getySlpySlp() );
    tds.set(acdIdTds,root.getTkrIndex(),
            root.getDoca(),root.getDocaErr(),root.getDocaRegion(),
            poca,paramsAtPoca);
  }

  void convert( const Event::AcdTkrPoint& tds, AcdTkrPoint& root) {
    static TVector3 point;
    const Event::TkrTrackParams& tkrParTds = tds.paramsAtPoint();
    const Point& tdsPoint = tds.point();
    TkrTrackParams params(tkrParTds.getxPosition(), tkrParTds.getxSlope(), 
                          tkrParTds.getyPosition(), tkrParTds.getySlope(),
                          tkrParTds.getxPosxPos(),  tkrParTds.getxPosxSlp(), 
                          tkrParTds.getxPosyPos(),  tkrParTds.getxPosySlp(),
                          tkrParTds.getxSlpxSlp(),  tkrParTds.getxSlpyPos(), 
                          tkrParTds.getxSlpySlp(),  tkrParTds.getyPosyPos(), 
                          tkrParTds.getyPosySlp(),  tkrParTds.getySlpySlp() );
    point.SetXYZ(tdsPoint.x(),tdsPoint.y(),tdsPoint.z());
    root.set(tds.arcLength(),tds.face(),point,params);
  }
  
  void convert( const AcdTkrPoint& root, Event::AcdTkrPoint& tds) {
    static Point point;
    const TkrTrackParams& paramsRoot = root.getParamsAtPoint();
    const TVector3& rootPoint = root.getPoint();
    point.set(rootPoint.X(),rootPoint.Y(),rootPoint.Z());
    Event::TkrTrackParams params(paramsRoot.getxPosition(), paramsRoot.getxSlope(), 
                                 paramsRoot.getyPosition(), paramsRoot.getySlope(),
                                 paramsRoot.getxPosxPos(),  paramsRoot.getxPosxSlp(), 
                                 paramsRoot.getxPosyPos(),  paramsRoot.getxPosySlp(),
                                 paramsRoot.getxSlpxSlp(),  paramsRoot.getxSlpyPos(), 
                                 paramsRoot.getxSlpySlp(),  paramsRoot.getyPosyPos(), 
                                 paramsRoot.getyPosySlp(),  paramsRoot.getySlpySlp() );
    tds.setLocalData(root.getArcLength(),root.getFace(),point,params);
  }

  void convert( const Event::AcdTkrPoint& tds, AcdTkrPointV2& root) {
    root.set(tds.getTrackIndex());
    convert( (const Event::AcdTkrLocalCoords&)tds, root ); 
  }
  
  void convert( const AcdTkrPointV2& root, Event::AcdTkrPoint& tds) {
    tds.set(root.getTrackIndex());
    convert( (const AcdTkrLocalCoordsV2&)root, tds ); 
  }

  void convert( const Event::AcdTkrHitPoca& tds, AcdTkrHitPoca& root) {
    static TMatrixD localCov(2,2);
    static float position[2];
    static AcdTkrLocalCoords local;
    static TVector3 poca;
    static TVector3 pocaVector;
    static AcdPocaData pocaData;

    position[0] = tds.getActiveX();
    position[1] = tds.getActiveY();
    localCov[0][0] = tds.getLocalXXCov();
    localCov[1][1] = tds.getLocalYYCov();
    localCov[0][1] = localCov[1][0] =  tds.getLocalXYCov();    
    local.set(position,tds.getPathLength(),tds.getCosTheta(),tds.getRegion(),localCov);
    poca.SetXYZ(tds.getPoca().x(),tds.getPoca().y(),tds.getPoca().z());
    pocaVector.SetXYZ(tds.getPocaVector().x(),tds.getPocaVector().y(),tds.getPocaVector().z());
    pocaData.set(tds.getArcLength(),tds.getDoca(),tds.getDocaErr(),poca,pocaVector);         
    AcdId acdId = convert(tds.getId());
    root.set(acdId,tds.trackIndex(),local,pocaData);
  }

  void convert( const AcdTkrHitPoca& root, Event::AcdTkrHitPoca& tds) {
    static HepSymMatrix localCov(2);
    static float position[2];
    static const float mips[2] = {0,0};
    static const unsigned short flags[2] = {0,0};
    static Point poca;
    static Vector pocaVector;
    position[0] = root.getActiveX();
    position[1] = root.getActiveY();
    localCov[0][0] = root.getLocalXXCov();
    localCov[1][1] = root.getLocalYYCov();
    localCov[0][1] = localCov[1][0] =  root.getLocalXYCov();    
    poca.set(root.getPoca().X(),root.getPoca().Y(),root.getPoca().Z());
    pocaVector.set(root.getPocaVector().X(),root.getPocaVector().Y(),root.getPocaVector().Z());
    idents::AcdId acdId = convert(root.getId());
    // set Top part of TDS object
    tds.set(acdId,root.getTrackIndex(),mips,0.,0.,0.,flags);
    // set PocaData part of TDS object
    tds.setPocaData(root.getArcLength(),root.getDoca(),root.getDocaErr(),poca,pocaVector);      
    // set LocalData part of TDS object
    tds.setLocalData(position,root.getPathLength(),root.getCosTheta(),root.getRegion(),localCov);
  }

  void convert( const Event::AcdTkrHitPoca& tds, AcdTkrHitPocaV2& root) {
    static Float_t mips[2];
    static UShort_t flags[2];
    AcdId acdId = convert(tds.getId());
    mips[0] = tds.mipsPmtA(); mips[1] = tds.mipsPmtB();     
    flags[0] = tds.flagsPmtA(); flags[1] = tds.flagsPmtB();     
    root.set(acdId,tds.trackIndex(),mips,
             tds.vetoSigmaHit(),tds.vetoSigmaProj(),tds.vetoSigmaProp(),flags ); 
    convert( (const Event::AcdTkrLocalCoords&)tds, root );
    convert( (const Event::AcdPocaData&)tds, root );  
  }
  void convert( const AcdTkrHitPocaV2& root, Event::AcdTkrHitPoca& tds) {
    static Float_t mips[2];
    static UShort_t flags[2];
    idents::AcdId acdId = convert(root.getId());
    mips[0] = root.mipsPmtA(); mips[1] = root.mipsPmtB();     
    flags[0] = root.flagsPmtA(); flags[1] = root.flagsPmtB();     
    tds.set(acdId,root.getTrackIndex(),mips,
            root.vetoSigmaHit(),root.vetoSigmaProj(),root.vetoSigmaProp(),flags); 
    convert( (const AcdTkrLocalCoordsV2&)root, tds );
    convert( (const AcdPocaDataV2&)root, tds );
  }

  void convert( const Event::AcdTkrGapPoca& tds, AcdTkrGapPoca& root) {
    static TMatrixD localCov(2,2);
    static float position[2];
    static AcdTkrLocalCoords local;
    static TVector3 poca;
    static TVector3 pocaVector;
    static AcdPocaData pocaData;
    position[0] = tds.getActiveX();
    position[1] = tds.getActiveY();
    localCov[0][0] = tds.getLocalXXCov();
    localCov[1][1] = tds.getLocalYYCov();
    localCov[0][1] = localCov[1][0] =  tds.getLocalXYCov();    
    local.set(position,tds.getPathLength(),tds.getCosTheta(),tds.getRegion(),localCov);
    AcdGapId gapId( tds.getId().asShort() );
    poca.SetXYZ(tds.getPoca().x(),tds.getPoca().y(),tds.getPoca().z());
    pocaVector.SetXYZ(tds.getPocaVector().x(),tds.getPocaVector().y(),tds.getPocaVector().z());
    pocaData.set(tds.getArcLength(),tds.getDoca(),tds.getDocaErr(),poca,pocaVector);         
    root.initialize(gapId,tds.trackIndex(),local,pocaData); 
  }

  void convert( const AcdTkrGapPoca& root, Event::AcdTkrGapPoca& tds) {
    static HepSymMatrix localCov(2);
    static float position[2];
    static Point poca;
    static Vector pocaVector;
    idents::AcdGapId acdId( root.getGapId().asShort() );
    position[0] = root.getActiveX();
    position[1] = root.getActiveY();
    localCov[0][0] = root.getLocalXXCov();
    localCov[1][1] = root.getLocalYYCov();
    localCov[0][1] = localCov[1][0] =  root.getLocalXYCov();    
    poca.set(root.getPoca().X(),root.getPoca().Y(),root.getPoca().Z());
    pocaVector.set(root.getPocaVector().X(),root.getPocaVector().Y(),root.getPocaVector().Z());
    // set Top part of TDS object
    tds.set(acdId,root.getTrackIndex(),0.,0.,0.); 
    // set PocaData part of TDS object    
    tds.setPocaData(root.getArcLength(),root.getDoca(),root.getDocaErr(),poca,pocaVector);      
    // set LocalData part of TDS object
    tds.setLocalData(position,root.getPathLength(),root.getCosTheta(),root.getRegion(),localCov);    
  }

  void convert( const Event::AcdTkrGapPoca& tds, AcdTkrGapPocaV2& root) {
    AcdGapId gapId( tds.getId().asShort() );
    root.initialize(gapId,tds.trackIndex(),
                    tds.vetoSigmaHit(),tds.vetoSigmaProj(),tds.vetoSigmaProp());
    convert( (const Event::AcdTkrLocalCoords&)tds, root );
    convert( (const Event::AcdPocaData&)tds, root );
  }

  void convert( const AcdTkrGapPocaV2& root, Event::AcdTkrGapPoca& tds) {    
    idents::AcdGapId gapId( root.getGapId().asShort() );
    tds.set(gapId,root.getTrackIndex(),
            root.vetoSigmaHit(),root.vetoSigmaProj(),root.vetoSigmaProp() ); 
    convert( (const AcdTkrLocalCoordsV2&)root, tds );
    convert( (const AcdPocaDataV2&)root, tds );    
  }

  void convert( const Event::AcdSplashVars& tds, AcdSplashVars& root) {
    static TVector3 point;
    static TVector3 vector;
    point.SetXYZ(tds.calEntryPoint().x(),tds.calEntryPoint().y(),tds.calEntryPoint().z());
    vector.SetXYZ(tds.calEntryVector().x(),tds.calEntryVector().y(),tds.calEntryVector().z());
    AcdId acdId = convert(tds.getId());
    root.set(acdId,tds.getTrackIndex(),point,vector,
             tds.tileSolidAngle(),tds.weightedTrackAngle(),tds.weightedPathlength());
  }

  void convert( const AcdSplashVars& root, Event::AcdSplashVars& tds) {
    static Point point;
    static Vector vector;
    point.set(root.calEntryPoint().X(),root.calEntryPoint().Y(),root.calEntryPoint().Z());
    vector.set(root.calEntryVector().X(),root.calEntryVector().Y(),root.calEntryVector().Z());
    idents::AcdId acdId = convert(root.getId());
    tds.set(acdId,root.getTrackIndex(),point,vector,
            root.tileSolidAngle(),root.weightedTrackAngle(),root.weightedPathlength());
  }
  


  void convert( const Event::AcdTkrLocalCoords& tds, AcdTkrLocalCoordsV2& root) {
    static Float_t active2d[2];
    static Float_t local[2];
    static TMatrixFSym localProj(2);
    static TMatrixFSym localProp(2);
    static TVector3 point;
    local[0] = tds.getLocalX();
    local[1] = tds.getLocalY();    
    active2d[0] = tds.getActiveX(); 
    active2d[1] = tds.getActiveY(); 
    const HepPoint3D& global = tds.getGlobalPosition();
    point.SetXYZ(global.x(),global.y(),global.z());
    convert( tds.getLocalCovProj(), localProj);
    convert( tds.getLocalCovProp(), localProp);
    root.set(tds.getLocalVolume(),tds.getArclengthToPlane(),tds.getCosTheta(),
             point,local,active2d,
             localProj,localProp);             
  }

  void convert( const AcdTkrLocalCoordsV2& root, Event::AcdTkrLocalCoords& tds) {
    static Float_t active2d[2];
    static Float_t local[2];
    static HepSymMatrix localProj(2,0);
    static HepSymMatrix localProp(2,0);
    static HepPoint3D point;
    local[0] = root.getLocalX();
    local[1] = root.getLocalY();    
    active2d[0] = root.getActiveX(); 
    active2d[1] = root.getActiveY(); 
    const TVector3& global = root.getGlobalPosition();
    point.set(global.X(),global.Y(),global.Z());
    convert( root.getLocalCovProj(), localProj);
    convert( root.getLocalCovProp(), localProp);
    tds.setLocalData(root.getVolume(),root.getArclengthToPlane(),root.getCosTheta(),
                     point,local,active2d,
                     localProj,localProp);
  }

  void convert( const Event::AcdPocaData& tds, AcdPocaDataV2& root) {
    static TVector3 poca;
    static TVector3 pocaVector;
    poca.SetXYZ(tds.getPoca().x(),tds.getPoca().y(),tds.getPoca().z());
    pocaVector.SetXYZ(tds.getPocaVector().x(),tds.getPocaVector().y(),tds.getPocaVector().z());
    root.set( tds.getVolume(), tds.getRegion(), tds.getArcLength(),
              tds.getDoca(), tds.getDocaErrProj(), tds.getDocaErrProp(), 
              poca, pocaVector );
  }

  void convert( const AcdPocaDataV2& root, Event::AcdPocaData& tds) {
    static Point poca;
    static Vector pocaVector;
    poca.set(root.getPoca().X(),root.getPoca().Y(),root.getPoca().Z());
    pocaVector.set(root.getPocaVector().X(),root.getPocaVector().Y(),root.getPocaVector().Z());
    tds.setPocaData( root.getVolume(), root.getRegion(), root.getArcLength(),
                     root.getDoca(), root.getDocaErrProj(), root.getDocaErrProp(), 
                     poca, pocaVector );
  }

  void convert( const HepSymMatrix& tds, TMatrixDSym& root) {
    //ADW: Make this function do something...
    //return;
    if ( tds.num_row() != root.GetNrows() ) return;
    int size = tds.num_row();
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        root[i][j] = root[j][i] = tds[i][j];
      }
    }
    return;
  }

  void convert( const TMatrixDSym& root, HepSymMatrix& tds) {    
    //ADW: Make this function do something...
    //return;
    if ( root.GetNrows() != tds.num_row() ) return;
    Int_t size = root.GetNcols();
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        tds[i][j] = tds[j][i] = root[i][j];
      }
    }
    return;
  }

  void convert( const HepSymMatrix& tds, TMatrixFSym& root) {
    if ( tds.num_row() != root.GetNrows() ) return;
    int size = tds.num_row();
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        root[i][j] = root[j][i] = tds[i][j];
      }
    }
    return;
  }

  void convert( const TMatrixFSym& root, HepSymMatrix& tds) {    
    if ( root.GetNrows() != tds.num_row() ) return;
    Int_t size = root.GetNcols();
    for (int i = 0; i < size; i++) {
      for (int j = i; j < size; j++) {
        tds[i][j] = tds[j][i] = root[i][j];
      }
    }
    return;
  }

}



