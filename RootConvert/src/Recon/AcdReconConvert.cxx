
#include <RootConvert/Recon/AcdReconConvert.h>
#include <RootConvert/Utilities/Toolkit.h>
#include <CLHEP/Matrix/Matrix.h>
#include <TMatrixD.h>

namespace RootPersistence {

void convert( const Event::AcdRecon & tdsAcdRec, AcdRecon & rootAcdRec ) {
    
  std::vector<AcdId> rootIdCol ;
  std::vector<idents::AcdId>::const_iterator tdsIdIt ;
  for ( tdsIdIt = tdsAcdRec.getIdCol().begin() ;
        tdsIdIt != tdsAcdRec.getIdCol().end() ; tdsIdIt++ ) {
    rootIdCol.push_back(
        AcdId(
            tdsIdIt->layer(), tdsIdIt->face(),
            tdsIdIt->row(), tdsIdIt->column()
        )
    ) ;
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
    
  const Event::AcdTkrIntersectionCol& tdsAcdTkrIntersects
   = tdsAcdRec.getAcdTkrIntersectionCol() ;

  TVector3 globalPosition;
  Double_t localPosition[2];
  TMatrixD covMatrix(2,2);

  for ( Event::AcdTkrIntersectionCol::const_iterator itr = tdsAcdTkrIntersects.begin() ;
        itr != tdsAcdTkrIntersects.end(); itr++ ) {
    const Event::AcdTkrIntersection* tdsAcdIntersect = *itr;
    const Point& glbPos = tdsAcdIntersect->getGlobalPosition();
    globalPosition.SetXYZ(glbPos.x(),glbPos.y(),glbPos.z());
    localPosition[0] = tdsAcdIntersect->getLocalX();
    localPosition[1] = tdsAcdIntersect->getLocalY();
    covMatrix[0][0] = tdsAcdIntersect->getLocalXXCov();
    covMatrix[1][1] = tdsAcdIntersect->getLocalYYCov();
    covMatrix[0][1] = tdsAcdIntersect->getLocalXYCov();
    covMatrix[1][0] = tdsAcdIntersect->getLocalXYCov();
    AcdId acdId( tdsAcdIntersect->getTileId().id() );
    AcdTkrIntersection rootAcdInter( acdId, tdsAcdIntersect->getTrackIndex(),
      globalPosition, 
      localPosition, covMatrix, 
      tdsAcdIntersect->getArcLengthToIntersection(),
      tdsAcdIntersect->getPathLengthInTile(),
      tdsAcdIntersect->tileHit());
    rootAcdRec.addAcdTkrIntersection(rootAcdInter);
  }

}
 
void convert( const AcdRecon & rootAcdRec, Event::AcdRecon & tdsAcdRec )
{
    std::vector<Event::AcdTkrIntersection*> acdTkrIntersections;
    Point globalPosition;
    HepMatrix covMatrix(2,2);
    double localPosition[2];
    int nInter = rootAcdRec.nAcdIntersections() ;
    for ( int iInter(0); iInter < nInter; iInter++ ) {

        const AcdTkrIntersection * rootAcdInter = rootAcdRec.getAcdTkrIntersection(iInter);      

        const TVector3& glbPos = rootAcdInter->getGlobalPosition();
        globalPosition.set(glbPos.X(),glbPos.Y(),glbPos.Z());
        localPosition[0] = rootAcdInter->getLocalX();
        localPosition[1] = rootAcdInter->getLocalY();
        covMatrix[0][0] = rootAcdInter->getLocalXXCov();
        covMatrix[1][1] = rootAcdInter->getLocalYYCov();
        covMatrix[0][1] = covMatrix[1][0] = rootAcdInter->getLocalXYCov();
        Event::AcdTkrIntersection* tdsAcdInter = new
            Event::AcdTkrIntersection(
                rootAcdInter->getTileId().getId(), rootAcdInter->getTrackIndex(),
                globalPosition,localPosition, covMatrix,
                rootAcdInter->getArcLengthToIntersection(),
                rootAcdInter->getPathLengthInTile(),
                rootAcdInter->tileHit()
                ) ;
        acdTkrIntersections.push_back(tdsAcdInter) ;
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
        rootAcdRec.getRibbonActiveDist(),
        tdsRibActDistId, 
        rootAcdRec.getCornerDoca()
        ) ;
}

}
