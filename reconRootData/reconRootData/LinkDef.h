/** @file LinkDef.h
    @brief for rootcint
 $Header$

*/
#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace reconRootData;
#pragma link C++ class ReconEvent+;

#pragma link C++ class AcdTkrIntersection+;
#pragma link C++ class AcdTkrPoca+;
#pragma link C++ class AcdHit+;
#pragma link C++ class AcdPocaData+;
#pragma link C++ class AcdTkrPoint+;
#pragma link C++ class AcdTkrLocalCoords+;
#pragma link C++ class AcdTkrHitPoca+;
#pragma link C++ class AcdTkrGapPoca+;
#pragma link C++ class AcdSplashVars+;

#pragma link C++ class AcdPocaDataV2+;
#pragma link C++ class AcdTkrLocalCoordsV2+;
#pragma link C++ class AcdTkrPointV2+;
#pragma link C++ class AcdTkrHitPocaV2+;
#pragma link C++ class AcdTkrGapPocaV2+;
#pragma link C++ class AcdAssoc+;
#pragma link C++ class AcdEventTopology+;

#pragma link C++ class AcdRecon+;
#pragma link C++ class AcdReconV2+;

#pragma link C++ class CalRecon+;
#pragma link C++ class CalXtalsParams+;
#pragma link C++ class CalParams+;
#pragma link C++ class CalMomParams+;
#pragma link C++ class CalFitParams+;
#pragma link C++ class CalMSTreeParams+;
#pragma link C++ class CalClassParams+;
#pragma link C++ class CalClusterLayerData+;
#pragma link C++ class CalCluster+;
#pragma link C++ class CalCorToolResult+;
#pragma link C++ class CalEventEnergy+;
#pragma link C++ class CalRangeRecData+;
#pragma link C++ class CalXtalRecData+;
#pragma link C++ class CalMipTrack+;
//CL 08/22/06:
#pragma link C++ class GcrXtal+;
#pragma link C++ class GcrTrack+;

#pragma link C++ class TkrRecon+;
#pragma link C++ class TkrCluster+;
#pragma link C++ class TkrTrack+;
#pragma link C++ class TkrTrackHit+;
#pragma link C++ class TkrTrackParams+;
#pragma link C++ class TkrVertex+;
#pragma link C++ class TkrVecPoint+;
#pragma link C++ class TkrVecPointsLink+;
#pragma link C++ class TkrVecNode+;
#pragma link C++ class TkrTree+;
#pragma link C++ class TkrFilterParams+;
#pragma link C++ class TkrTruncationData+;
#pragma link C++ class TkrEventParams+;
#pragma link C++ class TreeClusterRelation+;
#pragma link C++ class TkrDiagnostics+;

#pragma link C++ class reconRootData::TaggerCluster+;
#pragma link C++ class reconRootData::AdfRecon+;

#endif
