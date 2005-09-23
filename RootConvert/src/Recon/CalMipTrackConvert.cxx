
#include <RootConvert/Recon/CalMipTrackConvert.h>
#include <RootConvert/Recon/CalParamsConvert.h>
#include <RootConvert/Utilities/Toolkit.h>

namespace RootPersistence {

void convert( const Event::CalMipTrack& calMipTrackTds, CalMipTrack& calMipTrackRoot )
{
    Point    pointTds = calMipTrackTds.getPoint();
    Vector   dirTds   = calMipTrackTds.getDir();

    TVector3 pointRoot(pointTds.x(), pointTds.y(), pointTds.z());
    TVector3 dirRoot(dirTds.x(), dirTds.y(), dirTds.z());

    calMipTrackRoot.setPoint(pointRoot);
    calMipTrackRoot.setDir(dirRoot);

    double   d2CTds     = calMipTrackTds.getD2C();
    double   d2EdgeTds  = calMipTrackTds.getD2Edge();
    int      calEdgeTds = calMipTrackTds.getCalEdge();
    double   arcLenTds  = calMipTrackTds.getArcLen();
    double   ecorTds    = calMipTrackTds.getEcor();
    double   ecorRmsTds = calMipTrackTds.getEcorRms();
    double   chi2Tds    = calMipTrackTds.getChi2();   
    double   ermTds     = calMipTrackTds.getErm();   
    
    calMipTrackRoot.setD2C(d2CTds);    
    calMipTrackRoot.setD2Edge(d2EdgeTds); 
    calMipTrackRoot.setCalEdge(calEdgeTds);
    calMipTrackRoot.setArcLen(arcLenTds); 
    calMipTrackRoot.setEcor(ecorTds);   
    calMipTrackRoot.setEcorRms(ecorRmsTds);
    calMipTrackRoot.setChi2(chi2Tds);   
    calMipTrackRoot.setErm(ermTds);   
}
 
void convert( const CalMipTrack& calMipTrackRoot, Event::CalMipTrack& calMipTrackTds )
{
    TVector3 pointRoot = calMipTrackRoot.getPoint();
    TVector3 dirRoot   = calMipTrackRoot.getDir();

    Point    pointTds(pointRoot.x(), pointRoot.y(), pointRoot.z());
    Vector   dirTds(dirRoot.x(), dirRoot.y(), dirRoot.z());

    calMipTrackTds.setPoint(pointTds);
    calMipTrackTds.setDir(dirTds);

    double   d2CRoot     = calMipTrackRoot.getD2C();
    double   d2EdgeRoot  = calMipTrackRoot.getD2Edge();
    int      calEdgeRoot = calMipTrackRoot.getCalEdge();
    double   arcLenRoot  = calMipTrackRoot.getArcLen();
    double   ecorRoot    = calMipTrackRoot.getEcor();
    double   ecorRmsRoot = calMipTrackRoot.getEcorRms();
    double   chi2Root    = calMipTrackRoot.getChi2();   
    double   ermRoot     = calMipTrackRoot.getErm();   
    
    calMipTrackTds.setD2C(d2CRoot);    
    calMipTrackTds.setD2Edge(d2EdgeRoot); 
    calMipTrackTds.setCalEdge(calEdgeRoot);
    calMipTrackTds.setArcLen(arcLenRoot); 
    calMipTrackTds.setEcor(ecorRoot);   
    calMipTrackTds.setEcorRms(ecorRmsRoot);
    calMipTrackTds.setChi2(chi2Root);   
    calMipTrackTds.setErm(ermRoot);   
}

}
