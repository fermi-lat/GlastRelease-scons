//      $Header$
//
// Description:
//      TkrMakeClusterTable has the methods for making the clusters, 
//      and for setting the cluster flag.
//
// Author(s):
//      Tracy Usher     
//      Leon Rochester     



#include "src/Cluster/TkrMakeClusterTable.h"
#include "facilities/Util.h"
#include <algorithm>

using namespace Event;

TkrMakeClusterTable::TkrMakeClusterTable(const TkrClusterCol* pClus, 
                                         const TkrDigiCol* pDigi,
                                         ObjectList<
                                         Relation<TkrDigi,McPositionHit> >* pRelTab,
                                         RelTable<TkrCluster, McPositionHit>* pclustTab,
                                         ITkrGeometrySvc* tkrGeom)
{   
    m_tkrGeom = tkrGeom;

    RelTable <TkrDigi,McPositionHit>     digiHitsTab(pRelTab);

    RelTable <TkrCluster, McPositionHit> clustHitsTab = *pclustTab;

    TkrDigiCol::const_iterator itD = pDigi->begin();

    int clsSize = pClus->size();

    // go through all the clusters: these are in the same order as the digis
    for (int iclu=0; iclu<clsSize; iclu++) {
        TkrCluster* p_clu = (*pClus)[iclu];
        int firstStrip = p_clu->firstStrip();
        int lastStrip  = p_clu->lastStrip();
        // turn into strings for comparison
        std::string lowStrip;
        std::string highStrip;
        facilities::Util::itoa(firstStrip, lowStrip);
        while(lowStrip.size()<4) { lowStrip = " "+lowStrip;}
        facilities::Util::itoa(lastStrip, highStrip);
        while(highStrip.size()<4) { highStrip = " "+highStrip;}
        TkrDigi* p_digi = *itD;
        int order = digiOrder(p_clu);
             //std::cout << "Before loop (o/d) " << order << " " 
             //   << (int) *p_digi << std::endl;
        while (order!=*p_digi && itD!=pDigi->end()) {
             //std::cout << "In     loop (o/d) " << order << " " 
             //   << (int) *p_digi << std::endl;
            p_digi = *(++itD);
        }
             //std::cout << "After  loop (o/d) " << order << " " 
             //   << (int) *p_digi << std::endl;
        if (itD==pDigi->end()) return;
        // cluster and digi match; get the McHits
        std::vector<Relation<TkrDigi, McPositionHit> *> relsByDigi = digiHitsTab.getRelByFirst(p_digi);
        // collect the corresponding hits
        std::vector<McPositionHit*> mcHits;
        int relSize = relsByDigi.size();
        for(int irel=0;irel<relSize;irel++) {
            McPositionHit* theHit = relsByDigi[irel]->getSecond();
            // is this hit already in the list?
            if(std::find(mcHits.begin(),mcHits.end(), theHit)==mcHits.end()) {
                // if not, does the info contain any strip in the cluster?
                
                std::vector<std::string> info = relsByDigi[irel]->getInfos();
                std::vector<std::string>::const_iterator itI;
                
                for (itI=info.begin(); itI!=info.end(); itI++) {
                    std::string thisStrip= *itI;
                    if(thisStrip<lowStrip)  {continue;}
                    if(thisStrip>highStrip) {continue;}
                    // in range, add the hit, and get out
                    mcHits.push_back(theHit);
                    break;
                }
            }
        }

        // make a relation with the current cluster and each McHit
        std::vector<McPositionHit*>::const_iterator itH;
        for(itH=mcHits.begin(); itH!=mcHits.end();itH++) {
            McPositionHit* theHit = *itH;
            Relation<TkrCluster, McPositionHit>* rel 
                = new Relation<TkrCluster, McPositionHit>(p_clu, theHit);
            clustHitsTab.addRelation(rel);
        }
    }   
}

int TkrMakeClusterTable::digiOrder ( const TkrCluster* pClust) {
//    TkrCluster clust = *pClust;
//    return clust.v() + 2*(m_tkrGeom->reverseLayerNumber(clust.plane()))
//        + 64*clust.tower();
    return     pClust->getTkrId().getView() 
           + 2*pClust->getTkrId().getLayer() + 64*pClust->tower();
}
    

