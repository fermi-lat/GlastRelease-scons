#include "OnboardFilterTds/Obf_TFC_prjs.h"
#ifdef OBF_B1_3_3
# include "EDS/FFS.h"
#else 
# include "PBI/FFS.ih"
# define  FFS_mask  FFSL_mask
# define  FFS_eliminate FFSL_eliminate
# define  FFS  FFSL
#endif

#include <iostream>
#include <cstring>

namespace OnboardFilterTds{
    
Obf_TFC_prjs::Obf_TFC_prjs(const TFC_prjs* prjs)
{
    // How much memory to delete
    static int memSize = sizeof(TFC_prjs) - 1000*sizeof(TFC_prj) + sizeof(TFC_prj*);

    // Zero the memory locations
    memset(&m_maxCnt, 0, memSize);

    // Now fill from the input structure
    initialize(prjs);
}

Obf_TFC_prjs::~Obf_TFC_prjs()
{
    delete [] m_prjs;
}

void Obf_TFC_prjs::initialize(const TFC_prjs* prjs)
{
    // Copy basic information 
    m_maxCnt          = prjs->maxCnt;
    m_curCnt          = prjs->curCnt;
    m_twrMsk          = prjs->twrMsk;

    // Recover the tower mask
    unsigned int twrMask = m_twrMsk << 16;

    // Loop over towers with projects copying the prjDir info
    while(twrMask)
    {
        int towerId = FFS (twrMask);

        m_dir[towerId]          = prjs->dir[towerId];
      
        twrMask = FFS_eliminate(twrMask, towerId);
    }

    // DON'T copy the top/bot prjlink since these pointers are mostly invalid

    // Allocate memory for the projections
    m_prjs = new TFC_prj[prjs->curCnt];

    // Copy the projections in one operation...
    memcpy(m_prjs, prjs->prjs, prjs->curCnt * sizeof(TFC_prj));

    // Finally, clean out the same in the projections
    for (int idx = 0; idx < m_curCnt; idx++)
    {
        TFC_prj& prj = m_prjs[idx];

        prj.topNode.flnk = 0;
        prj.topNode.blnk = 0;
        prj.botNode.flnk = 0;
        prj.botNode.blnk = 0;
    }
}

} // namespace
