/**
 * @class TkrDigiHitToDigiAlg
 *
 * @brief Calls an user-chosen tool to convert tkr hits into tkr digis.
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#ifndef __TKRDIGIHITTODIGIALG_H__
#define __TKRDIGIHITTODIGIALG_H__

#include "../IHitToDigiTool.h"

#include "GaudiKernel/Algorithm.h"

class TkrDigiHitToDigiAlg : public Algorithm {

 public:

    TkrDigiHitToDigiAlg(const std::string&, ISvcLocator*);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

 private:

    /// Type of tool to run
    std::string m_type;
    /// Pointer to the tool
    IHitToDigiTool* m_tool;

};

#endif
