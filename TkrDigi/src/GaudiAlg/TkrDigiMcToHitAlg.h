/*
 * @class TkrDigiMcToHitAlg
 *
 * @brief Calls either Simple or Bari tool to convert MC hits into tkr hits.
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#ifndef __TKRDIGIMCTOHITALG_H__
#define __TKRDIGIMCTOHITALG_H__

#include "../IMcToHitTool.h"

#include "GaudiKernel/Algorithm.h"

class TkrDigiMcToHitAlg : public Algorithm {

 public:

    TkrDigiMcToHitAlg(const std::string&, ISvcLocator*);
    StatusCode initialize();
    StatusCode execute();
    StatusCode finalize();

 private:

    /// Type of tool to run
    std::string m_type;
    /// Pointer to the tool
    IMcToHitTool* m_tool;

};

#endif
