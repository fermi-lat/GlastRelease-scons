/*
 * @file TkrDigi_load.cxx
 *
 * @brief This is needed for forcing the linker to load all components
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(TkrDigi)
{
  DECLARE_ALGORITHM(TkrDigiAlg);
  DECLARE_ALGORITHM(TkrDigiMcToHitAlg);
  DECLARE_ALGORITHM(TkrDigiNoiseAlg);
  DECLARE_ALGORITHM(TkrDigiHitToDigiAlg);
  DECLARE_TOOL     (BariMcToHitTool);
  DECLARE_TOOL     (SimpleMcToHitTool);
  DECLARE_TOOL     (GeneralNoiseTool);
  DECLARE_TOOL     (GeneralHitToDigiTool);
  DECLARE_TOOL     (TkrDigiRandom); 
}
