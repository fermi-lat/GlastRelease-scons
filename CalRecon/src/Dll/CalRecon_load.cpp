/** 
* @file CalRecon_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*/

#include "GaudiKernel/DeclareFactoryEntries.h"
//#include "GaudiKernel/IToolFactory.h"

DECLARE_FACTORY_ENTRIES( CalRecon ) {
    DECLARE_SERVICE( CalReconSvc ) ;
    DECLARE_ALGORITHM( CalClustersAlg ) ;
    DECLARE_ALGORITHM( CalEventEnergyAlg ) ;
    DECLARE_ALGORITHM( CalDisplay ) ;
    DECLARE_ALGORITHM( PropertiesCheckAlg ) ;
    DECLARE_TOOL( CalSingleClustering ) ;
    DECLARE_TOOL( CalSimpleClustering ) ;
    DECLARE_TOOL( CalRawEnergyCorr );
    DECLARE_TOOL( CalLastLayerCorr );
    DECLARE_TOOL( CalTkrLikelihoodCorr );
    DECLARE_TOOL( CalProfileCorr );
    DECLARE_TOOL( CalValsCorr );
    DECLARE_TOOL( CalTransvOffsetCorr );
} 




