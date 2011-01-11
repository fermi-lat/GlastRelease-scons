
// Mainpage for doxygen

/** @mainpage package CalRecon
 *
 * @authors A.Chekhtman, R.Terrier, J.A.Hernando, G. Musat, B. Giebels, L. Baldini, J. Bregeon, C. Sgro'
 *
 * @section intro Introduction
 *
 *  CalRecon package reconstructs the energy and direction of incident particle from
 *   the calorimeter information. 
 *   
 *   The package contains 6 algorithms:
 *   - CalClustersAlg: groups nearby hit cristals into clusters, different algorithms are available
 *   - CalClassifyAlg: classifies clusters according to their topology (gamma-like, mip-like...)
 *   - CalEventEnergyAlg: controls and applies the various energy correction tools used to determine the final event energy
 *   - CalMipFinderAlg: looks for minimum ionizing particle tracks in the calorimeter (usually turned off)
 *   - PropertiesCheckAlg: internal algorithm to check jobOption properties ?
 *   - CalDisplay: extracts the reconstructed calorimeter data from TDS classes and draw the graphic representation
 *
 *   The  control flow diagram of CalRecon package, including also the CalDigi package,
 *   is given on the following WEB page:
 * <A HREF="http://www-glast.slac.stanford.edu/software/CAL/meetings/calrecondiagram.pdf"> CalRecon diagram </A> 
 *
 * @section algs Description of the main algorithms
 *   CalDisplay provides the display of reconstructed data.
 *
 *   CalClustersAlg calculates the energy, position and direction for
 *   calorimeter clusters and applies energy corrections.
 *   Actually there is no real clustering algorithm implemented,
 *   CalClusterAlg now considers, that there is only one cluster including
 *   all hitted calorimeter crystals. For this "cluster" the algorithm
 *   calculates the following parameters and stores them in the 
 *   Event::CalClusterCol object:
 *     - total energy sum 
 *     - energy corrected by profile fitting method
 *     - energy corrected by last layer correlation method
 *     - energy sum for each layer
 *     - average position 
 *     - average position for each layer
 *     - direction of incident particle
 *    
 *   CalClassifyAlg
 *
 *    
 *   CalEventEnergyAlg
 *
 *    
 *   CalMipFinderAlg
 *
 *
 * @section tools Tools
 *
 * The original clustering and energy corrections (profile and last-layer
 * correlation methods) have been recast as Gaudi tools.
 *
 * The clustering is now based on the ICalClustering and Clustering base classes, 
 * while the leakage corrections derive from ICalEnergyCorr and AlgTool.
 * The clustering tools currently available are CalSingleClustering and
 * CalSimpleClustering (FuzzyClustering needs additional tuning). The last
 * layer leakage tool is CalLastLayerLikelihoodTool. The tracker hit correlation
 * tool is CalTkrLikelihoodTool. The last two tools derive from CaLikelihood.
 * The profile tool is CalProfileTool.
 *
 * The calo FuzzyClustering retrieves the generic FuzzyCluster Gaudi
 * tool to perform fuzzy clustering on all cal hits. However, if there are
 * less than one hit point per cluster, it applies one single cluster
 * calculations, like the CalSingleClustering. 
 *
 * CalClustersAlg calls all 3 tools so far.
 *
 * <hr>
 *
 * @section jobOptions jobOptions
 * @param CalClustersAlg.callNumber
 *        this parameter is used to distinguish multiple calls to 
 *        CalClustersAlg (for example, before and after TkrRecon).
 *        The default value is 0 .
 *
 * @param CalClustersAlg.clusteringToolName
 *        name of tool performing clustering. Default is CalSingleClustering
 *        If set to FuzzyClustering, the following param is mandatory:
 * @param ToolSvc.FuzzyClusterTool.FuzzyCluster.command
 *        the command param is mandatory for the FuzzyCluster Gaudi tool
 *        which is retrieved as a private tool by FuzzyClusterTool, that
 *        in its turn is a shared tool. Example:
 *        ToolSvc.FuzzyClusterTool.FuzzyCluster.command="cini -r2:4|fcm|picK";
 * @param ToolSvc.FuzzyClusterTool.clusterSetNo
 *        clusterSetNo is a param of FuzzyClusterTool that specifies
 *        the index of the partition that should be considered; its default
 *        value is zero.
 * @param ApplicationMgr.DLLs += {"FuzzyCluster"};
 *        The FuzzyCluster library should be loaded if one uses 
 *        FuzzyClustering.
 *
 * @param CalClustersAlg.corrToolNames
 *        name of tools performing corrections
 *
 *
 * <hr>
 * @section notes release notes
 * @include release.notes
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 *  
 *
 * @todo move out hardwired constants describing the energy correction fit
 *       from CalClustersAlg::Leak function into a data file;
 *       update this constants to reflect the actual detector geometry
*/

