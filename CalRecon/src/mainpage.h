// $Header$
// Mainpage for doxygen

/** @mainpage package CalRecon
 *
 * @authors A.Chekhtman, R.Terrier, J.A.Hernando, G. Musat, B. Giebels
 *
 * @section intro Introduction
 *
 *  CalRecon package reconstructs the energy and direction  of  incident particle from
 *   the calorimeter information. 
 *   
 *   Package contains 3 algorithms: CalXtalRecAlg, CalClustersAlg
 *   and CalDisplay.
 *
 *   The  control flow diagram of CalRecon package, including also the CalDigi package,
 *   is given on the following WEB page:
 * <A HREF="http://www-glast.slac.stanford.edu/software/CAL/meetings/calrecondiagram.pdf"> CalRecon diagram </A> 
 *
 *   CalXtalRecAlg takes the digitized calorimeter information from CalDigiCol
 *   as input, calculates the energy and position in each hitted crystal
 *   and stores this data into CalXtalRecCol.  CalXtalResponse package is 
 *   used for the estimation of energy & position from digi info.  See
 *   documentation in CalXtalResponse for details.
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
 *
 *   CalDisplay algorithm provides the display of reconstructed data.
 *
 * @section Tools Tools
 *
 * The original clustering and energy corrections (profile and last-layer
 * correlation methods) have been recast as Gaudi tools.
 *
 * The clustering is now based on the CalIClusteringTool and ClusteringTool base classes, 
 * while the leakage corrections derive from IEnergyCorr and EnergyCorr.
 * The clustering tools currently available are CalSingleClusteringTool and
 * CalSimpleClusteringTool (FuzzyClusteringTool needs additional tuning). The last
 * layer leakage tool is LastLayerCorrTool, and the profile tool is
 * ProfileTool. The tool using energy correlation with the number of hits in the
 * tracker is CalTkrLikelihoodTool
 *
 * The calo FuzzyClusteringTool retrieves the generic FuzzyCluster Gaudi
 * tool to perform fuzzy clustering on all cal hits. However, if there are
 * less than one hit point per cluster, it applies one single cluster
 * calculations, like the CalSingleClusteringTool. 
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
 *        name of tool performing clustering. Default is CalSingleClusteringTool
 *        If set to FuzzyClusteringTool, the following param is mandatory:
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
 *        FuzzyClusteringTool.
 *
 * @param CalClustersAlg.corrToolNames
 *        name of tools performing corrections
 *
 *
 * @param CalXtalRecAlg.xtalEneTool
 *        name of CalXtalResponse/IXtalEneTool based tool performing xtal digi->energy conversion (default is "XtalEneTool")
 * @param CalXtalRecAlg.xtalPosTool
 *        name of CalXtalResponse/IXtalPosTool based tool performing xtal digi->position conversion (default is "XtalPosTool")
 *
 * <hr>
 * @section notes release notes
 * @include release.notes
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 *  
 * @todo implement real clustering in CalClustersAlg to determine the 
 *       energies of electron and positron produced
 *       by low energy photon
 *
 * @todo move out hardwired constants describing the energy correction fit
 *       from CalClustersAlg::Leak function into a data file;
 *       update this constants to reflect the actual detector geometry
*/

