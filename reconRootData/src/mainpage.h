
/** @mainpage package reconRootData
*
* @section intro Introduction
*
* ROOT classes used to store GLAST reconstruction data. 
* The reconstruction data consist of TKR, CAL, and ACD data.
* ReconEvent is the main branch on the tree.
*
* @image html reconTree.gif
* 
* @section ReconEvent ReconEvent
* The top level object in the Recon Tree is ReconEvent, which contains:
* - Run Number
* - Event Id
* - A pointer to the ACD reconstruction object.
* - A pointer to the CAL reconstruction object.
* - A pointer to the TKR reconstruction object.
* - Event Flags - current a single bit packed UInt_t
*
* @section AcdRecon AcdRecon
* The AcdRecon object provides access to some simple ACD reconstructed quantities.
* -# Total Reconstructed Energy detected
* -# Number of ACD tiles found to be above veto threshold
* -# Minimum DOCA (Distance of Closest Approach)
* -# Minimum Active Distance
* -# Collection of minimum DOCA values for top and side rows.
* -# Collection of minimum Active Distance values for top and side rows.
* -# Gamma DOCA which is computed using all ACD tiles and the reconstructed gamma
*
* @section CalRecon CalRecon
* The CalRecon object provides access to two collections:
* -# CalXtalRecData which are reconstructed CAL crystals
* -# CalCluster which are clusters of crystals
*
* @section TkrRecon TkrRecon
* The TkrRecon object provides access to four collections of reconstructed 
* TKR objects:
* -# TkrCluster which are clusters of silicon hits
* -# TkrCandTrack which are pattern recognized tracks
* -# TkrTrack which are reconstructed tracks
* -# TkrVertex which are vertices
*
* <hr>
* @section notes release notes
* release.notes
* <hr>
* @section requirements requirements
* @verbinclude requirements
*
*/

