// Mainpage for doxygen

/** @mainpage package CalUtil
 *
 * @authors R.Dubois
 *
 * @section description Description
 *
 * This package provides utilities for common use in the CAL.     
 *
 * @section CalFailureModeSvc CalFailureModeSvc
 * CalFailureModeSvc creates a list of large-scale failures in the CAL, 
 * and utilities
 * to search the lists to allow digi and recon algorithms to ignore hits based
 * on those lists.
 *
 * It can take lists of towers, (tower,AFEE) or (tower, Controller) pairs
 * to create the lists of dead objects. It provides a method to see if a 
 * given CalXtalId is contained in the lists.
 *
 * @section jobOptions jobOptions
 *
 * @param CalFailureModeSvc.towerList
 * Provide a list of towers that will be made dead.
 * Format: "tower"
 * @param CalFailureModeSvc.towerAfeeList
 * Provide a list of (tower, AFEE) pairs that will be made dead.
 * Format: "tower_afee". Afee runs 0-3: x+,y+,x-,y-
 * @param CalFailureModeSvc.towerControllerList
 * Provide a list of (tower, controller) pairs that will be made dead.
 * Format: "tower_controller". Controller runs 0-3: x+,y+,x-,y-.
 *
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 * <hr>
 * @section requirements requirements
 * @verbinclude requirements
 * <hr>
 *
 */

