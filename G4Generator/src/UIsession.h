#include "G4UIsession.hh"
#include "G4coutDestination.hh"
#include "globals.hh"

/** 
 * @class UIsession
 *
 * @brief A dummy Geant4 user interface 
 *
 * Geant4 has some strange behavior with respect to information streaming on the
 * standard output during initialization phase. Most information can be silenced
 * by choosing proper values for verbosity levels; netherless Geant4 print on
 * the std output a lot of information that can be annoing (and not very
 * useful). To avoid this information to interfere with the standard output of
 * the Gaudi initialization phase this small class shut up the Geant4
 * output. Note that this has been set in such a way that only the
 * initialization phase is silenced; it is still possible, with the jobOptions
 * mechanism, to set the interface of Geant4 to a particular verbose level to
 * obtain information during the run
 *  
 * @author R.Giannitrapani
 *
 * $Header$ 
 */
class UIsession : public G4UIsession
{  
 public:
  /// This dummy method silence the abnormal Geant4 output
  virtual G4int ReceiveG4cout(G4String coutString){return 0;};
};
