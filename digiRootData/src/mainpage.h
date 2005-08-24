/** @mainpage package digiRootData

  @section intro Introduction
  This package defines the ROOT classes used to store GLAST digitization data.
  The main branch is defined as the DigiEvent class.  
  
  @image html digiTree.gif

  The main branch DigiEvent contains:
  - Run Id
  - Event Id
  - Time Stamp in seconds
  - Live Time
  - Flag denoting if this data is simulated or not
  - L1T (Level One Trigger)
  - Collection of CalDigi objects.
    -# Mode:  BESTRANGE or ALLRANGE
    -# CalXtalId - packed UInt_t storing tower, layer, column
    -# Collection of one (BESTRANGE) or four (ALLRANGE) CalXtalReadout objects\n
       Each readout contains:  range and ADC values for each end (POS and NEG).
  - Collection of TkrDigi objects
    -# Tower and layer ids
    -# List of hit strips
    -# ToT values
  - Collection of AcdDigi objects
    -# AcdId
    -# 2 PMTs each contains: PHA and 3 discriminators:  low(accept), veto(hit), high (CNO)
    -# MC energy (MeV) deposited as a check on PHA values
  - Some items are optionally populated if this is real data
    - PPC Time Base, stored as two words denoting the 2 32bit PPC registers
    - Original clock from LDF at 30 Hz, stored as two words sec and nano-seconds
    - EventSummaryData (real data)
    -# Collection of CalDiagnosticData trigger primitives
    -# Collection of TkrDiagnosticData trigger primitives
    -# Gem (GLT Electronics Module)
    -# ErrorData

  <hr>
  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements

  @authors Dan Flath and Heather Kelly
*/

