/** @mainpage package digiRootData

  @section intro Introduction
  This package defines the ROOT classes used to store GLAST digitization data.
  The main branch is defined as the DigiEvent class.  This class contains:
  - Run Id
  - Event Id
  - Collection of CalDigi objects.
    -# Mode:  BESTRANGE or ALLRANGE
    -# CalXtalId - packed UInt_t storing tower, layer, column
    -# Collection of one (BESTRANGE) or four (ALLRANGE) CalXtalReadout objects\n
       Each readout contains:  range and ADC values for each end (POS and NEG).
       

  <hr>
  @section notes release notes
  @include release.notes
  @section requirements requirements
  @include requirements

  @authors Dan Flath and Heather Kelly
*/

