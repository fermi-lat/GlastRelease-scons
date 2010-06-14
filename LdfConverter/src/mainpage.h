xx/** @mainpage package LdfConverter

  @section intro Introduction
  Handles ebf input, either raw or "FITSified", storing digis in the usual
  place in the TDS.  It depends on the ldfReader package to parse the
  input file.
  <hr>
  @section jobOptions jobOptions
  @param StorageType
  REQUIRED. String which must be one of supported types, currently 
  "EBFFILE" for raw ebf or "EBFFITS" for (slightly) FITSified ebf or
  "NONE" for no real input at all.
  @param InputList
  REQUIRED String array of input data file specs.  Currently (May, 2003)
  only the first is used.
  @param  EvtMax
  REQUIRED  Integer. Max. number of events to be processed
  @param  Instrument
  String. Unused for now (with external EBF library v3 and later)
  @param EbfDebugLevel
  Integer.  Defaults to 0 (no debug output).  Anything else enables 
  debug output
  <hr>
  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements

  @author Joanne Bogart and Heather Kelly

*/

