xx/** @mainpage package LdfConverter

  @section intro Introduction
  Handles ldf or ebf input, either raw or "FITSified", storing digis in the usual
  place in the TDS.  It depends on the ldfReader package to parse the
  input file.

  The LDF Online's library EventSummary::marker routine is used to determine if
  the retrieved event is a data event or not.  A marker == 0 denotes a data
  event (from Jim Panette 6/22/2004) non-zero makrers are non-data events..where
  marker == 5 is a "sweep event"...whatever that means..
  <hr>
  @section jobOptions jobOptions
  @param LdfEventSelector.StorageType
  REQUIRED. String which must be one of supported types, currently 
  "CCSDSFILE" for LSF
  "LDFFILE" for raw ldf or  "LDFFITS" for FITSified ldf or
  "EBFFILE" for raw ebf or "EBFFITS" for (slightly) FITSified ebf or
  "NONE" for no real input at all.
  @param InputList
  REQUIRED String array of input data file specs.  Currently (May, 2003)
  only the first is used.
  @param  EvtMax
  REQUIRED  Integer. Max. number of events to be processed
  @param  LdfEventSelector.Instrument
  String. Unused for now (with external EBF library v3 and later)
  @param LdfEventSelector.StartEventIndex
  Integer.  Allows user to request a specific event identified by its index
  (placement) within the input.
  @param LdfEventSelector.StartEventNumber
  Integer.  Allows user to request a specific event number (as reflected in 
  the data) as the stating event.
  @param LdfEventSelector.EbfDebugLevel
  Integer.  Defaults to 0 (no debug output).  Anything else enables 
  debug output
  @param LdfEventSelector.SweepEventSearch
  Integer.  Defaults to 1 (search for first sweep event).  Setting to zero
  skips this search.
  @param LdfEventSelector.GemZeroCheck
  Integer.  Deafults to 0 (skip check for GEM condition summary == 0).
  Setting to one performs this check and prints warnings messages for the
  first five events where Condition Summary is zero.
  @param LdfEventSelector.SocketConnection
  Integer.  Defaults to 0.  If set to non-zero, a UDP client is created to 
  accept incoming LDF data.
  @param LdfEventSelector.SocketServer
  Integer.  Defaults to 60.  If SocketConnection is set non-zero, the 
  SocketServer JO parameter detemines the server number used for listening.
  The beamtest uses a default of 60.  This value can be in [1,99].
  @param LdfEventSelector.AcdRemapFile
  String.  Defaults to "" (null string).  If set to a non-null string, the
  conversion of the ACD PMTs will look to the specified text file to 
  determine the remapping to be used for the tileIds.  This has been added
  strictly from beamtest2006 and is not meant to be used for typical LDF
  processing.
  @param LdfEventSelector.IgnoreSegFault
  Integer.  Defaults to 0 (false).  If set to non-zero value, this will set 
  the signal handles for a segmentation fault to ignore to allow processing
  to continue...not to be used lightly.
  @param LdfEventSelector.OldStyleRunId
  Integer.  Defaults to 0 (false).  If set to non-zero value, this will cause
  the runId to be obtained the same way done for I&T (via LSEReader::runid())
  rather than use the method for flight using ctx.run.startedAt
  <hr>
  @section notes release notes
  release.notes
  <hr>
  @section requirements requirements
  @verbinclude requirements

  @author Joanne Bogart and Heather Kelly

*/

