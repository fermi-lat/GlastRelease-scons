/** @mainpage package eventFile
 *
 * @authors Bryson Lee <blee@slac.stanford.edu>
 * 
 * @section intro Introduction
 * The eventFile package provides a layer of abstraction between the eventRet
 * package and ldfReader / LatIntegration / Gleam.  It has minimal dependencies
 * on external libraries, and can be compiled as a single-threaded library.
 * It supports reading a serialized stream of context+metainfo+event blocks
 * from files created using eventRet
 * 
 * @section depends Dependencies
 * - The \c ldf external package which provides the LDF library.
 *
 * <hr>
 * @section notes release.notes
 * release.notes
 *
 * <hr>
 * @todo implement handling for charge-injection (LSEC) data.
 *
 */
