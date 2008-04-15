#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "eventFile/LSEReader.h"
#include "eventFile/LSEWriter.h"
#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/LPA_Handler.h"
#include "eventFile/EBF_Data.h"
#include "eventFile/LSE_Keys.h"

struct EvtIdx {
  std::string tag;
  unsigned startedAt;
  unsigned long long sequence;
  unsigned long long elapsed;
  unsigned long long livetime;
  double evtUTC;
  unsigned scid;
  unsigned apid;
  unsigned datagrams;
  unsigned dgmevt;
  std::string oaction;
  std::string caction;
  double dgmUTC;
  off_t fileofst;
  std::string evtfile;

  EvtIdx( const std::string& idxline ) {
    // make the text line into an input stream
    std::istringstream iss( idxline );

    // extract the components
    iss >> tag;
    iss >> startedAt;
    iss >> sequence;
//     iss >> elapsed;
//     iss >> livetime;
//     iss >> evtUTC;
//     iss >> scid;
    iss >> apid;
    iss >> datagrams;
//     iss >> dgmevt;
    iss >> oaction;
    iss >> caction;
//     iss >> dgmUTC;
    iss >> fileofst;
    iss >> evtfile;
  };
};

typedef std::map< std::string, eventFile::LSEReader* > lser_map;
typedef lser_map::iterator lser_iter;

struct cleanup {
  void operator() ( lser_map::value_type x ) { delete x.second; }
};

int main( int argc, char* argv[] )
{
  // get the input and output file names
  if ( argc < 4 ) {
    std::cout << "writeMerge: not enough input arguments" << std::endl;
    exit( EXIT_FAILURE );
  }
  std::string idxfile( argv[1] );
  std::string evtfile( argv[2] );  // this contains format specifiers
  int downlinkID = atoi( argv[3] );

  // add support for configurable output file size
  int maxEvents = -1;
  if ( argc >= 5 ) {
    maxEvents = atoi( argv[4] );
  }
  int currMax = maxEvents;

  // make the output-file-size scaling externally tunable
  double WRITEMERGE_CHUNKSCALE = 0.90;
  char* envbuf = getenv( "WRITEMERGE_CHUNKSCALE" );
  if ( envbuf ) {
    WRITEMERGE_CHUNKSCALE = atof( envbuf );
  }
  double WRITEMERGE_CHUNKFLOOR = 0.50;
  envbuf = getenv( "WRITEMERGE_CHUNKFLOOR" );
  if ( envbuf ) {
    WRITEMERGE_CHUNKFLOOR = atof( envbuf );
  }

  // add support for overriding translated LATC master key
  unsigned long overrideLATC = 0xffffffff;
  if ( argc >= 6 ) {
    overrideLATC = strtoul( argv[5], NULL, 0 );
    std::cout << "writeMerge: overriding LATC key to " << overrideLATC << std::endl;
  } else {
    std::cout << "writeMerge: no LATC key override" << std::endl;
  }

  // declare the file-output object pointer
  int eventsOut = 0;
  eventFile::LSEWriter* pLSEW = NULL;

  // create a container for the chunk-evt input files
  lser_map mapLSER;

  // declare object to receive the event information
  eventFile::LSE_Context ctx;
  eventFile::EBF_Data    ebf;
  eventFile::LSE_Info::InfoType infotype;
  eventFile::LPA_Info    pinfo;
  eventFile::LCI_ACD_Info ainfo;
  eventFile::LCI_CAL_Info cinfo;
  eventFile::LCI_TKR_Info tinfo;
  eventFile::LSE_Keys::KeysType ktype;
  eventFile::LPA_Keys     pakeys;
  eventFile::LCI_Keys     cikeys;

  // read the index file and parse the entries, retrieving the 
  // requeseted events as we go
  std::string idxline;
  std::ifstream idx( idxfile.c_str() );
  if ( idx.fail() ) {
    std::cout << "writeMerge: failed to open " << idxfile;
    std::cout << " (" << errno << ":" << strerror(errno) << ")" << std::endl;
    exit( EXIT_FAILURE );
  }
  while ( getline( idx, idxline ) ) {
    
    // skip non-event records
    if ( idxline.find( "EVT:" ) != 0 ) continue;

    // make an event-index object
    EvtIdx edx( idxline );

    // get an LSEReader for the file
    lser_iter it = mapLSER.find( edx.evtfile );
    if ( it == mapLSER.end() ) {
      try {
	it = mapLSER.insert( std::pair< std::string, eventFile::LSEReader* >( edx.evtfile, new eventFile::LSEReader( edx.evtfile ) ) ).first;
      }	catch ( std::runtime_error& e ) {
	std::cout << e.what() << std::endl;
	exit( EXIT_FAILURE );
      }
    }

    // read the event at the specified location
    bool bevtread = false;
    try {
      it->second->seek( edx.fileofst );
      bevtread = it->second->read( ctx, ebf, infotype, pinfo, ainfo, cinfo, tinfo, ktype, pakeys, cikeys );
    } catch( std::runtime_error e ) {
      std::cout << e.what() << std::endl;
      exit( EXIT_FAILURE );
    }
    if ( !bevtread ) {
      std::cout << "no event read from " << edx.fileofst << " of " << edx.evtfile << std::endl;
      exit( EXIT_FAILURE);
    }

    // open an output file if necessary
    if ( !pLSEW ) {
      // create the output filename from the user-supplied template
      char ofn[512];
#ifndef _FILE_OFFSET_BITS
      _snprintf( ofn, 512, evtfile.c_str(), ctx.run.startedAt, ctx.scalers.sequence );
#else
      snprintf( ofn, 512, evtfile.c_str(), ctx.run.startedAt, ctx.scalers.sequence );
#endif

      // open the output file
      try {
	pLSEW = new eventFile::LSEWriter( std::string( ofn ), downlinkID );
      } catch ( std::runtime_error& e ) {
	std::cout << e.what() << std::endl;
	exit( EXIT_FAILURE );
      }
      std::cout << "writeMerge: created output file " << pLSEW->name() << std::endl;
    }

    // write the event to the merged file
    try {
      switch( infotype ) {
      case eventFile::LSE_Info::LPA:
	if ( overrideLATC != 0xffffffff ) {
	  pakeys.LATC_master = overrideLATC;
	}
	pLSEW->write( ctx, ebf, pinfo, pakeys );
	break;
      case eventFile::LSE_Info::LCI_ACD:
	pLSEW->write( ctx, ebf, ainfo, cikeys );
	break;
      case eventFile::LSE_Info::LCI_CAL:
	pLSEW->write( ctx, ebf, cinfo, cikeys );
	break;
      case eventFile::LSE_Info::LCI_TKR:
	pLSEW->write( ctx, ebf, tinfo, cikeys );
	break;
      default:
	std::ostringstream ess;
	ess << "writeMerge: unknown LSE_Info type " << infotype;
	ess << " at offset "            << edx.fileofst;
	ess << " in "                   << edx.evtfile;
	throw std::runtime_error( ess.str() );
	break;
      }
    } catch ( std::runtime_error& e ) {
      std::cout << e.what() << std::endl;
      exit( EXIT_FAILURE );
    }

    // check to see if the output file is full
    if ( currMax > 0 && ++eventsOut >= currMax ) {
      // close the current file and reset the event counter
      std::cout << "writeMerge: wrote " << pLSEW->evtcnt() << " events to " << pLSEW->name() << std::endl;
      delete pLSEW; pLSEW = NULL;
      eventsOut = 0;

      // rescale the max-event count for the next file
      currMax = ( currMax <= WRITEMERGE_CHUNKFLOOR * maxEvents ) ? maxEvents : WRITEMERGE_CHUNKSCALE * currMax;
    }
  }
  std::for_each( mapLSER.begin(), mapLSER.end(), cleanup() );
  idx.close();
  if ( pLSEW ) {
    std::cout << "writeMerge: wrote " << pLSEW->evtcnt() << " events to " << pLSEW->name() << std::endl;
    delete pLSEW;
  }

  // all done
  return 0;
}
