#include <stdio.h>
#include <string.h>

#include <iostream>
#include <stdexcept>

#include "eventFile/LPA_Handler.h"
#include "eventFile/EBF_Data.h"

#include "lsfData/LSFReader.h"
#include "lsfData/LsfCcsds.h"
#include "lsfData/LsfDatagramInfo.h"
#include "lsfData/LsfGemScalers.h"
#include "lsfData/LsfGemTime.h"
#include "lsfData/LsfRunInfo.h"
#include "lsfData/LsfTime.h"
#include "lsfData/LsfTimeTone.h"
#include "lsfData/LsfMetaEvent.h"

int main( int argc, char* argv[] )
{
  // create the LPA_File object from which input will be read
  lsfData::LSFReader* pLSF = NULL;
  std::string lsefile( "$(EVENTFILEROOT)/src/test/events.lpa" );
  if ( argc >= 2 ) {
    lsefile = argv[1];
  }
  try {
    pLSF = new lsfData::LSFReader( lsefile );
  } catch( std::runtime_error e ) {
    std::cout << e.what() << std::endl;
  }

  // print the header summary information
  printf( "retrieved %llu events for run %09u\n", pLSF->evtcnt(), pLSF->runid() );
  printf( "first Spacecraft clock value = %u\n", pLSF->begSec() );
  printf( "last  Spacecraft clock value = %u\n", pLSF->endSec() );
  printf( "first GEM sequence counter value = %llu\n", pLSF->begGEM() );
  printf( "last  GEM sequence counter value = %llu\n", pLSF->endGEM() );
  for ( int i=0; i<LSEHEADER_MAX_APIDS; i++ ) {
    std::pair<unsigned, unsigned> err = pLSF->seqErr( i );
    if ( err.first > 0 ) {
      printf( "apid %04u had %10u sequence errors\n", err.first, err.second );
    }
  }
  for ( int i=0; i<LSEHEADER_MAX_APIDS; i++ ) {
    std::pair<unsigned, unsigned> err = pLSF->dfiErr( i );
    if ( err.first > 0 ) {
      printf( "apid %04u had %10u DFI errors\n", err.first, err.second );
    }
  }
  printf( "\n" );

  // declare objects to receive the event information
  lsfData::LsfCcsds lccsds;
  lsfData::MetaEvent lmeta;
  eventFile::EBF_Data ebf;

  // retrieve each event in turn
  bool bmore = true;
  do {
    try {
      bmore = pLSF->read( lccsds, lmeta, ebf );
    } catch( std::runtime_error e ) {
      std::cout << e.what() << std::endl;
      break;
    }
    if ( !bmore ) break;

    // print out the context
    printf( "==========================================\n" );
    printf( "\nEvent %lld CCSDS:", lmeta.scalers().sequence() );
    printf( "\n------------------\n" );
    lccsds.print();

    printf( "\nEvent %lld info:", lmeta.scalers().sequence() );
    printf( "\n---------------\n" );
    lmeta.run().print();
    lmeta.datagram().print();
    lmeta.scalers().print();
    lmeta.time().print();
    lmeta.configuration()->print();
    lmeta.keys()->print();

    // iterate over the contributions and print them out
    printf( "\nEvent %lld data:", lmeta.scalers().sequence() );
    printf( "\n---------------\n" );
    printf( "%d bytes of EBF\n", ebf.size() );
    printf( "\n" );

  } while ( true );
  delete pLSF;

  // all done
  return 0;
}
