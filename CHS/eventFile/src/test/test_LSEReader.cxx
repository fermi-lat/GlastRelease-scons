#include <stdio.h>
#include <string.h>

#include <iostream>
#include <stdexcept>

#include "eventFile/LSEReader.h"

#include "eventFile/LSE_Context.h"
#include "eventFile/LSE_Info.h"
#include "eventFile/EBF_Data.h"

#include "LATdatagramIterator.h"
#include "EBFeventIterator.h"

class EBFevent;
class MyLATcomponentIterator;
class MyEBFeventIterator : public EBFeventIterator
{
public:
  MyEBFeventIterator();
  virtual ~MyEBFeventIterator();

  virtual int handleError(EBFevent* evt, unsigned code, unsigned p1=0, unsigned p2=0) const;
  virtual int process(EBFevent*);
private:
  MyLATcomponentIterator* _lci;
};

int main( int argc, char* argv[] )
{
  // create the LPA_File object from which input will be read
  eventFile::LSEReader* pLSE = NULL;
  std::string lsefile( "$(EVENTFILEROOT)/src/test/events.lpa" );
  if ( argc >= 2 ) {
    lsefile = argv[1];
  }
  try {
    pLSE = new eventFile::LSEReader( lsefile );
  } catch( std::runtime_error e ) {
    std::cout << e.what() << std::endl;
  }

  // print the header summary information
  printf( "retrieved %llu events for run %09u\n", pLSE->evtcnt(), pLSE->runid() );
  printf( "first Spacecraft clock value = %u\n", pLSE->begSec() );
  printf( "last  Spacecraft clock value = %u\n", pLSE->endSec() );
  printf( "first GEM sequence counter value = %llu\n", pLSE->begGEM() );
  printf( "last  GEM sequence counter value = %llu\n", pLSE->endGEM() );
  for ( int i=0; i<LSEHEADER_MAX_APIDS; i++ ) {
    std::pair<unsigned, unsigned> err = pLSE->seqErr( i );
    if ( err.first > 0 ) {
      printf( "apid %04u had %10u sequence errors\n", err.first, err.second );
    }
  }
  for ( int i=0; i<LSEHEADER_MAX_APIDS; i++ ) {
    std::pair<unsigned, unsigned> err = pLSE->dfiErr( i );
    if ( err.first > 0 ) {
      printf( "apid %04u had %10u DFI errors\n", err.first, err.second );
    }
  }
  printf( "\n" );

  // declare objects to receive the event information
  eventFile::LSE_Context ctx;
  eventFile::EBF_Data    ebf;
  eventFile::LSE_Info::InfoType infotype;
  eventFile::LPA_Info    pinfo;
  eventFile::LCI_ACD_Info ainfo;
  eventFile::LCI_CAL_Info cinfo;
  eventFile::LCI_TKR_Info tinfo;

  // retrieve each event in turn
  bool bmore = true;
  do {
    try {
      bmore = pLSE->read( ctx, ebf, infotype, pinfo, ainfo, cinfo, tinfo );
    } catch( std::runtime_error e ) {
      std::cout << e.what() << std::endl;
      break;
    }
    if ( !bmore ) break;

    // print out the context
    printf( "==========================================\n" );
    printf( "\nEvent %lld context:", ctx.scalers.sequence );
    printf( "\n------------------\n" );
    ctx.dump();

    // print out the event-type-specific information
    printf( "\nEvent %lld info:", ctx.scalers.sequence );
    printf( "\n---------------\n" );
    switch ( infotype ) {
    case eventFile::LSE_Info::LPA:
      pinfo.dump();
      break;
    case eventFile::LSE_Info::LCI_ACD:
      ainfo.dump();
      break;
    case eventFile::LSE_Info::LCI_CAL:
      cinfo.dump();
      break;
    case eventFile::LSE_Info::LCI_TKR:
      tinfo.dump();
      break;
    }

    // iterate over the contributions and print them out
    printf( "\nEvent %lld data:", ctx.scalers.sequence );
    printf( "\n---------------\n" );
    printf( "%d bytes of EBF\n", ebf.size() );
//     MyEBFeventIterator eei;
//     eei.iterate( const_cast< EBFevent* >( ebf.start() ), 
// 		 const_cast< EBFevent* >( ebf.end() ) );
    printf( "\n" );

  } while ( true );
  delete pLSE;

  // all done
  return 0;
}
