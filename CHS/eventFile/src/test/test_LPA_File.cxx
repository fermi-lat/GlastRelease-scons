#include <stdio.h>
#include <string.h>

#include <iostream>
#include <stdexcept>

#include "eventFile/LPA_File.h"

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
  eventFile::LPA_File* pLPA = NULL;
  std::string lpafile( "$(EVENTFILEROOT)/src/test/events.lpa" );
  if ( argc >= 2 ) {
    lpafile = argv[1];
  }
  try {
    pLPA = new eventFile::LPA_File( lpafile, eventFile::LPA_File::Read );
  } catch( std::runtime_error e ) {
    std::cout << e.what() << std::endl;
  }

  std::cout << "Data file for run " << pLPA->runid() << std::endl;

  // declare objects to receive the event information
  eventFile::LSE_Context ctx;
  eventFile::EBF_Data    ebf;
  eventFile::LPA_Info    info;

  // retrieve each event in turn
  bool bmore = true;
  do {
    try {
      bmore = pLPA->read( ctx, info, ebf );
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
    info.dump();

    // iterate over the contributions and print them out
    printf( "\nEvent %lld data:", ctx.scalers.sequence );
    printf( "\n---------------\n" );
    MyEBFeventIterator eei;
    eei.iterate( const_cast< EBFevent* >( ebf.start() ), 
		 const_cast< EBFevent* >( ebf.end() ) );
    printf( "\n" );

  } while ( true );
  delete pLPA;

  // all done
  return 0;
}
