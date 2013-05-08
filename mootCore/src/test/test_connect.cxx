// $Header$

// Try out readonly and read/write type of connection.  

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include "mootCore/MoodConnection.h"
#include "facilities/commonUtilities.h"

int main(int, char**)    {

  facilities::commonUtilities::setupEnvironment();
  MOOT::MoodConnection* moodCon = 
    new MOOT::MoodConnection(false);  // read only


  if (!(moodCon) ) {
    std::cerr << std::endl << "unable to connect to MOOD for read "
              << std::endl << std::endl;
    std::cerr.flush();
    return 1;
  }
  else if (!(moodCon->getRdb())) {
    
    std::cerr << std::endl << "unable to connect to MOOD for read  "
              << std::endl << std::endl;
    std::cerr.flush();
    return 1;
  }

  else {
    std::cout << std::endl << "Successfully connected "
              << " for read " << std::endl << std::endl 
              << std::endl << std::endl;
  }
  
  delete moodCon;
  moodCon = 0;

    
  moodCon = new MOOT::MoodConnection("true");

  if (!(moodCon) ) {
    std::cerr << std::endl << "unable to connect to MOOD for write" 
              << std::endl << std::endl;
    std::cerr.flush();
    return 1;
  }
  else if (!(moodCon->getRdb())) {
    std::cerr << std::endl << "unable to connect to MOOD for write" 
              << std::endl << std::endl;
    std::cerr.flush();
    return 1;
  }

  else {
    std::cout << std::endl << "Successfully connected " 
              << " for write" << std::endl << std::endl 
              << std::endl << std::endl;
  }
  return 0;
}
