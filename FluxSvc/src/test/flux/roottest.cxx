// Flux test program that generates a ROOT macro to plot the flux
//

#include "rootplot.h"


main(int argc, char** argv)
{
  std::vector<char*> arguments;
  arguments.push_back("chime");
  arguments.push_back("-min");
  arguments.push_back("10");
  rootplot abc(arguments);
   //return 0;
}
