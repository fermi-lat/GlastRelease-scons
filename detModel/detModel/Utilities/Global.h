#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <cstdlib>

namespace detModel{
  /// Some global methods for detModel

  /// This methods can be used to abort detModel with a message (s)
  inline void detAbort(std::string s){
    std::cout << s << std::endl;
    std::cout << "#################" << std::endl;
    std::cout << "Aborting detModel" << std::endl;
    exit(0);
  }

}
#endif // GLOBAL_H

