#ifndef PURGE_H
#define PURGE_H
#include <algorithm>

namespace detModel{

  template<class Seq> void purge(Seq & c){
    typename Seq::iterator i;
    
    for(i = c.begin(); i != c.end(); i++)
      {
	delete *i;
	(*i) = 0;
      }
  }

}
#endif // PURGE_H
