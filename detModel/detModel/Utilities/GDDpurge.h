#ifndef GDDPURGE_H
#define GDDPURGE_H
#include <algorithm>

template<class Seq> void purge(Seq & c){
  typename Seq::iterator i;

  for(i = c.begin(); i != c.end(); i++)
    {
      delete *i;
      *i = 0;
    }
}

#endif // GDDPURGE_H
