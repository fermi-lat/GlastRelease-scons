// $Header $

/** \file Test program for IdKey class */

// Could use some work.  For one thing, there is no interface to
// the outside world except the debugger.

#include  "xmlUtil/id/IdKey.h"
#include  <vector>
#include  <map>

int main(int argc, char* argv[]) {

  std::vector<unsigned int> a, b, c, d, e;

  for (int i = 0; i < 15; i++) {
    a.push_back(i);
    d.push_back(i);
    b.push_back(15-i);
    c.push_back(0);
  }

  d.push_back(0);
  e.push_back(20);

  xmlUtil::IdKey* aKey = new xmlUtil::IdKey(a);
  xmlUtil::IdKey* bKey = new xmlUtil::IdKey(b);
  xmlUtil::IdKey* cKey = new xmlUtil::IdKey(c);
  xmlUtil::IdKey* dKey = new xmlUtil::IdKey(d);
  xmlUtil::IdKey* eKey = new xmlUtil::IdKey(e);

  // Now make a couple keys from deque
  std::deque<unsigned int> deqy, deqz;

  deqy.push_back(2);
  deqy.push_back(3);
  deqy.push_front(17);

  deqz.push_back(0);
  deqz.push_back(0);
  deqz.push_back(0);
  deqz.push_back(0);
  deqz.push_back(0);
  deqz.push_back(0);

  xmlUtil::IdKey* deqyKey = new xmlUtil::IdKey(deqy);
  xmlUtil::IdKey* deqzKey = new xmlUtil::IdKey(deqz);


  std::map<xmlUtil::IdKey, const char*, xmlUtil::IdKey::ltkey> idMap;

  const char *aName="aKey";
  const char *bName="bKey";
  const char *cName="cKey";
  const char *dName="dKey";
  const char *eName="eKey";
  const char *deqyName="deqy";
  const char *deqzName="deqz";

  idMap[*aKey] = aName;
  idMap[*bKey] = bName;
  idMap[*cKey] = cName;
  idMap[*dKey] = dName;
  idMap[*eKey] = eName;
  idMap[*deqyKey] = deqyName;
  idMap[*deqzKey] = deqzName;

  return 1;
}
