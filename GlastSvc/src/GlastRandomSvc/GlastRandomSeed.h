#ifndef _GlastRandomSeed_H
#define _GlastRandomSeed_H

#include <iostream>

struct GlastRandomSeed {
  long m_run;
  long m_seqNo;
};

inline std::istream& operator >> (std::istream& in, GlastRandomSeed& seed)
{
  in >> seed.m_run >> seed.m_seqNo;
  return in;
}

inline std::ostream& operator << (std::ostream& out, 
				  const GlastRandomSeed& seed)
{
  out << seed.m_run << " " << seed.m_seqNo;
  return out;
}

#endif
