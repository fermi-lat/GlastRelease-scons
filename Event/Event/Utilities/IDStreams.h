
#ifndef IDSTREAMS_H
#define IDSTREAMS_H



#include "GaudiKernel/StreamBuffer.h"
#include "idents/VolumeIdentifier.h"

/** 
 * The following stream operators are used for serialization of the class
 * idents::VolumeIdentifier into and from a presistent data store
 *
 * @author Marco Frailis
 *    
 * $Header$
 */


/** 
 * Output operator for the class idents::VolumeIdentifier. Used for serialization
 * in the persistent data store
 */
inline StreamBuffer& operator<< ( StreamBuffer& s, const idents::VolumeIdentifier& id)    
{

  idents::VolumeIdentifier::int64 value = id.getValue();
  unsigned int sup, inf;
  
  // The int64 value is splitted into two unsigned int of 32 bits
  sup = value >> 32;

  // This is a 32 bits mask with all bits set to 1
  unsigned int mask = 0xffffffff;
  inf = value & mask;

  unsigned int size = id.size();

  return s  << size << inf << sup;
}

/** 
 * Input operator for the class idents::VolumeIdentifier. Used for serialization
 * in the persistent data store
 */
inline StreamBuffer& operator>> ( StreamBuffer& s, idents::VolumeIdentifier& id )         
{  
  unsigned int size, inf, sup;;
 
  s >> size >> inf >> sup;

  idents::VolumeIdentifier::int64 sup64 = sup;
  
  idents::VolumeIdentifier::int64 value = inf; 
  value |= (sup << 32); 
  id.init(value, size);

  return s;
}



#endif
