#ifndef Event_Ebf_H
#define Event_Ebf_H

#include <iostream>

#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/StreamBuffer.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/ObjectVector.h"

#include "Event/TopLevel/Definitions.h"

/**
 * @class Event
 *
 * @brief TDS for storing ebf format of an event
 *
 * The data is stored as one continuos string of bytes
 * No attempt is made to verify that the data stored is correctly
 * formated ebf.
 * $Header$
 */

extern const CLID& CLID_Ebf;

namespace EbfWriterTds{
  class Ebf : public DataObject{
  public:
    Ebf();
    Ebf(char *newData,unsigned int dataLength);
    virtual ~Ebf();
    static const CLID& classID() {return CLID_Ebf;}
    virtual const CLID& clID() const {return classID();}

    ///Retrieve pointer to the ebf data.
    char *get(unsigned int &dataLength) const;

    ///Store the provided ebf pointer in and delete any previous ones
    void set(char *newData, unsigned int dataLength);

    virtual std::ostream& fillStream(std::ostream &s) const;
    friend std::ostream& operator << (std::ostream &s, const Ebf& obj);
  private:
    ///Pointer to the ebf data
    char *data;
    ///Number of bytes that are stored in data pointer
    unsigned int length;
  };
};
#endif
