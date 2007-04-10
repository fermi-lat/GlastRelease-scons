/** @file HFCFilterOutput.h

* @class HFCFilterOutput
*
* @brief This class handles the End of Event and End of Run "output" processing for the Gamma Filter 
*
* last modified 12/04/2006
*
* @authors T. Usher
*
* $Header$
*/

#ifndef __HFCFilterOutput_H
#define __HFCFilterOutput_H

#include "OutputRtn.h"

class HFCFilterOutput : virtual public OutputRtn
{
public:
    HFCFilterOutput(int offset, bool passThrough=false);
    virtual ~HFCFilterOutput() {}

    // This defines the method called for end of event processing
    void eovProcessing(void* callBackParm, EDS_fwIxb* ixb);

    // This for end of run processing
    void eorProcessing(MsgStream& log);
private:
    // Local functions

    int m_offset;         // Offset into ixb event desriptor block for this information
    bool m_passThrough;   // Running filter in pass through mode

    int m_vetoBits[17];   //array to count # of times each veto bit was set
    int m_statusBits[15]; //array to count # of times each veto bit was set
};

#endif // __ObfInterface_H
