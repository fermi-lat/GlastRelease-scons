// $Header$


#include "enums/GemConditionSummary.h"
#include "enums/EventFlags.h"
#include "enums/TriggerBits.h"

#include <iostream>

int main() 
{
    unsigned int m_gemCondSummary  = enums::GEM_mask;

    unsigned m_triggerWord = 0;
    int ibit = enums::number_of_trigger_bits;
    while(ibit--) { m_triggerWord |= 1<<ibit; }

    unsigned int m_eventFlag = 31;

    std::cout << "Test EventFlag Enum " << std::endl;
    std::cout << "Test EvtSeq Bit " << (m_eventFlag & enums::EVTSEQ) << " "
        << "Test TkrRecon bit " << (m_eventFlag & enums::TKRRECON) << std::endl;
    std::cout << "Test PacketError Bit " << (m_eventFlag & enums::PACKETERROR) 
        << ", SummaryError Bit " << (m_eventFlag & enums::SUMMARYERROR)
        << ", ParityError Bit " << ( m_eventFlag & enums::TRGPARITYERROR)
        << std::endl << std::endl;

    std::cout << "Test Trigger Word Enum " << std::endl;
    std::cout << "Test ACDL bit " << (m_triggerWord & enums::b_ACDL) << " "
        << "Test 3-in-a-row-bit " << (m_triggerWord & enums::b_Track) << std::endl;
    std::cout << "Test Low CAL bit " << (m_triggerWord & enums::b_LO_CAL) << " "
        << "Test Hi CAL bit " << (m_triggerWord & enums::b_HI_CAL) << std::endl;
    std::cout << "Test ACD CNO bit " << (m_triggerWord & enums::b_ACDH) << " "
        << "Test Ritz Throttle bit " << (m_triggerWord & enums::b_THROTTLE) 
        << std::endl << std::endl;

    std::cout << "Gem Condition Summary Bits" << std::endl;
    std::cout << "Test ROI bit " << (m_gemCondSummary & enums::ROI) << " "
        << "Test TKR bit " << (m_gemCondSummary & enums::TKR) << std::endl;
    std::cout << "Test CALLE " << (m_gemCondSummary & enums::CALLE) << " "
        << "Test CALHE " << (m_gemCondSummary & enums::CALHE) << std::endl;
    std::cout << "Test CNO bit " << (m_gemCondSummary & enums::CNO) << " "
        << "Test Periodic bit " << (m_gemCondSummary & enums::PERIODIC) << std::endl;
    std::cout << "Test Solicited bit " << (m_gemCondSummary & enums::SOLICITED) << " "
        << "Test External bit " << (m_gemCondSummary & enums::EXTERNAL) << std::endl;

  return 0;
}
