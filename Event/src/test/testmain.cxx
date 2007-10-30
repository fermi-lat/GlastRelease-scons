#include "Event/TopLevel/EventModel.h"
#include "geometry/Vector.h"
#include "Event/Utilities/TimeStamp.h"
#include "Event/Digi/GltDigi.h"

#include <iostream>

/// test GltDigi/CalTrigVector
/// @author Zachary Fewtrell
/// @return non zero on failure
int test_GltDigi_CalTrigVector() {
  // basically populate a few trigger bits & check that the correct bits are set and not set

  Event::GltDigi gltDigi;

  // expected FLE trigger vector
  unsigned short fle_answer_vector = 0;
  // expected FHE trigger vector
  unsigned short fhe_answer_vector = 0;

  // set FLE T0 L0 C0 F0
  gltDigi.setCALLOtrigger(idents::CalXtalId(0,0,0,0));
  // set expected tower bit
  fle_answer_vector |= 1<<0;

  // set FHE T7 L5 C7 F1
  gltDigi.setCALHItrigger(idents::CalXtalId(7,5,7,1));
  // set expected tower bit
  fhe_answer_vector |= 1<<7;

  // set FLE T15 L7 C11 F1
  gltDigi.setCALLOtrigger(idents::CalXtalId(15,7,11,1));
  // set expected tower bit
  fle_answer_vector |= 1<<15;

  if (fle_answer_vector != gltDigi.getCALLETriggerVector()) {
    std::cerr << "Cal FLE trigger vector, unexpected answer" << std::endl;
    return -1;
  }

  if (fhe_answer_vector != gltDigi.getCALHETriggerVector()) {
    std::cerr << "Cal FLE trigger vector, unexpected answer" << std::endl;
    return -1;
  }

  return 0;
}

/// return non-zero on success
int main(){
  return test_GltDigi_CalTrigVector();
}
