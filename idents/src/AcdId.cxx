// File and Version information
// $Header$
//
// ClassName:   AcdId
//  
// Description: ID class for ACD components
//
//              
// Author:  H. Kelly     20 Oct 2004

// Include files
#include "idents/AcdId.h"
#include "facilities/Util.h"
#include <stdexcept>
#include <iostream>
#include <ios>

using namespace idents; 

/** Constructor from VolumeIdentifier.  Assuming a particular form for
    fields within VolumeIdentifier (correct as of Oct 20, 2004):
*/
AcdId::AcdId(const VolumeIdentifier& vId) : m_id(0) {
  constructorGuts(vId);
}

void AcdId::constructorGuts(const VolumeIdentifier& volId)   {

    na(0);
    if (!checkVolId(volId)) throw std::invalid_argument("VolumeIdentifier");

    if (volId[2] == tileVolId) {
        face(volId[1]);
        /* as of revised ACD geometry (nov ?? 2002) row and column
        are always in the right order, but there is an intervening
        field which has to be thrown away. For now it should always
        have value 40 (tile).  When ribbons are sensitive, will also see
        41 if volume refers to a ribbon.  */
        row(volId[3]);
        column(volId[4]);
    } else if (volId[2] == ribbonVolId) {
        ribbonNum(volId[4]);
        ribbonOrientation(6-volId[3]);
    }

}
AcdId::AcdId(const std::string& base10) {
  // Following may throw an exception
  int b10 = facilities::Util::stringToInt(base10);

  m_id = 0;
  na(b10/1000);
  b10 = b10 % 1000;
  face(b10/100);
  b10 = b10 % 100;
  if (na() ) column(b10);
  else {
    row(b10/10);
    b10 = b10 % 10;
    column(b10);
  }
}

  





