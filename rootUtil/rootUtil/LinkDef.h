/** @file LinkDef.h
    @brief for rootcint
*/

#include "BranchGroup.h"
#include "BgDataHandleBase.h"
#include "BgDataHandle.h"
#include "CompositeEventList.h"
#include "CelEventIDs.h"
#include "CelEventComponent.h"
#include "CelEventEntry.h"
#include "CelFileAndTreeSet.h"
#include "CelFileAndTreeOffset.h"
#include "CelEventIDs.h"
#include "CelEventLink.h"
//#include "CelIndex.h"
#include "CelUtil.h"
#include "ComponentsInfo.h"
#include "TestData.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class BgDataHandleBase+ ;
#pragma link C++ class BranchGroup+ ;

#pragma link C++ class BgDataHandle<char *>+ ;
#pragma link C++ class BgDataHandle<Bool_t>+ ;
#pragma link C++ class BgDataHandle<Char_t>+ ;
#pragma link C++ class BgDataHandle<UChar_t>+ ;
#pragma link C++ class BgDataHandle<UShort_t>+ ;
#pragma link C++ class BgDataHandle<Short_t>+ ;
#pragma link C++ class BgDataHandle<UShort_t>+ ;
#pragma link C++ class BgDataHandle<Int_t>+ ;
#pragma link C++ class BgDataHandle<UInt_t>+ ;
#pragma link C++ class BgDataHandle<Long64_t>+ ;
#pragma link C++ class BgDataHandle<ULong64_t>+ ;
#pragma link C++ class BgDataHandle<Float_t>+ ;
#pragma link C++ class BgDataHandle<Double_t>+ ;

#pragma link C++ class CelEventIDs+ ;
#pragma link C++ class CelEventLink+ ;
#pragma link C++ class CelEventEntry+ ;
#pragma link C++ class CelFileAndTreeSet+ ;
#pragma link C++ class CelFileAndTreeOffset+ ;
//#pragma link C++ class CelIndex+ ;
#pragma link C++ class CelEventComponent+ ;
#pragma link C++ class CompositeEventList+ ;

#pragma link C++ class TestAbstractData+ ;
#pragma link C++ class TestDigi+ ;
#pragma link C++ class TestRecon+ ;

#pragma link C++ namespace CelUtil+ ;

#endif
