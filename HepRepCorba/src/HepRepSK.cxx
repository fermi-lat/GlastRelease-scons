// This file is generated by omniidl (C++ backend)- omniORB_4_1. Do not edit.

#include "HepRep.hh"
#include <omniORB4/IOP_S.h>
#include <omniORB4/IOP_C.h>
#include <omniORB4/callDescriptor.h>
#include <omniORB4/callHandle.h>
#include <omniORB4/objTracker.h>


OMNI_USING_NAMESPACE(omni)

static const char* _0RL_library_version = omniORB_4_1;



void
HepRepAttDef::operator>>= (cdrStream &_n) const
{
  _n.marshalString(name,0);
  _n.marshalString(desc,0);
  _n.marshalString(category,0);
  _n.marshalString(extra,0);

}

void
HepRepAttDef::operator<<= (cdrStream &_n)
{
  name = _n.unmarshalString(0);
  desc = _n.unmarshalString(0);
  category = _n.unmarshalString(0);
  extra = _n.unmarshalString(0);

}

void
HepRepAttValue::operator>>= (cdrStream &_n) const
{
  _n.marshalString(name,0);
  (const ::CORBA::Any&) value >>= _n;
  showLabel >>= _n;

}

void
HepRepAttValue::operator<<= (cdrStream &_n)
{
  name = _n.unmarshalString(0);
  (::CORBA::Any&)value <<= _n;
  (::CORBA::Long&)showLabel <<= _n;

}

void
HepRepPoint::operator>>= (cdrStream &_n) const
{
  x >>= _n;
  y >>= _n;
  z >>= _n;
  (const HepRepAttValueList&) attValues >>= _n;

}

void
HepRepPoint::operator<<= (cdrStream &_n)
{
  (::CORBA::Double&)x <<= _n;
  (::CORBA::Double&)y <<= _n;
  (::CORBA::Double&)z <<= _n;
  (HepRepAttValueList&)attValues <<= _n;

}

void
HepRepInstance::operator>>= (cdrStream &_n) const
{
  _n.marshalString(typeName,0);
  (const _CORBA_Unbounded_Sequence< HepRepInstance > &) instances >>= _n;
  (const HepRepPointList&) points >>= _n;
  (const HepRepAttValueList&) attValues >>= _n;

}

void
HepRepInstance::operator<<= (cdrStream &_n)
{
  typeName = _n.unmarshalString(0);
  (_CORBA_Unbounded_Sequence< HepRepInstance > &)instances <<= _n;
  (HepRepPointList&)points <<= _n;
  (HepRepAttValueList&)attValues <<= _n;

}

void
HepRepTreeID::operator>>= (cdrStream &_n) const
{
  _n.marshalString(name,0);
  _n.marshalString(version,0);

}

void
HepRepTreeID::operator<<= (cdrStream &_n)
{
  name = _n.unmarshalString(0);
  version = _n.unmarshalString(0);

}

void
HepRepInstanceTree::operator>>= (cdrStream &_n) const
{
  (const HepRepTreeID&) id >>= _n;
  (const HepRepTreeID&) typeTreeID >>= _n;
  (const HepRepTreeIDList&) instanceTreeIDs >>= _n;
  (const HepRepInstanceList&) instances >>= _n;

}

void
HepRepInstanceTree::operator<<= (cdrStream &_n)
{
  (HepRepTreeID&)id <<= _n;
  (HepRepTreeID&)typeTreeID <<= _n;
  (HepRepTreeIDList&)instanceTreeIDs <<= _n;
  (HepRepInstanceList&)instances <<= _n;

}

void
HepRepType::operator>>= (cdrStream &_n) const
{
  _n.marshalString(name,0);
  _n.marshalString(desc,0);
  _n.marshalString(infoURL,0);
  (const _CORBA_Unbounded_Sequence< HepRepType > &) types >>= _n;
  (const HepRepAttDefList&) attDefs >>= _n;
  (const HepRepAttValueList&) attValues >>= _n;

}

void
HepRepType::operator<<= (cdrStream &_n)
{
  name = _n.unmarshalString(0);
  desc = _n.unmarshalString(0);
  infoURL = _n.unmarshalString(0);
  (_CORBA_Unbounded_Sequence< HepRepType > &)types <<= _n;
  (HepRepAttDefList&)attDefs <<= _n;
  (HepRepAttValueList&)attValues <<= _n;

}

void
HepRepTypeTree::operator>>= (cdrStream &_n) const
{
  (const HepRepTreeID&) id >>= _n;
  (const HepRepTypeList&) types >>= _n;

}

void
HepRepTypeTree::operator<<= (cdrStream &_n)
{
  (HepRepTreeID&)id <<= _n;
  (HepRepTypeList&)types <<= _n;

}

void
HepRepAction::operator>>= (cdrStream &_n) const
{
  _n.marshalString(name,0);
  _n.marshalString(expression,0);

}

void
HepRepAction::operator<<= (cdrStream &_n)
{
  name = _n.unmarshalString(0);
  expression = _n.unmarshalString(0);

}

HepRep_ptr HepRep_Helper::_nil() {
  return ::HepRep::_nil();
}

::CORBA::Boolean HepRep_Helper::is_nil(::HepRep_ptr p) {
  return ::CORBA::is_nil(p);

}

void HepRep_Helper::release(::HepRep_ptr p) {
  ::CORBA::release(p);
}

void HepRep_Helper::marshalObjRef(::HepRep_ptr obj, cdrStream& s) {
  ::HepRep::_marshalObjRef(obj, s);
}

HepRep_ptr HepRep_Helper::unmarshalObjRef(cdrStream& s) {
  return ::HepRep::_unmarshalObjRef(s);
}

void HepRep_Helper::duplicate(::HepRep_ptr obj) {
  if( obj && !obj->_NP_is_nil() )  omni::duplicateObjRef(obj);
}

HepRep_ptr
HepRep::_duplicate(::HepRep_ptr obj)
{
  if( obj && !obj->_NP_is_nil() )  omni::duplicateObjRef(obj);
  return obj;
}

HepRep_ptr
HepRep::_narrow(::CORBA::Object_ptr obj)
{
  if( !obj || obj->_NP_is_nil() || obj->_NP_is_pseudo() ) return _nil();
  _ptr_type e = (_ptr_type) obj->_PR_getobj()->_realNarrow(_PD_repoId);
  return e ? e : _nil();
}


HepRep_ptr
HepRep::_unchecked_narrow(::CORBA::Object_ptr obj)
{
  if( !obj || obj->_NP_is_nil() || obj->_NP_is_pseudo() ) return _nil();
  _ptr_type e = (_ptr_type) obj->_PR_getobj()->_uncheckedNarrow(_PD_repoId);
  return e ? e : _nil();
}

HepRep_ptr
HepRep::_nil()
{
#ifdef OMNI_UNLOADABLE_STUBS
  static _objref_HepRep _the_nil_obj;
  return &_the_nil_obj;
#else
  static _objref_HepRep* _the_nil_ptr = 0;
  if( !_the_nil_ptr ) {
    omni::nilRefLock().lock();
    if( !_the_nil_ptr ) {
      _the_nil_ptr = new _objref_HepRep;
      registerNilCorbaObject(_the_nil_ptr);
    }
    omni::nilRefLock().unlock();
  }
  return _the_nil_ptr;
#endif
}

const char* HepRep::_PD_repoId = "IDL:HepRep:1.0";


_objref_HepRep::~_objref_HepRep() {
  
}


_objref_HepRep::_objref_HepRep(omniIOR* ior, omniIdentity* id) :
   omniObjRef(::HepRep::_PD_repoId, ior, id, 1)
   
   
{
  _PR_setobj(this);
}

void*
_objref_HepRep::_ptrToObjRef(const char* id)
{
  if( id == ::HepRep::_PD_repoId )
    return (::HepRep_ptr) this;
  
  if( id == ::CORBA::Object::_PD_repoId )
    return (::CORBA::Object_ptr) this;

  if( omni::strMatch(id, ::HepRep::_PD_repoId) )
    return (::HepRep_ptr) this;
  
  if( omni::strMatch(id, ::CORBA::Object::_PD_repoId) )
    return (::CORBA::Object_ptr) this;

  return 0;
}

// Proxy call descriptor class. Mangled signature:
//  _cHepRepInstanceTree_i_cstring_i_cstring
class _0RL_cd_05fda001c74a4a06_00000000
  : public omniCallDescriptor
{
public:
  inline _0RL_cd_05fda001c74a4a06_00000000(LocalCallFn lcfn,const char* op_,size_t oplen,_CORBA_Boolean upcall=0):
     omniCallDescriptor(lcfn, op_, oplen, 0, _user_exns, 0, upcall)
  {
    
  }
  
  void marshalArguments(cdrStream&);
  void unmarshalArguments(cdrStream&);

  void unmarshalReturnedValues(cdrStream&);
  void marshalReturnedValues(cdrStream&);
  
  
  static const char* const _user_exns[];

  ::CORBA::String_var arg_0_;
  const char* arg_0;
  ::CORBA::String_var arg_1_;
  const char* arg_1;
  HepRepInstanceTree_var result;
};

void _0RL_cd_05fda001c74a4a06_00000000::marshalArguments(cdrStream& _n)
{
  _n.marshalString(arg_0,0);
  _n.marshalString(arg_1,0);

}

void _0RL_cd_05fda001c74a4a06_00000000::unmarshalArguments(cdrStream& _n)
{
  arg_0_ = _n.unmarshalString(0);
  arg_0 = arg_0_.in();
  arg_1_ = _n.unmarshalString(0);
  arg_1 = arg_1_.in();

}

void _0RL_cd_05fda001c74a4a06_00000000::marshalReturnedValues(cdrStream& _n)
{
  (const HepRepInstanceTree&) result >>= _n;

}

void _0RL_cd_05fda001c74a4a06_00000000::unmarshalReturnedValues(cdrStream& _n)
{
  result = new HepRepInstanceTree;
  (HepRepInstanceTree&)result <<= _n;

}

const char* const _0RL_cd_05fda001c74a4a06_00000000::_user_exns[] = {
  0
};

// Local call call-back function.
static void
_0RL_lcfn_05fda001c74a4a06_10000000(omniCallDescriptor* cd, omniServant* svnt)
{
  _0RL_cd_05fda001c74a4a06_00000000* tcd = (_0RL_cd_05fda001c74a4a06_00000000*)cd;
  _impl_HepRep* impl = (_impl_HepRep*) svnt->_ptrToInterface(HepRep::_PD_repoId);
  tcd->result = impl->getInstanceTreeTop(tcd->arg_0, tcd->arg_1);


}

HepRepInstanceTree* _objref_HepRep::getInstanceTreeTop(const char* instanceTreeName, const char* instanceTreeVersion)
{
  _0RL_cd_05fda001c74a4a06_00000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_10000000, "getInstanceTreeTop", 19);
  _call_desc.arg_0 = instanceTreeName;
  _call_desc.arg_1 = instanceTreeVersion;

  _invoke(_call_desc);
  return _call_desc.result._retn();


}
// Proxy call descriptor class. Mangled signature:
//  _cHepRepTypeTree_i_cstring_i_cstring
class _0RL_cd_05fda001c74a4a06_20000000
  : public omniCallDescriptor
{
public:
  inline _0RL_cd_05fda001c74a4a06_20000000(LocalCallFn lcfn,const char* op_,size_t oplen,_CORBA_Boolean upcall=0):
     omniCallDescriptor(lcfn, op_, oplen, 0, _user_exns, 0, upcall)
  {
    
  }
  
  void marshalArguments(cdrStream&);
  void unmarshalArguments(cdrStream&);

  void unmarshalReturnedValues(cdrStream&);
  void marshalReturnedValues(cdrStream&);
  
  
  static const char* const _user_exns[];

  ::CORBA::String_var arg_0_;
  const char* arg_0;
  ::CORBA::String_var arg_1_;
  const char* arg_1;
  HepRepTypeTree_var result;
};

void _0RL_cd_05fda001c74a4a06_20000000::marshalArguments(cdrStream& _n)
{
  _n.marshalString(arg_0,0);
  _n.marshalString(arg_1,0);

}

void _0RL_cd_05fda001c74a4a06_20000000::unmarshalArguments(cdrStream& _n)
{
  arg_0_ = _n.unmarshalString(0);
  arg_0 = arg_0_.in();
  arg_1_ = _n.unmarshalString(0);
  arg_1 = arg_1_.in();

}

void _0RL_cd_05fda001c74a4a06_20000000::marshalReturnedValues(cdrStream& _n)
{
  (const HepRepTypeTree&) result >>= _n;

}

void _0RL_cd_05fda001c74a4a06_20000000::unmarshalReturnedValues(cdrStream& _n)
{
  result = new HepRepTypeTree;
  (HepRepTypeTree&)result <<= _n;

}

const char* const _0RL_cd_05fda001c74a4a06_20000000::_user_exns[] = {
  0
};

// Local call call-back function.
static void
_0RL_lcfn_05fda001c74a4a06_30000000(omniCallDescriptor* cd, omniServant* svnt)
{
  _0RL_cd_05fda001c74a4a06_20000000* tcd = (_0RL_cd_05fda001c74a4a06_20000000*)cd;
  _impl_HepRep* impl = (_impl_HepRep*) svnt->_ptrToInterface(HepRep::_PD_repoId);
  tcd->result = impl->getTypeTree(tcd->arg_0, tcd->arg_1);


}

HepRepTypeTree* _objref_HepRep::getTypeTree(const char* typeTreeName, const char* typeTreeVersion)
{
  _0RL_cd_05fda001c74a4a06_20000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_30000000, "getTypeTree", 12);
  _call_desc.arg_0 = typeTreeName;
  _call_desc.arg_1 = typeTreeVersion;

  _invoke(_call_desc);
  return _call_desc.result._retn();


}
// Proxy call descriptor class. Mangled signature:
//  _cHepRepInstanceTree_i_cstring_i_cstring_i_cStringArray
class _0RL_cd_05fda001c74a4a06_40000000
  : public omniCallDescriptor
{
public:
  inline _0RL_cd_05fda001c74a4a06_40000000(LocalCallFn lcfn,const char* op_,size_t oplen,_CORBA_Boolean upcall=0):
     omniCallDescriptor(lcfn, op_, oplen, 0, _user_exns, 0, upcall)
  {
    
  }
  
  void marshalArguments(cdrStream&);
  void unmarshalArguments(cdrStream&);

  void unmarshalReturnedValues(cdrStream&);
  void marshalReturnedValues(cdrStream&);
  
  
  static const char* const _user_exns[];

  ::CORBA::String_var arg_0_;
  const char* arg_0;
  ::CORBA::String_var arg_1_;
  const char* arg_1;
  StringArray_var arg_2_;
  const StringArray* arg_2;
  HepRepInstanceTree_var result;
};

void _0RL_cd_05fda001c74a4a06_40000000::marshalArguments(cdrStream& _n)
{
  _n.marshalString(arg_0,0);
  _n.marshalString(arg_1,0);
  (const StringArray&) *arg_2 >>= _n;

}

void _0RL_cd_05fda001c74a4a06_40000000::unmarshalArguments(cdrStream& _n)
{
  arg_0_ = _n.unmarshalString(0);
  arg_0 = arg_0_.in();
  arg_1_ = _n.unmarshalString(0);
  arg_1 = arg_1_.in();
  arg_2_ = new StringArray;
  (StringArray&)arg_2_ <<= _n;
  arg_2 = &arg_2_.in();

}

void _0RL_cd_05fda001c74a4a06_40000000::marshalReturnedValues(cdrStream& _n)
{
  (const HepRepInstanceTree&) result >>= _n;

}

void _0RL_cd_05fda001c74a4a06_40000000::unmarshalReturnedValues(cdrStream& _n)
{
  result = new HepRepInstanceTree;
  (HepRepInstanceTree&)result <<= _n;

}

const char* const _0RL_cd_05fda001c74a4a06_40000000::_user_exns[] = {
  0
};

// Local call call-back function.
static void
_0RL_lcfn_05fda001c74a4a06_50000000(omniCallDescriptor* cd, omniServant* svnt)
{
  _0RL_cd_05fda001c74a4a06_40000000* tcd = (_0RL_cd_05fda001c74a4a06_40000000*)cd;
  _impl_HepRep* impl = (_impl_HepRep*) svnt->_ptrToInterface(HepRep::_PD_repoId);
  tcd->result = impl->getInstances(tcd->arg_0, tcd->arg_1, *tcd->arg_2);


}

HepRepInstanceTree* _objref_HepRep::getInstances(const char* instanceTreeName, const char* instanceTreeVersion, const ::StringArray& typeNames)
{
  _0RL_cd_05fda001c74a4a06_40000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_50000000, "getInstances", 13);
  _call_desc.arg_0 = instanceTreeName;
  _call_desc.arg_1 = instanceTreeVersion;
  _call_desc.arg_2 = &(::StringArray&) typeNames;

  _invoke(_call_desc);
  return _call_desc.result._retn();


}
// Proxy call descriptor class. Mangled signature:
//  _cHepRepInstanceTree_i_cstring_i_cstring_i_cStringArray_i_cHepRepActionList_i_cboolean_i_cboolean_i_cboolean_i_cStringArray
class _0RL_cd_05fda001c74a4a06_60000000
  : public omniCallDescriptor
{
public:
  inline _0RL_cd_05fda001c74a4a06_60000000(LocalCallFn lcfn,const char* op_,size_t oplen,_CORBA_Boolean upcall=0):
     omniCallDescriptor(lcfn, op_, oplen, 0, _user_exns, 0, upcall)
  {
    
  }
  
  void marshalArguments(cdrStream&);
  void unmarshalArguments(cdrStream&);

  void unmarshalReturnedValues(cdrStream&);
  void marshalReturnedValues(cdrStream&);
  
  
  static const char* const _user_exns[];

  ::CORBA::String_var arg_0_;
  const char* arg_0;
  ::CORBA::String_var arg_1_;
  const char* arg_1;
  StringArray_var arg_2_;
  const StringArray* arg_2;
  HepRepActionList_var arg_3_;
  const HepRepActionList* arg_3;
  ::CORBA::Boolean arg_4;
  ::CORBA::Boolean arg_5;
  ::CORBA::Boolean arg_6;
  StringArray_var arg_7_;
  const StringArray* arg_7;
  HepRepInstanceTree_var result;
};

void _0RL_cd_05fda001c74a4a06_60000000::marshalArguments(cdrStream& _n)
{
  _n.marshalString(arg_0,0);
  _n.marshalString(arg_1,0);
  (const StringArray&) *arg_2 >>= _n;
  (const HepRepActionList&) *arg_3 >>= _n;
  _n.marshalBoolean(arg_4);
  _n.marshalBoolean(arg_5);
  _n.marshalBoolean(arg_6);
  (const StringArray&) *arg_7 >>= _n;

}

void _0RL_cd_05fda001c74a4a06_60000000::unmarshalArguments(cdrStream& _n)
{
  arg_0_ = _n.unmarshalString(0);
  arg_0 = arg_0_.in();
  arg_1_ = _n.unmarshalString(0);
  arg_1 = arg_1_.in();
  arg_2_ = new StringArray;
  (StringArray&)arg_2_ <<= _n;
  arg_2 = &arg_2_.in();
  arg_3_ = new HepRepActionList;
  (HepRepActionList&)arg_3_ <<= _n;
  arg_3 = &arg_3_.in();
  arg_4 = _n.unmarshalBoolean();
  arg_5 = _n.unmarshalBoolean();
  arg_6 = _n.unmarshalBoolean();
  arg_7_ = new StringArray;
  (StringArray&)arg_7_ <<= _n;
  arg_7 = &arg_7_.in();

}

void _0RL_cd_05fda001c74a4a06_60000000::marshalReturnedValues(cdrStream& _n)
{
  (const HepRepInstanceTree&) result >>= _n;

}

void _0RL_cd_05fda001c74a4a06_60000000::unmarshalReturnedValues(cdrStream& _n)
{
  result = new HepRepInstanceTree;
  (HepRepInstanceTree&)result <<= _n;

}

const char* const _0RL_cd_05fda001c74a4a06_60000000::_user_exns[] = {
  0
};

// Local call call-back function.
static void
_0RL_lcfn_05fda001c74a4a06_70000000(omniCallDescriptor* cd, omniServant* svnt)
{
  _0RL_cd_05fda001c74a4a06_60000000* tcd = (_0RL_cd_05fda001c74a4a06_60000000*)cd;
  _impl_HepRep* impl = (_impl_HepRep*) svnt->_ptrToInterface(HepRep::_PD_repoId);
  tcd->result = impl->getInstancesAfterAction(tcd->arg_0, tcd->arg_1, *tcd->arg_2, *tcd->arg_3, tcd->arg_4, tcd->arg_5, tcd->arg_6, *tcd->arg_7);


}

HepRepInstanceTree* _objref_HepRep::getInstancesAfterAction(const char* instanceTreeName, const char* instanceTreeVersion, const ::StringArray& typeNames, const ::HepRepActionList& actions, ::CORBA::Boolean getPoints, ::CORBA::Boolean getDrawAtts, ::CORBA::Boolean getNonDrawAtts, const ::StringArray& invertAtts)
{
  _0RL_cd_05fda001c74a4a06_60000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_70000000, "getInstancesAfterAction", 24);
  _call_desc.arg_0 = instanceTreeName;
  _call_desc.arg_1 = instanceTreeVersion;
  _call_desc.arg_2 = &(::StringArray&) typeNames;
  _call_desc.arg_3 = &(::HepRepActionList&) actions;
  _call_desc.arg_4 = getPoints;
  _call_desc.arg_5 = getDrawAtts;
  _call_desc.arg_6 = getNonDrawAtts;
  _call_desc.arg_7 = &(::StringArray&) invertAtts;

  _invoke(_call_desc);
  return _call_desc.result._retn();


}
// Proxy call descriptor class. Mangled signature:
//  _cStringArray
class _0RL_cd_05fda001c74a4a06_80000000
  : public omniCallDescriptor
{
public:
  inline _0RL_cd_05fda001c74a4a06_80000000(LocalCallFn lcfn,const char* op_,size_t oplen,_CORBA_Boolean upcall=0):
     omniCallDescriptor(lcfn, op_, oplen, 0, _user_exns, 0, upcall)
  {
    
  }
  
  
  void unmarshalReturnedValues(cdrStream&);
  void marshalReturnedValues(cdrStream&);
  
  
  static const char* const _user_exns[];

  StringArray_var result;
};

void _0RL_cd_05fda001c74a4a06_80000000::marshalReturnedValues(cdrStream& _n)
{
  (const StringArray&) result >>= _n;

}

void _0RL_cd_05fda001c74a4a06_80000000::unmarshalReturnedValues(cdrStream& _n)
{
  result = new StringArray;
  (StringArray&)result <<= _n;

}

const char* const _0RL_cd_05fda001c74a4a06_80000000::_user_exns[] = {
  0
};

// Local call call-back function.
static void
_0RL_lcfn_05fda001c74a4a06_90000000(omniCallDescriptor* cd, omniServant* svnt)
{
  _0RL_cd_05fda001c74a4a06_80000000* tcd = (_0RL_cd_05fda001c74a4a06_80000000*)cd;
  _impl_HepRep* impl = (_impl_HepRep*) svnt->_ptrToInterface(HepRep::_PD_repoId);
  tcd->result = impl->getLayerOrder();


}

StringArray* _objref_HepRep::getLayerOrder()
{
  _0RL_cd_05fda001c74a4a06_80000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_90000000, "getLayerOrder", 14);


  _invoke(_call_desc);
  return _call_desc.result._retn();


}
// Proxy call descriptor class. Mangled signature:
//  _cstring
class _0RL_cd_05fda001c74a4a06_a0000000
  : public omniCallDescriptor
{
public:
  inline _0RL_cd_05fda001c74a4a06_a0000000(LocalCallFn lcfn,const char* op_,size_t oplen,_CORBA_Boolean upcall=0):
     omniCallDescriptor(lcfn, op_, oplen, 0, _user_exns, 0, upcall)
  {
    
  }
  
  
  void unmarshalReturnedValues(cdrStream&);
  void marshalReturnedValues(cdrStream&);
  
  
  static const char* const _user_exns[];

  ::CORBA::String_var result;
};

void _0RL_cd_05fda001c74a4a06_a0000000::marshalReturnedValues(cdrStream& _n)
{
  _n.marshalString(result,0);

}

void _0RL_cd_05fda001c74a4a06_a0000000::unmarshalReturnedValues(cdrStream& _n)
{
  result = _n.unmarshalString(0);

}

const char* const _0RL_cd_05fda001c74a4a06_a0000000::_user_exns[] = {
  0
};

// Local call call-back function.
static void
_0RL_lcfn_05fda001c74a4a06_b0000000(omniCallDescriptor* cd, omniServant* svnt)
{
  _0RL_cd_05fda001c74a4a06_a0000000* tcd = (_0RL_cd_05fda001c74a4a06_a0000000*)cd;
  _impl_HepRep* impl = (_impl_HepRep*) svnt->_ptrToInterface(HepRep::_PD_repoId);
  tcd->result = impl->checkForException();


}

char* _objref_HepRep::checkForException()
{
  _0RL_cd_05fda001c74a4a06_a0000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_b0000000, "checkForException", 18);


  _invoke(_call_desc);
  return _call_desc.result._retn();


}
_pof_HepRep::~_pof_HepRep() {}


omniObjRef*
_pof_HepRep::newObjRef(omniIOR* ior, omniIdentity* id)
{
  return new ::_objref_HepRep(ior, id);
}


::CORBA::Boolean
_pof_HepRep::is_a(const char* id) const
{
  if( omni::ptrStrMatch(id, ::HepRep::_PD_repoId) )
    return 1;
  
  return 0;
}

const _pof_HepRep _the_pof_HepRep;

_impl_HepRep::~_impl_HepRep() {}


::CORBA::Boolean
_impl_HepRep::_dispatch(omniCallHandle& _handle)
{
  const char* op = _handle.operation_name();

  if( omni::strMatch(op, "getInstanceTreeTop") ) {

    _0RL_cd_05fda001c74a4a06_00000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_10000000, "getInstanceTreeTop", 19, 1);
    
    _handle.upcall(this,_call_desc);
    return 1;
  }

  if( omni::strMatch(op, "getTypeTree") ) {

    _0RL_cd_05fda001c74a4a06_20000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_30000000, "getTypeTree", 12, 1);
    
    _handle.upcall(this,_call_desc);
    return 1;
  }

  if( omni::strMatch(op, "getInstances") ) {

    _0RL_cd_05fda001c74a4a06_40000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_50000000, "getInstances", 13, 1);
    
    _handle.upcall(this,_call_desc);
    return 1;
  }

  if( omni::strMatch(op, "getInstancesAfterAction") ) {

    _0RL_cd_05fda001c74a4a06_60000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_70000000, "getInstancesAfterAction", 24, 1);
    
    _handle.upcall(this,_call_desc);
    return 1;
  }

  if( omni::strMatch(op, "getLayerOrder") ) {

    _0RL_cd_05fda001c74a4a06_80000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_90000000, "getLayerOrder", 14, 1);
    
    _handle.upcall(this,_call_desc);
    return 1;
  }

  if( omni::strMatch(op, "checkForException") ) {

    _0RL_cd_05fda001c74a4a06_a0000000 _call_desc(_0RL_lcfn_05fda001c74a4a06_b0000000, "checkForException", 18, 1);
    
    _handle.upcall(this,_call_desc);
    return 1;
  }


  return 0;
}

void*
_impl_HepRep::_ptrToInterface(const char* id)
{
  if( id == ::HepRep::_PD_repoId )
    return (::_impl_HepRep*) this;
  
  if( id == ::CORBA::Object::_PD_repoId )
    return (void*) 1;

  if( omni::strMatch(id, ::HepRep::_PD_repoId) )
    return (::_impl_HepRep*) this;
  
  if( omni::strMatch(id, ::CORBA::Object::_PD_repoId) )
    return (void*) 1;
  return 0;
}

const char*
_impl_HepRep::_mostDerivedRepoId()
{
  return ::HepRep::_PD_repoId;
}

POA_HepRep::~POA_HepRep() {}

