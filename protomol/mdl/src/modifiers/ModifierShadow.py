# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _ModifierShadow

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types


class Modifier(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Modifier, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Modifier, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C++ ProtoMol::Modifier instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __del__(self, destroy=_ModifierShadow.delete_Modifier):
        try:
            if self.thisown: destroy(self)
        except: pass

    def execute(*args): return _ModifierShadow.Modifier_execute(*args)
    def isInternal(*args): return _ModifierShadow.Modifier_isInternal(*args)
    def order(*args): return _ModifierShadow.Modifier_order(*args)
    def enable(*args): return _ModifierShadow.Modifier_enable(*args)
    def disable(*args): return _ModifierShadow.Modifier_disable(*args)
    def isEnabled(*args): return _ModifierShadow.Modifier_isEnabled(*args)
    def __lt__(*args): return _ModifierShadow.Modifier___lt__(*args)
    def initialize(*args): return _ModifierShadow.Modifier_initialize(*args)
    def ProtoMol_Modifier_print(*args): return _ModifierShadow.Modifier_ProtoMol_Modifier_print(*args)

class ModifierPtr(Modifier):
    def __init__(self, this):
        _swig_setattr(self, Modifier, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, Modifier, 'thisown', 0)
        _swig_setattr(self, Modifier,self.__class__,Modifier)
_ModifierShadow.Modifier_swigregister(ModifierPtr)

class ModifierShadow(Modifier):
    __swig_setmethods__ = {}
    for _s in [Modifier]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, ModifierShadow, name, value)
    __swig_getmethods__ = {}
    for _s in [Modifier]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, ModifierShadow, name)
    def __repr__(self):
        return "<%s.%s; proxy of C++ ProtoMol::ModifierShadow instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    def __init__(self, *args):
        _swig_setattr(self, ModifierShadow, 'this', _ModifierShadow.new_ModifierShadow(*args))
        _swig_setattr(self, ModifierShadow, 'thisown', 1)
    def isInternal(*args): return _ModifierShadow.ModifierShadow_isInternal(*args)
    def calcShadow(*args): return _ModifierShadow.ModifierShadow_calcShadow(*args)
    def resetHistory(*args): return _ModifierShadow.ModifierShadow_resetHistory(*args)
    def __del__(self, destroy=_ModifierShadow.delete_ModifierShadow):
        try:
            if self.thisown: destroy(self)
        except: pass


class ModifierShadowPtr(ModifierShadow):
    def __init__(self, this):
        _swig_setattr(self, ModifierShadow, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, ModifierShadow, 'thisown', 0)
        _swig_setattr(self, ModifierShadow,self.__class__,ModifierShadow)
_ModifierShadow.ModifierShadow_swigregister(ModifierShadowPtr)


