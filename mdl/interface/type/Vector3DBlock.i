%module Vector3DBlock
%{
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include "ndarrayobject.h"
using namespace ProtoMol;
%}

%include <protomol/type/Vector3DBlock.h>

%extend ProtoMol::Vector3DBlock {


   Real __getitem__(int index){
	return self->c[index];
   }

   void __setitem__(int index, Real val) {
	self->c[index] = val;
   }


   void printC() {
      cout << "C: " << self->c << " MEMBER C: " << (*self)[0].c << endl;
   }

   void setC(PyObject* rhs) {
      import_array();
      self->resize(((PyArrayObject*)rhs)->dimensions[0]/3);
      if (!(PyArrayObject*)(self->c)) {
	delete self->c;
      }
      self->c = (Real*)(((PyArrayObject*)rhs)->data); 
   }

   PyObject* getC() {
	int mySize = self->size()*3;
        npy_intp dims[1] = {mySize};
	import_array1(NULL);
	PyObject* rhs = PyArray_SimpleNewFromData(1,dims,PyArray_DOUBLE,(char*)(self->c));
        return rhs;
   }

};