// -*- c++ -*-

// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details


%module Quaternions

 // Quiet warnings about overloaded operators being ignored.
#pragma SWIG nowarn=362,389,401,509
%include <typemaps.i>
%include <stl.i>

// %{
//   #include "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
//   #include <vector>
// %}

// %init %{
//   import_array();
// %}

// %typemap(out) std::vector<int> {
//   npy_intp result_size = $1.size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = $1[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<int>& {
//   npy_intp result_size = $1->size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = (*$1)[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<std::vector<int> >& {
//   npy_intp result_size = $1->size();
//   npy_intp result_size2 = (result_size>0 ? (*$1)[0].size() : 0);
//   npy_intp dims[2] = { result_size, result_size2 };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { for (size_t j = 0; j < result_size2; ++j) { dat[i*result_size2+j] = (*$1)[i][j]; } }
//   $result = PyArray_Return(npy_arr);
// }

// %typemap(out) std::vector<double> {
//   npy_intp result_size = $1.size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = $1[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<double>& {
//   npy_intp result_size = $1->size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = (*$1)[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<std::vector<double> >& {
//   npy_intp result_size = $1->size();
//   npy_intp result_size2 = (result_size>0 ? (*$1)[0].size() : 0);
//   npy_intp dims[2] = { result_size, result_size2 };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { for (size_t j = 0; j < result_size2; ++j) { dat[i*result_size2+j] = (*$1)[i][j]; } }
//   $result = PyArray_Return(npy_arr);
// }

%include "docs/Quaternions_Doc.i"

///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////
%exception {
  try {
    $action;
  } catch(int i) {
    if(i==0) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Not yet implemented.");
    } else if(i==1) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Index out of bounds.");
    } else if(i==2) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Infinitely many solutions.");
    } else if(i==3) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Input vector size mismatch.");
    } else if(i==4) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Cannot extrapolate quaternions.");
    } else if(i==5) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Matrix size mismatch.");
    } else if(i==6) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Matrix size is assumed to be 3x3 in this function.");
    } else if(i==7) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Quaternion constructor's vector size not understood; should be 3 or 4.");
    } else if(i==8) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Waveform is missing requested l,m component.");
    } else if(i==9) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Bad file name.");
    } else if(i==10) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Not enough points to take a derivative.");
    } else if(i==11) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Empty intersection requested.");
    } else if(i==12) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Failed system call.");
    } else if(i==13) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Wrong FrameType for this operation.  Maybe you forgot to `SetFrameType`?");
    } else if(i==14) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: GSL failed.");
    } else if(i==15) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Bad Waveform information.");
    } else if(i==16) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Bad switches; we should not have gotten here.");
    } else if(i==17) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Bad value.");
    } else  {
      PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
    }
    return NULL;
  }
}

/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <iomanip>
  #include <complex>
  #include "Quaternions.hpp"
%}


//////////////////////////////////////////////////////////////////////
//// The following translates between c++ and python types nicely ////
//////////////////////////////////////////////////////////////////////
//// This lets me use numpy.array in the code below
%pythoncode %{
  import numpy;
  %}
//// Make sure std::strings are dealt with appropriately
%include <std_string.i>
//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>
// namespace std {
//   %template(complexd) complex<double>; // Don't use this line!!!
// };
//// Make sure std::vectors are dealt with appropriately
%include <std_vector.i>
namespace Quaternions {
  class Quaternion;
 };
namespace std {
  %template(vectori) vector<int>;
  %template(vectorvectori) vector<vector<int> >;
  %template(vectord) vector<double>;
  %template(vectorvectord) vector<vector<double> >;
  %template(vectorc) vector<std::complex<double> >;
  %template(vectorvectorc) vector<vector<std::complex<double> > >;
  %template(vectorq) vector<Quaternions::Quaternion>;
  %template(vectors) vector<string>;
  %template(vectorvectors) vector<vector<std::string> >;
};


/////////////////////////////////////
//// Import the quaternion class ////
/////////////////////////////////////
%ignore Quaternions::Quaternion::operator=;
%rename(__getitem__) Quaternions::Quaternion::operator [](const unsigned int) const;
%rename(__setitem__) Quaternions::Quaternion::operator [](const unsigned int);
%include "Quaternions.hpp"
%extend Quaternions::Quaternion {
  // This function is called when printing a Quaternion object
  const char* __str__() {
    std::stringstream S;
    S << std::setprecision(14) << "["
      << $self->operator[](0) << ", "
      << $self->operator[](1) << ", "
      << $self->operator[](2) << ", " 
      << $self->operator[](3) << "]";
    const std::string& tmp = S.str();
    const char* cstr = tmp.c_str();
    return cstr;
  }
  // This prints the Quaternion nicely at the prompt and allows nucer manipulations
  %pythoncode{
    def __repr__(self):
        return 'Quaternion('+repr(self[0])+', '+repr(self[1])+', '+repr(self[2])+', '+repr(self[3])+')'
    def __pow__(self, P) :
        return self.pow(P)
    __radd__ = __add__
    def __rsub__(self, t) :
        return -self+t
    __rmul__ = __mul__
    def __rdiv__(self, t) :
        return self.inverse()*t
  };
 };


/// Add utility functions that are specific to python.  Note that
/// these are defined in the Quaternions namespace.
%insert("python") %{


%}
