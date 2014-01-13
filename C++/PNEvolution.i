// -*- c++ -*-

// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

%module PNEvolution

// Quiet warnings about overloaded operators being ignored.
#pragma SWIG nowarn=362,389,401,509


/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
  #define SWIG_FILE_WITH_INIT
  #include <vector>
  #include <iostream>
  #include <iomanip>
  #include <complex>
  #include "Quaternions.hpp"
  #include "IntegrateAngularVelocity.hpp"
  #include "PNEvolution.hpp"
  #include "PNWaveformModes.hpp"
%}

///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////
%exception {
  try {
    $action;
  } catch(int i) {
    if(i==0) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Index out of bounds.");
    } else if(i==1) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Infinitely many solutions.");
    } else if(i==2) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Not enough points to take a derivative.");
    } else if(i==3) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Vector size not understood.");
    } else if(i==4) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Vector size inconsistent with another vector's size.");
    } else if(i==5) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Cannot extrapolate quaternions.");
    } else if(i==6) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Failed call to GSL.");
    } else if(i==7) {
      PyErr_SetString(PyExc_RuntimeError, "Quaternions: Unknown exception.");
    } else  {
      PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
    }
    return NULL;
  }
}


%include "Quaternions_typemaps.i"
%include "Quaternions.hpp"
// %include "IntegrateAngularVelocity.hpp"
// #if defined(SWIGPYTHON_BUILTIN)
// %feature("python:slot", "sq_length", functype="lenfunc") Quaternions::Quaternion::__len__;
// %feature("python:slot", "mp_subscript", functype="binaryfunc") Quaternions::Quaternion::__getitem__;
// %feature("python:slot", "mp_ass_subscript", functype="objobjargproc") Quaternions::Quaternion::__setitem__;
// %feature("python:slot", "tp_str",  functype="reprfunc") Quaternions::Quaternion::__str__;
// %feature("python:slot", "tp_repr", functype="reprfunc") Quaternions::Quaternion::__repr__;
// #endif // SWIGPYTHON_BUILTIN
// %extend Quaternions::Quaternion {
//   unsigned int __len__() const {
//     return 4;
//   }
//   inline double __getitem__(const unsigned int i) const {
//     return (*$self)[i];
//   }
//   inline void __setitem__(const unsigned int i, const double a) {
//     (*$self)[i] = a;
//   }
//   const char* __str__() {
//     std::stringstream S;
//     S << std::setprecision(15) << "["
//       << $self->operator[](0) << ", "
//       << $self->operator[](1) << ", "
//       << $self->operator[](2) << ", "
//       << $self->operator[](3) << "]";
//     const std::string& tmp = S.str();
//     const char* cstr = tmp.c_str();
//     return cstr;
//   }
//   const char* __repr__() {
//     std::stringstream S;
//     S << std::setprecision(15) << "Quaternion("
//       << $self->operator[](0) << ", "
//       << $self->operator[](1) << ", "
//       << $self->operator[](2) << ", "
//       << $self->operator[](3) << ")";
//     const std::string& tmp = S.str();
//     const char* cstr = tmp.c_str();
//     return cstr;
//   }
//  };


// Return the values by reference as python
ARGOUT_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, t, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, v, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, chi1, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, chi2, NPY_DOUBLE)
%apply std::vector<Quaternions::Quaternion>& Quaternion_argout { std::vector<Quaternions::Quaternion>& R_frame };
ARGOUT_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, Phi, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, L, NPY_DOUBLE)
OUT_TYPEMAP_STD_VECTOR_OF_STD_VECTOR_OF_PRIMITIVES(double, NPY_DOUBLE)

//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>
// namespace std {
//   %template(complexd) complex<double>; // Don't use this line!!!
// };
//// Make sure std::vectors are dealt with appropriately
%include <std_vector.i>
namespace std {
  %template(vectorc) vector<std::complex<double> >;
  %template(vectorvectorc) vector<vector<std::complex<double> > >;
};


//////////////////////////////////////
//// Import the PNEvolution class ////
//////////////////////////////////////
%include "PNEvolution.hpp"
%include "PNWaveformModes.hpp"


/// Add utility functions that are specific to python.  Note that
/// these are defined in the Quaternions namespace.
%insert("python") %{

%}
