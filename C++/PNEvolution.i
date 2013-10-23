// -*- c++ -*-

// Copyright (c) 2013, Michael Boyle
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

namespace std {
  %template(vectord) vector<double>;
};

// Return the values by reference as python
ARGOUT_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, t, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, v, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, chi1, NPY_DOUBLE)
ARGOUT_TYPEMAP_STD_VECTOR_OF_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, chi2, NPY_DOUBLE)
%typemap (in,numinputs=0) std::vector<Quaternions::Quaternion>& R_frame (std::vector<Quaternions::Quaternion> vec_temp) {
  $1 = &vec_temp;
}
%typemap(argout) std::vector<Quaternions::Quaternion>& R_frame {
  npy_intp size = $1->size();
  PyArrayObject *npy_arr = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNew(1, &size, NPY_OBJECT));
  PyObject** data = static_cast<PyObject**>(PyArray_DATA(npy_arr));
  for(npy_intp i=0; i<size; ++i) {
    PyObject* qobj = SWIG_NewPointerObj((new Quaternions::Quaternion((*$1)[i])),
					SWIGTYPE_p_Quaternions__Quaternion, SWIG_POINTER_OWN);
    if(!qobj) {SWIG_fail;}
    Py_INCREF(qobj);
    data[i] = qobj;
  }
  %append_output(PyArray_Return(npy_arr));
}
ARGOUT_TYPEMAP_STD_VECTOR_OF_PRIMITIVES(double, DOUBLE, Phi, NPY_DOUBLE)



/////////////////////////////////////
//// Import the quaternion class ////
/////////////////////////////////////
// %ignore Quaternions::Quaternion::operator=;
%include "PNEvolution.hpp"


/// Add utility functions that are specific to python.  Note that
/// these are defined in the Quaternions namespace.
%insert("python") %{

%}
