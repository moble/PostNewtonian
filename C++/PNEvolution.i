// -*- c++ -*-
// vim: filetype=cpp

// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

%module PNEvolution

// Quiet warnings about overloaded operators being ignored.
#pragma SWIG nowarn=362,389,401,509

// Unless otherwise provided, just let SWIG document any functions.
%feature("autodoc", "1");

#ifndef SWIGIMPORTED

////////////////////////////////////////////////////////
//// This lets me use numpy.array in the code below ////
////////////////////////////////////////////////////////
%{
  #define SWIG_FILE_WITH_INIT
%}
%include "Quaternions/numpy.i"
%init %{
  import_array();
%}


///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////

// The following will appear in the header of the `_wrap.cpp` file.
%{
  const char* const PostNewtonianErrors[] = {
    "This function is not yet implemented.",
    "Failed system call.",
    "Bad file name.",
    "Failed GSL call.",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Bad value.",
    "Bad switches; we should not have gotten here.",
    "Index out of bounds.",
    "Unknown exception",
    "Unknown exception",
    "Vector size mismatch.",
    "Matrix size mismatch.",
    "Matrix size is assumed to be 3x3 in this function.",
    "Not enough points to take a derivative.",
    "Empty intersection requested.",
    "Waveform is missing requested (ell,m) component.",
    "Wrong frame type for this operation.",
    "Bad Waveform information."
  };
  const int PostNewtonianNumberOfErrors = 23;
  PyObject* const PostNewtonianExceptions[] = {
    PyExc_NotImplementedError, // Not implemented
    PyExc_SystemError, // Failed system call
    PyExc_IOError, // Bad file name
    PyExc_RuntimeError, // GSL failed
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_ValueError, // Bad value
    PyExc_ValueError, // Bad switches
    PyExc_IndexError, // Index out of bounds
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_AssertionError, // Mismatched vector size
    PyExc_AssertionError, // Mismatched matrix size
    PyExc_AssertionError, // 3x3 matrix assumed
    PyExc_AssertionError, // Not enough points for derivative
    PyExc_AssertionError, // Empty intersection
    PyExc_IndexError, // Waveform missing ell,m
    PyExc_AssertionError, // Bad frame type
    PyExc_ValueError, // Bad Waveform information
  };
%}

// This will go inside every python wrapper for any function I've
// included; the code of the function itself will replace `$action`.
// It's a good idea to try to keep this part brief, just to cut down
// the size of the wrapper file.
%exception {
  try {
    $action;
  } catch(int i) {
    std::stringstream s;
    if(i>-1 && i<PostNewtonianNumberOfErrors) { s << "PostNewtonian exception: " << PostNewtonianErrors[i]; }
    else  { s << "PostNewtonian: Unknown exception number {" << i << "}"; }
    PyErr_SetString(PostNewtonianExceptions[i], s.str().c_str());
    return NULL;
  }
}

#endif // SWIGIMPORTED


/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
  #include <vector>
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <iomanip>
  #include <complex>
  #include "Quaternions.hpp"
  #include "IntegrateAngularVelocity.hpp"
  #include "PNEvolution.hpp"
  #include "PNWaveformModes.hpp"
%}


//////////////////////////////////////////////////////////////////////
//// The following translates between c++ and python types nicely ////
//////////////////////////////////////////////////////////////////////

%include "Quaternions/Quaternions.i"
// %include "Quaternions/Quaternions_typemaps.i"
%include "Quaternions/vector_typemaps.i"

%include <typemaps.i>
%include <stl.i>
//// Make sure std::strings are dealt with appropriately
%include <std_string.i>
//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>



//////////////////////////////////////////////////////////////////////////////////
//// This ensures that PNEvolution has local copies of these modules imported ////
//////////////////////////////////////////////////////////////////////////////////
%pythoncode %{
  ## We have to be able to import numpy
  import numpy;

  ## We have to be able to import Quaternions
  import Quaternions
%}


///////////////////////////////
//// Import the PN classes ////
///////////////////////////////
%include "PNEvolution.hpp"
%include "PNWaveformModes.hpp"


/// Add utility functions that are specific to python here.  Note that
/// these are defined in the PNEvolution namespace.
%insert("python") %{

%}
