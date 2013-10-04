Quaternions
===========

Quaternion library for C++, with python bindings via SWIG.

This code provides a simple class for calculating easily with
quaternions.  Operator overloads are provided for all the basic
operations such as addition, multiplication, etc.  The standard
functions like `exp` and `log` (which are what make quaternions so
uniquely powerful) are also available.

Additionally, support is included for time series of quaternions, as
well as operations such as differentiation, interpolation (linear and
quadratic).  Several algorithms also provide capabilities for finding
the minimal-rotation frame that tracks a certain vector, or finds the
frame that has a particular angular velocity.



Installing the python module
============================

Though this code can be included as a library in other code, it can
also be used on its own as a python module.  Just run

    python setup.py install --user

The `--user` flag installs the module to the user's home directory,
which means that no root permissions are necessary.  If this succeeds,
just open an python session and type

    import Quaternions

In ipython, you can then see your options using tab completion by typing

    Quaternions.

and then hitting tab.  Help is available on many functions in ipython
by typing a question mark after the function name.  For example:

    Quaternions.Squad?

In plain python, the same thing can be achieved by entering
`help(Quaternions.Squad)`.
