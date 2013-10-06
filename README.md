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


Compiling with other code
=========================

This main library is fairly simple, and should be trivial to compile;
only the header `Quaternions.hpp` needs to be included, and only the
file `Quaternions.cpp` needs to be compiled.

The second pair of files, `IntegrateAngularVelocity.{ch}pp` contains
just one function, but depends on GSL for ODE integration.  So GSL
must be installed separately, compiled as a shared library, and the
`-I` and `-L` flags variables set appropriately on whatever
compilation is done.

For python, compilation is done automatically, and assumes the
presence of GSL.


Installing the python module
============================

Though this code can be included as a library in other code, it can
also be used on its own as a python module.  Just run

    python setup.py install --user

The `--user` flag installs the module to the user's home directory,
which means that no root permissions are necessary.

As mentioned above, GSL must be installed as a shared library, and the
`IncDirs` and `LibDirs` variables set appropriately in `setup.py`.
Sensible defaults are given, so this may only need to be done if the
compilation fails.

If the build succeeds, just open an python session and type

    import Quaternions

In ipython, you can then see your options using tab completion by typing

    Quaternions.

and then hitting tab.  Help is available on many functions in ipython
by typing a question mark after the function name.  For example:

    Quaternions.Squad?

In plain python, the same thing can be achieved by entering
`help(Quaternions.Squad)`.  But if you're using plain python
interactively, you should really give ipython a try.
