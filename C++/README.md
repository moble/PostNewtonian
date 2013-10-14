Collect C++ files for constructing PN orbits and waveforms

This subdirectory contains C++ files created automatically by the main
ipython notebooks in higher directories.  Compilation should be fairly
straightforward, as long as the two dependencies below are satisfied.


Dependencies
============

There are two dependencies for the C++ code: GSL, and the
`Quaternions` package.

GSL
---

The [GNU Scientific Library](http://www.gnu.org/software/gsl/), is
used for ODE integrations.  This needs to be compiled separately
(preferably as a shared library), and the libraries and headers need
to be accessible to the `Makefile`.

Quaternions
-----------

The submodule `Quaternions` refers to another git repository hosted
[here](https://github.com/MOBle/Quaternions), which should exist as a
subdirectory of this directory.  It may already be present, in which
case you do not need to do anything more.  However, if the directory
is empty, you will need to run two commands:

    git submodule init
    git submodule update

If you create new code that is compiled outside of this directory, you
will need to ensure that the header and code files can be found by
your compilation process.
