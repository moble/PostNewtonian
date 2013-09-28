PostNewtonian
=============

Sympy code to collect and process terms for various post-Newtonian
(PN) equations, and to generate working python/C/C++ code to evolve PN
systems

Introduction
============

The current state of PN literature is somewhat confusing to anyone
wishing to implement PN evolutions.  There is a huge array of sources,
each of which provides partial or outdated listings of the various PN
expressions.  For example, a recent paper describing a new spin-orbit
contribution to the gravitational-wave flux might give the complete
expression for all spin-orbit terms in the flux, but omit non-spinning
and spin-spin terms.  This leaves anyone wishing to collect the most
accurate and up-to-date PN expressions to search through a vast
literature, sorting out various conventions.  Even worse, this is
generally done in isolation by one person, whose work is then lost to
other people.

This project aims to correct that situation by providing a simple
framework for collecting PN expressions, combining them to
automatically calculate the various TaylorTn approximants, for
example, centralizing the results, and sharing them in a way that can
be used by as many people as possible -- for both analytical work and
computational work.  In particular, the computational side is intended
to support people working in pure python, C, C++, and Mathematica.

This repository contains python notebooks which collect various PN
expressions, allow for their symbolic manipulation using the , and
generate C/C++ code which allows for efficient evaluation.

- Collect various PN expressions, with annotations describing where
  they come from

- Manipulate and combine partial expressions from different sources
  into complete expressions appropriate for actual use

- Generate python code for simple (but possibly slow) evaluation of PN
  systems and waveforms with, e.g., precession effects, or
  neutron-star tidal effects, etc.

- Generate C/C++ code for more efficient evaluation of those systems

- Export expressions to Mathematica

- Export expressions to LaTeX


Contributing
============

Contributions are enthusiastically welcomed.  Github has [very useful
features](https://help.github.com/articles/be-social) for enabling
collaboration.  The basic process is to
[fork](https://help.github.com/articles/fork-a-repo) this repository,
make the changes in your fork, and then submit a [pull
request](https://help.github.com/articles/using-pull-requests) for the
author to pull changes back into this main repository.  This is easier
than it might sound.  The links listed here give more than enough
information to do this.


Getting Started
===============

The code in this directory is contained primarily in ipython notebooks
(which are much like Mathematica notebooks, but use python).  To start
it up, run

    ipython notebook --pylab=inline

Your browser should open automatically, and you should see a list of
notebooks in this directory.  (If not, see the
[Installation](#Installation) section below.)  Click one of those
notebooks, which should open in a new tab.  To run code, just put your
cursor in any code cell and hit Shift-Enter, as with Mathematica.


Installation
============

All of the code here uses python and various python packages (though
C/C++ code is generated).  So, you need an up-to-date installation of
python, as well as various python packages.  Most package managers
(apt, macports, homebrew, etc.) can install these packages for you.
For manual installation, install python and pip (python's package
manager), and then run

```Shell
pip install numpy
pip install matplotlib
pip install pandas
pip install sympy
pip install scipy
pip install ipython[notebook]
```

(For tcsh users, the brackets in the last command will need to be
escaped with backslashes.)


