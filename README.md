PostNewtonian
=============

Sympy code to collect terms for various PN equations, and to generate
working C++ code to evolve PN systems


Getting Started
===============

The code in this directory is contained primarily in ipython notebooks
(which are much like Mathematica notebooks, but use python).  To start
it up, run

    ipython notebook --pylab=inline

You should see a list of notebooks in this directory.  (If not, see
the [Installation](#Installation) section below.)


Installation
============

All of the code here uses python and various python packages (though
C/C++ code is generated).  So, you need an up-to-date installation of
python, as well as various python packages.  Most package managers
(apt, macports, homebrew, etc.) can install these packages for you.
For manual installation, install python and pip (python's package
manager), and then run

    pip install numpy
    pip install matplotlib
    pip install pandas
    pip install sympy
    pip install scipy
    pip install ipython[notebook]

(For tcsh users, the brackets in the last command will need to be
escaped with backslashes.)
