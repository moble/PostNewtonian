PostNewtonian
=============

Sympy code to collect and process terms for various post-Newtonian
(PN) equations, and to generate working python/C/C++ code to evolve PN
systems


Quick Start
===========

From the command line, change your directory to this code directory
(with all the `.ipynb` files), and run

    ipython notebook --pylab=inline

Your browser should open automatically, and you should see a list of
notebooks.  (If not, see the [Installation](#Installation) section
below.)  Click one of those notebooks, which should open in a new tab.
To run code, just put your cursor in any code cell and hit
Shift-Enter, as with Mathematica.

Introduction
============

The aim of this project is to provide a simple centralized framework
for collecting PN expressions, combining them to automatically
calculate, e.g., TaylorTn expressions, and sharing them in a way that
can be used by as many people as possible---through LaTeX, C, C++, or
Mathematica code.

To be more explicit, the tasks performed by code in this module will:

- Collect various PN expressions, with annotations describing where
  they come from

- Manipulate and combine partial expressions from different sources
  into complete expressions appropriate for actual use

- Generate python code for simple (but possibly slow) evaluation of PN
  systems and waveforms with, e.g., precession effects, or
  neutron-star tidal effects, etc.

- Export C/C++ code for more efficient evaluation of those systems

- Export expressions to Mathematica

- Export expressions to LaTeX

One important feature of the `ipython notebook` is that it allows easy
description of code, meaning that we can give very explicit citations
and other explanations for where to find the collected terms in the
literature.  The notebooks depend on each other, and build up more
complicated expressions.  The rough order of this dependency is as
follows.

- `Variables.ipynb`: Define the fundamental variables, and write all
  the non-fundamental variables in terms of them.  Centralizing these
  definitions reduces mistakes.

- `EnergyAbsorption.ipynb`, `Flux.ipynb`, `OrbitalEnergy.ipynb`:
  Collect the PN expressions for these various quantities, classified
  by their type and PN order.

- `Precession.ipynb`: Define the precession system.

- `OrbitalFrequencyEvolution.ipynb`: Derive the TaylorTn expressions
  from the three notebooks above, and generate C/C++ code to
  evolve or evaluate them.

- [More to follow...]


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


Installation
============

All of the code here uses python and various python packages (though
C/C++ code is generated).  So, you need an up-to-date installation of
python (version 2.7 or greater), as well as various python packages.
Most package managers (apt, macports, homebrew, etc.) can install
these packages for you.  For manual installation, install python and
pip (python's package manager), and then run

```Shell
pip install numpy
pip install matplotlib
pip install pandas
pip install sympy
pip install scipy
pip install ipython[notebook]
```

If these fail, complaining about permissions, simply add `--user` to
each command.  For tcsh users, the brackets in the last command will
need to be escaped with backslashes.


