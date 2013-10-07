#! /usr/bin/env python

"""
Installation file for python code associated with the paper "Angular
velocity of gravitational radiation from precessing binaries and the
corotating frame".

To build this code and run it in place, run
    python setup.py build_ext --inplace
then open python and type 'import Quaternions' in
the current directory.

To install in the user's directory, run
    python setup.py install --user
Now, 'import Quaternions' may be run from a python
instance started in any directory on the system.
"""

## Check for `--no-GSL` option
from sys import argv
if '--no-GSL' in argv:
    GSL=False
    argv.remove('--no-GSL')
else:
    GSL=True

## If PRD won't let me keep a subdirectory, make one
from os.path import exists
from os import makedirs
if not exists('Quaternions') :
    makedirs('Quaternions')
# from shutil import copyfile
# if not exists('Quaternions/plot.py') :
#     copyfile('plot.py', 'Quaternions/plot.py')

## distutils doesn't build swig modules in the correct order by
## default -- the python module is installed first.  This will pop
## 'build_ext' to the beginning of the command list.
from distutils.command.build import build
build.sub_commands = sorted(build.sub_commands, key=lambda sub_command: int(sub_command[0]!='build_ext'))

## We also need to copy the SWIG-generated python script Quaternions.py
## to Quaternions/__init__.py so that it gets installed correctly.
from distutils.command.build_ext import build_ext as _build_ext
from distutils.file_util import copy_file
class build_ext(_build_ext):
    """Specialized Python source builder for moving SWIG module."""
    def run(self):
        _build_ext.run(self)
        copy_file('Quaternions.py', 'Quaternions/__init__.py')

## Now import the basics
from distutils.core import setup, Extension
from subprocess import check_output, CalledProcessError
from os import devnull, environ

## See if GSL_HOME is set; if so, use it
if GSL :
    SourceFiles = ['Quaternions.cpp',
                   'IntegrateAngularVelocity.cpp',
                   'Quaternions.i']
    Dependencies = ['Quaternions.hpp',
                    'IntegrateAngularVelocity.hpp']
    Libraries = ['gsl', 'gslcblas']
    if "GSL_HOME" in environ :
        IncDirs = [environ["GSL_HOME"]+'/include', '/opt/local/include']
        LibDirs = [environ["GSL_HOME"]+'/lib', '/opt/local/lib']
    else :
        IncDirs = ['/opt/local/include']
        LibDirs = ['/opt/local/lib']
else :
    SourceFiles = ['Quaternions.cpp',
                   'Quaternions.i']
    Dependencies = ['Quaternions.hpp']
    Libraries = []
    IncDirs = []
    LibDirs = []

## Remove a compiler flag that doesn't belong there for C++
import distutils.sysconfig as ds
cfs=ds.get_config_vars()
for key, value in cfs.iteritems() :
    if(type(cfs[key])==str) :
        cfs[key] = value.replace('-Wstrict-prototypes', '')

# ## Try to determine an automatic version number for this
# try :
#     with open(devnull, "w") as fnull :
#         GitRev = check_output('git rev-parse HEAD', shell=True, stderr=fnull)[:-1]
#         CodeRevision = '"{0}"'.format(GitRev)
#         PackageVersion = GitRev[:9]
# except (NameError, CalledProcessError) :
#     CodeRevision = '"PaperVersion3"'
#     PackageVersion = '3'

## Read in the license
try :
    with open('LICENSE', 'r') as myfile :
        License=myfile.read()
except IOError :
    License = 'See LICENSE file in the source code for details.'

## This does the actual work
setup(name="Quaternions",
      # version=PackageVersion,
      description='Quaternion library for C++, with python bindings via SWIG.',
      #long_description=""" """,
      author='Michael Boyle',
      author_email='boyle@astro.cornell.edu',
      url='https://github.com/MOBle',
      license=License,
      packages = ['Quaternions'],
      # py_modules = ['Quaternions'],
      # scripts = ['Scripts/RunExtrapolations.py', 'Scripts/ConvertGWDatToH5.py'],
      ext_modules = [
        Extension('_Quaternions',
                  sources=SourceFiles,
                  depends=Dependencies,
                  include_dirs=IncDirs,
                  library_dirs=LibDirs,
                  libraries=Libraries,
                  #define_macros = [('CodeRevision', CodeRevision)],
                  language='c++',
                  swig_opts=['-globals', 'constants', '-c++'],
                  extra_link_args=['-lgomp', '-fPIC'],
                  # extra_link_args=['-lgomp', '-fPIC', '-Wl,-undefined,error'], # `-undefined,error` tells the linker to fail on undefined symbols
                  extra_compile_args=['-Wno-deprecated']
                  # extra_compile_args=['-fopenmp', '-ffast-math'] # DON'T USE fast-math!!!  It makes it impossible to detect NANs
                  )
        ],
      # classifiers = ,
      # distclass = ,
      # script_name = ,
      # script_args = ,
      # options = ,
      # license = ,
      # keywords = ,
      # platforms = ,
      # cmdclass = ,
      cmdclass={'build_ext': build_ext},
      # data_files = ,
      # package_dir =
      )
