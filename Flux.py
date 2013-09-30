# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Always run this first
# NOTE: Do not define new variables in this notebook or in Variables.py; 
#       define them in Variables.ipynb.  Variables.py will be automatically
#       recreated when Variables.ipynb is saved.

%run -i Variables.py
from __future__ import division # This needs to be here, even though it's in Variables.py

# <headingcell level=1>

# Individual flux terms

# <markdowncell>

# In this notebook, every term will be multiplied by the following coefficient.

# <codecell>

FluxCoefficient = (32/5)*nu**2*v**10

# <markdowncell>

# We will define a dictionary of dictionaries.  The top-level keys will be different types of terms.  The dictionaries they contain will give those terms at different PN orders. 

# <codecell>

FluxTerms = {'Nonspinning':{},
             'SpinSpin':{},
             'SpinOrbit':{},
             'NSTides':{},
             }

# <codecell>

FluxTerms['Nonspinning'][0] = 1
FluxTerms['Nonspinning'][1] = 0
FluxTerms['Nonspinning'][2] = -(1247/336) - (35/12)*nu
FluxTerms['Nonspinning'][3] = 4*pi
FluxTerms['Nonspinning'][4] = 0
FluxTerms['Nonspinning'][5] = 0
FluxTerms['Nonspinning'][6] = 0
FluxTerms['Nonspinning'][7] = 0
FluxTerms['Nonspinning'][8] = 0
FluxTerms['Nonspinning'][9] = 0

# <codecell>

print FluxTerms

# <headingcell level=1>

# Collected flux terms

# <codecell>

def Flux(SpinTerms=True, PrecessingSpinTerms=False, NSTidalTerms=False, EMRITerms=False) :
    """
    Return an expression for the GW flux with the given options.
    
    """
    FluxSum = 0
    
    return FluxCoefficient*FluxSum

