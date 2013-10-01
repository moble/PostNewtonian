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

# TaylorT1 and TaylorT4

# <markdowncell>

# These two very similar approximants are the simplest in construction, and most widely applicable.  In particular, they can both be applied to precessing systems.  Each gives rise to the same system of ODEs that need to be integrated in time, except that the right-hand side for $dv/dt$ is expanded as a series in $v$ and truncated for TaylorT4.

# <headingcell level=1>

# TaylorT2 and TaylorT3

# <markdowncell>

# These two approximants are also closely related to each other.

# <codecell>


