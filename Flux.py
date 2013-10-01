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

# In this notebook, every term will be multiplied by the following coefficient.  [Note that fractions need to be entered as, e.g., `frac(32,5)` so that they are not converted to finite-precision decimals.]

# <codecell>

FluxCoefficient = frac(32,5)*nu**2*v**10

# <markdowncell>

# We will define a dictionary of dictionaries to contain all the flux terms.  The top-level keys will be different types of terms: nonspinning, spin-spin, etc.  The dictionaries they contain will give those terms at different PN orders. 

# <codecell>

FluxTerms = {'Nonspinning':{},
             'ExtraEMRI':{},
             'SpinSpin':{},
             'SpinOrbit':{},
             'NSTidal':{},
             }

# <markdowncell>

# The nonspinning flux terms are complete to 3.5pN order.  These terms are given by Eq. (231) of [Blanchet (2006)](http://www.livingreviews.org/lrr-2006-4):

# <codecell>

FluxTerms['Nonspinning'][0] = 1
FluxTerms['Nonspinning'][1] = 0
FluxTerms['Nonspinning'][2] = -frac(1247,336) - frac(35,12)*nu
FluxTerms['Nonspinning'][3] = 4*pi
FluxTerms['Nonspinning'][4] = -frac(44711,9072) + frac(9271,504)*nu + frac(65,18)*nu**2
FluxTerms['Nonspinning'][5] = (-frac(8191,672) - frac(583,24)*nu)*pi
FluxTerms['Nonspinning'][6] = frac(6643739519,69854400) + frac(16,3)*pi**2 - EulerGamma*frac(1712,105) \
                              - frac(1712,105)*ln(4*v) + (-frac(134543,7776) + frac(41,48)*pi**2)*nu \
                              - frac(94403,3024)*nu**2 - frac(775,324)*nu**3
FluxTerms['Nonspinning'][7] = (-frac(16285,504) + frac(214745,1728)*nu + frac(193385,3024)*nu**2)*pi

# <markdowncell>

# This term is the 4.0pN nonspinning term in the limit as the mass ratio goes to zero.  This is given by Eq. (43) of [Tagoshi and Sasaki (1994)](http://arxiv.org/abs/gr-qc/9405062v1):

# <codecell>

FluxTerms['ExtraEMRI'][8] = -frac(323105549467,3178375200) - frac(1369,126)*pi**2 + frac(232597,4410)*EulerGamma \
                            + frac(39931,294)*ln(2) - frac(47385,1568)*ln(3) + frac(232597,4410)*ln(v)

# <markdowncell>

# The spin-spin terms in the flux are known to ...

# <codecell>

FluxTerms['SpinSpin'][0] = 
FluxTerms['SpinSpin'][1] = 
FluxTerms['SpinSpin'][2] = 
FluxTerms['SpinSpin'][3] = 
FluxTerms['SpinSpin'][4] = 
FluxTerms['SpinSpin'][5] = 
FluxTerms['SpinSpin'][6] = 
FluxTerms['SpinSpin'][7] = 
FluxTerms['SpinSpin'][8] = 

# <markdowncell>

# The spin-orbit terms in the flux are now known to 4.0pN.  These terms come from Eq. (4.9) of [Marsat et al. (2013)](http://arxiv.org/abs/1307.6793v1):

# <codecell>

FluxTerms['SpinOrbit'][0] = 0
FluxTerms['SpinOrbit'][1] = 0
FluxTerms['SpinOrbit'][2] = 0
FluxTerms['SpinOrbit'][3] = 4*S_l - frac(5,4)*delta*Sigma_l
FluxTerms['SpinOrbit'][4] = 0
FluxTerms['SpinOrbit'][5] = (-frac(9,2)+frac(272,9)*nu)*S_l + (-frac(13,16)+frac(43,4)*nu)*delta*Sigma_l
FluxTerms['SpinOrbit'][6] = (-16*pi)*S_l + (-frac(31,6)*pi)*delta*Sigma_l
FluxTerms['SpinOrbit'][7] = (frac(476645,6804)+frac(6172,189)*nu-frac(2810,27)*nu**2)*S_l \
                            + (frac(9535,336)+frac(1849,126)*nu-frac(1501,36)*nu**2)*delta*Sigma_l
FluxTerms['SpinOrbit'][8] = (-frac(3485,96)+frac(13879,72)*nu)*pi*S_l \
                            + (-frac(7163,672)+frac(130583,2016)*nu)*pi*delta*Sigma_l

# <markdowncell>

# The tidal-coupling terms come in to the energy at relative 5pN order, and are known to 6pN order.
# 
# These terms come from Eq. (3.6) of [Vines et al. (2011)](http://prd.aps.org/abstract/PRD/v83/i8/e084051).  Note their unusual convention for mass ratios, where $\chi_1 = m_1/m$ in their notation; in particular, $\chi$ is not a spin parameter.  Also note that $\hat{\lambda} = \lambda_2 v^{10}/(m_1+m_2)^5$, and we need to add the coupling terms again with $1 \leftrightarrow 2$.  Finally, note the normalization difference, where a different overall factor is used, leading to a sign difference.

# <codecell>

FluxTerms['NSTidal'][0] = 0
FluxTerms['NSTidal'][1] = 0
FluxTerms['NSTidal'][2] = 0
FluxTerms['NSTidal'][3] = 0
FluxTerms['NSTidal'][4] = 0
FluxTerms['NSTidal'][5] = 0
FluxTerms['NSTidal'][6] = 0
FluxTerms['NSTidal'][7] = 0
FluxTerms['NSTidal'][8] = 0
FluxTerms['NSTidal'][9] = 0
FluxTerms['NSTidal'][10] = (12-18/m2)*lambda2 + (12-18/m1)*lambda1
FluxTerms['NSTidal'][11] = 0
FluxTerms['NSTidal'][12] = (704+1803*m2-4501*m2**2+2170*m2**3)*lambda2/(28*m2) \
                           + (704+1803*m1-4501*m1**2+2170*m1**3)*lambda1/(28*m1)

# Note that the above terms should be divided by (m1+m2)**5, and m1 and m2 should be divided by (m1+m2),
# except that here we use units with m1+m2=1

# <codecell>

print FluxTerms

# <codecell>

FluxCoefficient*sum([FluxTerms['Nonspinning'][i]*v**i for i in range(8)])

# <headingcell level=1>

# Collected flux terms

# <codecell>

def Flux(SpinTerms=True, PrecessingSpinTerms=False, NSTidalTerms=False, ExtraEMRITerms=False) :
    """
    Return an expression for the GW flux with the given options.
    
    """
    FluxSum = sum([ for ])
    for
    # This is not yet implemented (the above should break evaluation of this cell)
    return FluxCoefficient*FluxSum

