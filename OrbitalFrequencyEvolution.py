# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# Always run this first
# NOTE: Do not define new variables in this notebook or in Variables.py; 
#       define them in Variables.ipynb.  Variables.py will be automatically
#       recreated when Variables.ipynb is saved.

from __future__ import division # This needs to be here, even though it's in Variables.py
execfile('ExecNotebook.ipy')
execnotebook('Variables.ipynb')
execnotebook('CodeOutput.ipynb')

# <headingcell level=1>

# TaylorT1 and TaylorT4

# <markdowncell>

# These two very similar approximants are the simplest in construction, and most widely applicable.  In particular, they can both be applied to precessing systems.  Each gives rise to the same system of ODEs that need to be integrated in time, except that the right-hand side for $dv/dt$ is expanded as a series in $v$ and truncated for TaylorT4.

# <codecell>

execnotebook('Flux.ipynb')
execnotebook('OrbitalEnergy.ipynb')
execnotebook('EnergyAbsorption.ipynb')
execnotebook('Precession.ipynb')

# Evaluate the flux and energy
Flux = FluxSum(IncompleteNonspinningTerms=True).subs(log(v), logv).subs(Pow(nu,3), nu__3).subs(Pow(nu,2), nu__2)
Energy = EnergySum(IncompleteNonspinningTerms=True).subs(log(v), logv).subs(Pow(nu,3), nu__3).subs(Pow(nu,2), nu__2)

# Differentiate Energy WRT v, treating log terms as constant
dEnergydv = horner(diff(Energy, v))

# Evaluate the energy absorption by the BHs, and make substitutions so that the Horner form is nice
Absorption = AbsorptionCoefficient*AlviTerm
Absorption = Absorption.subs(BasicSubstitutions[sqrt1Mchi1chi1], sqrt1Mchi1chi1)
Absorption = Absorption.subs(BasicSubstitutions[sqrt1Mchi2chi2], sqrt1Mchi2chi2)
Absorption = Absorption.subs(log(v), logv).subs(Pow(nu,3), nu__3)
Absorption = Absorption.subs(Pow(nu,2), nu__2)

# <codecell>

PNOrder = frac(7,2)
dvdt = series( series(- (Flux + Absorption), x=v, x0=0, n=10+2*PNOrder).removeO()
               / series(dEnergydv, x=v, x0=0, n=1+2*PNOrder).removeO(),
               x=v, x0=0, n=9+2*PNOrder).removeO()

# <codecell>

print dvdt

# <codecell>

print("""  
  int TaylorT1RHS(double t, const double* y, double* dydt) {
    // Stop integrating if v is greater than or equal to 1.0
    if(y[0]>=1.0) { return GSL_ETOLX; }
    
    // Fundamental variables (except for the Quaternions)
    const double v = y[0];
    const double Phi = y[1];
    const double chi1_x = y[2];
    const double chi1_y = y[3];
    const double chi1_z = y[4];
    const double chi2_x = y[5];
    const double chi2_y = y[6];
    const double chi2_z = y[7];
    const double Lhat_N_x = y[8];
    const double Lhat_N_y = y[9];
    const double Lhat_N_z = y[10];
    const double gamma = y[11];
    const GWFrames::Quaternion Rax = 
        GWFrames::sqrtOfRotor(-GWFrames::normalized(GWFrames::Quaternion(0., Lhat_N_x, Lhat_N_y, Lhat_N_z))*zHat);
    const GWFrames::Quaternion R = Rax * GWFrames::exp(((gamma+Phi)/2.)*zHat);
    const GWFrames::Quaternion nhatQ = R*xHat*R.conjugate();
    const double nhat_x = nhatQ[1];
    const double nhat_y = nhatQ[2];
    const double nhat_z = nhatQ[3];""")
print(CCodeOutput(['dEnergydv',
                   'Flux',
                   'Absorption',
                   'gamma_Lhat_N',
                   'a_ell_overvcubed',
                   'const double Omega_Lhat_N = gamma_Lhat_N*a_ell_overvcubed;',
                   'Omega1',
                   'Omega2']))
print("""    
    // Calculate derivatives of the fundamental variables
    dydt[0] = - (Flux + Absorption) / dEnergydv;
    if(dydt[0]<0.0) { return GSL_ETOLF; } // Stop integrating if v is decreasing
    dydt[1] = v*v*v;
    double Omega_spin1[3] = {Omega1*Lhat_N_x, Omega1*Lhat_N_y, Omega1*Lhat_N_z};
    double Omega_spin2[3] = {Omega2*Lhat_N_x, Omega2*Lhat_N_y, Omega2*Lhat_N_z};
    double Omega_prec[3]  = {Omega_Lhat_N*nhat_x, Omega_Lhat_N*nhat_y, Omega_Lhat_N*nhat_z};
    cross(&Omega_spin1[0], &y[2], &dydt[2]);
    cross(&Omega_spin2[0], &y[5], &dydt[5]);
    cross(&Omega_prec[0],  &y[8], &dydt[8]);
    const GWFrames::Quaternion adot(0., dydt[8], dydt[9], dydt[10]);
    const GWFrames::Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[10]))*adot*zHat - (dydt[10]/(2+2*y[10]))*Rax );
    const GWFrames::Quaternion dgammadt = 2*(Rax.conjugate() * Raxdot * zHat);
    dydt[11] = dgammadt[0];
    
    // GSL expects this to be the returned quantity if everything went well
    return GSL_SUCCESS;
  }
""")

# <headingcell level=1>

# TaylorT2 and TaylorT3

# <markdowncell>

# These two approximants are also closely related to each other.

# <codecell>


