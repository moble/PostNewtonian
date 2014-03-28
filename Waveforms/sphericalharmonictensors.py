"""\
Functions for generating spherical-harmonic tensors, etc.

[Thorne (1980)](http://link.aps.org/doi/10.1103/RevModPhys.52.299)
gives a nice review, along with the following formula for
$\mathcal{Y}^{\ell,m}_L$.  Also note the Eqs. (74) and the footnote on
page 32 of [Blanchet's Living Review
(2013)](http://arxiv.org/abs/1310.1528), which explains the necessary
normalizations for getting the metric perturbation modes from the
tensors.

"""

from __future__ import division
import sympy
from sympy import *
from sympy import Rational as frac
import simpletensors
from simpletensors import Vector, xHat, yHat, zHat
from simpletensors import TensorProduct, SymmetricTensorProduct, Tensor
init_printing()

var('vartheta, varphi')

DefaultOrthogonalRightHandedBasis=[xHat(t), yHat(t), zHat(t)]

def C(ell,m):
    return (-1)**abs(m) * sympy.sqrt( frac(2*ell+1,4) * frac(factorial(ell-m), factorial(ell+m)) / sympy.pi )

def a(ell,m,j):
    return frac((-1)**j, 2**ell * factorial(j) * factorial(ell-j)) * frac(factorial(2*ell-2*j), factorial(ell-m-2*j))

def YlmTensor(ell, m, OrthogonalRightHandedBasis=DefaultOrthogonalRightHandedBasis):
    if(ell<0 or abs(m)>ell):
        raise ValueError("YlmTensor({0},{1}) is undefined.".format(ell,m))
    from sympy import prod
    xHat, yHat, zHat = OrthogonalRightHandedBasis
    if(m<0):
        mVec = Tensor(SymmetricTensorProduct(xHat), SymmetricTensorProduct(yHat,coefficient=-sympy.I))
        #mVec = VectorFactory('mBarVec', [1,-sympy.I,0])(t)
    else:
        mVec = Tensor(SymmetricTensorProduct(xHat), SymmetricTensorProduct(yHat,coefficient=sympy.I))
        #mVec = VectorFactory('mVec', [1,sympy.I,0])(t)
    def TensorPart(ell,m,j):
        return sympy.prod((mVec,)*m) * SymmetricTensorProduct(*((zHat,)*(ell-2*j-m))) \
            * sympy.prod([sum([SymmetricTensorProduct(vHat, vHat) for vHat in OrthogonalRightHandedBasis]) for i in range(j)])
    if(m<0):
        Y = sum([TensorPart(ell,-m,j) * (C(ell,-m) * a(ell,-m,j))
                 for j in range(floor(frac(ell+m,2))+1) ]) * (-1)**(-m)
    else:
        Y = sum([TensorPart(ell,m,j) * (C(ell,m) * a(ell,m,j))
                 for j in range(floor(frac(ell-m,2))+1) ])
    try:
        Y.compress()
    except AttributeError:
        pass
    return Y

def YlmTensorConjugate(ell, m, OrthogonalRightHandedBasis=DefaultOrthogonalRightHandedBasis):
    return YlmTensor(ell, -m, OrthogonalRightHandedBasis) * (-1)**abs(m)

# This is Blanchet's version of the above
def alphalmTensor(ell, m, OrthogonalRightHandedBasis=DefaultOrthogonalRightHandedBasis):
    return YlmTensor(ell, -m, OrthogonalRightHandedBasis) * ( (-1)**abs(m) * (4*pi*factorial(ell)) / factorial2(2*ell+1) )

NVec = Vecto

r('NVec', r'\hat{N}', [sympy.sin(vartheta)*sympy.cos(varphi),
                                   sympy.sin(vartheta)*sympy.sin(varphi),
                                   sympy.cos(vartheta)])(t)

def NTensor(ell):
    return SymmetricTensorProduct(*((NVec,)*ell))

# These give the SWSH modes components from the given radiative tensors
def Ulm(U_L, m):
    ell = U_L.rank
    return (alphalmTensor(ell, m) | U_L) * (frac(4,factorial(ell))*sqrt(frac((ell+1)*(ell+2), 2*ell*(ell-1))))
def Vlm(V_L, m):
    ell = V_L.rank
    return (alphalmTensor(ell, m) | V_L) * (frac(-8,factorial(ell))*sqrt(frac(ell*(ell+2), 2*(ell+1)*(ell-1))))
def hlm(U_L, V_L, m):
    return ( -Ulm(U_L,m) + sympy.I * Vlm(V_L,m) ) / sympy.sqrt(2)
