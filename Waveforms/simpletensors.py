from __future__ import division
from sympy import simplify, var, Function
import sympy
from sympy import *
from sympy import Rational as frac

class _VectorFunctionOfTime(Function):
    """\
    This is just a base class for deriving other vectors from.

    You probably won't need to use this class directly; it is just the
    base class for the class created in VectorFunctionFactory below.

    This used to be more important because there were other subclasses
    of this, but now it's mostly just for convenient separation of a
    few methods, to clarify what's going on in the factory.

    """
    components = None
    def __or__(self,other):
        """
        In keeping with the notation of various other packages
        (most notably sympy.galgebra), contraction is denoted
        by the bitwise `or` operator.  So the contraction of
        vectors `v` and `w` is just `v|w`.  This notation will
        be used in the tensor classes also.
        """
        return sum( s*o for s,o in zip(self, other) )
    def __iter__(self):
        for c in self.__class__.components: yield c
    def diff(self, *args, **kwargs):
        return self._eval_derivative(*args, **kwargs)

def VectorFunctionFactory(Name, ComponentFunctions, DerivativeFunction=None):
    """Create a new vector function

    This function creates a class that is a subclass of
    `_VectorFunctionOfTime`, which is itself a subclass of
    `sympy.Function`.  Thus, the resulting quantity should be
    differentiable.  The class is created, and a class variable set to
    store the input components.  Then, if no derivative has been
    defined, one is defined for you, assuming that each component
    should be differentiated.  Alternatively, you can pass in a lambda
    function as the final argument.  It should take as inputs `self,
    *args, **kwargs`, and evaluate the value of the derivative from
    that information.  Finally, the class is renamed to carry the
    input name, so that sympy output looks nice, etc.
    """
    class Vector(_VectorFunctionOfTime):
        components = ComponentFunctions
    if(DerivativeFunction):
        Vector._eval_derivative = DerivativeFunction
    else:
        Vector._eval_derivative = lambda self, *args, **kwargs: \
                                  VectorFunctionFactory(Name+'Dot',
                                                        [c._eval_derivative(args[0])
                                                         for c in ComponentFunctions])(self.args[0])
    Vector.__name__ = Name
    return Vector



class TensorProduct(object):
    LaTeXProductString = r'\otimes'

    def __init__(self, *vectors, **kwargs):
        if(len(vectors)==1 and isinstance(vectors[0], TensorProduct)) :
            self.coefficient = vectors[0].coefficient
            self.vectors = vectors[0].vectors
        else:
            if('coefficient' in kwargs):
                self.coefficient = kwargs['coefficient']
            else:
                self.coefficient = 1
            self.vectors = list(vectors)

    def rank(self):
        return len(self.vectors)

    def __iter__(self):
        for v in self.vectors: yield v

    def has_same_basis_element(self, B):
        return self.vectors == B.vectors

    def __or__(self,B):
        if(B.rank() != self.rank()):
            raise ValueError("Cannot contract rank-{0} tensor with rank-{1} tensor.".format(self.rank(), B.rank()))
        from sympy import prod
        if(isinstance(B, TensorProduct)):
            return (self.coefficient*B.coefficient)*prod([v|w for v,w in zip(self, B)])
        else:
            try:
                if(B._is_tensor):
                    return sum( [(self|t_p) for t_p in B] )
            except AttributeError:
                raise ValueError("Don't know how to contract TensorProduct with '{0}'".format(type(B)))

    def trace(self, j, k):
        coefficient = sympy.simplify(self.coefficient * (self.vectors[j]|self.vectors[k]))
        if(self.rank()==2):
            return coefficient
        if(coefficient==0): return 0
        new = type(self)()
        new.vectors = list(v for i,v in enumerate(self) if (i!=j and i!=k))
        new.coefficient = coefficient
        return new

    def __mul__(self, B):
        """
        Return the scalar or tensor product
        """
        #print('TensorProduct.__mul__')
        #print(type(B), type(self))
        try:
            if(B._is_tensor):
                # Fall back to Tensor.__rmul__ by:
                return NotImplemented
        except AttributeError:
            pass
        if(isinstance(B, TensorProduct)):
            # Do tensor product
            new = type(self)(self)
            new.vectors = self.vectors + B.vectors
            new.coefficient = sympy.simplify( self.coefficient * B.coefficient )
            #print('TensorProduct.__mul__ return 1')
            return new
        elif(isinstance(B, _VectorFunctionOfTime)):
            new = type(self)(self)
            new.vectors = self.vectors + [B,]
            #print('TensorProduct.__mul__ return 2')
            return new
        else:
            # Otherwise, try scalar multiplication
            if(sympy.simplify(B)==0): return 0
            new = type(self)(self)
            new.coefficient = sympy.simplify( new.coefficient * B )
            #print('TensorProduct.__mul__ return 3')
            return new

    def __rmul__(self, B):
        """
        Return the scalar or tensor product
        """
        #print('TensorProduct.__rmul__')
        #print(type(B), type(self))
        try:
            if(B._is_tensor):
                # Fall back to Tensor.__mul__ by:
                return NotImplemented
        except AttributeError:
            pass
        if(isinstance(B, TensorProduct)):
            new = type(self)(self)
            new.vectors = B.vectors + self.vectors
            new.coefficient = sympy.simplify( self.coefficient * B.coefficient )
            #print('TensorProduct.__rmul__ return 1')
            return new
        elif(isinstance(B, _VectorFunctionOfTime)):
            new = type(self)(self)
            new.vectors = [B,] + self.vectors
            #print('TensorProduct.__rmul__ return 2')
            return new
        else:
            # Otherwise, try scalar multiplication
            if(sympy.simplify(B)==0): return 0
            new = type(self)(self)
            new.coefficient = sympy.simplify( new.coefficient * B )
            #print('TensorProduct.__rmul__ return 3')
            return new

    # These two will be defined below, once we have defined Tensor
    # objects:
    # TensorProduct.__add__ = lambda self, T: Tensor(self)+T
    # TensorProduct.__radd__ = TensorProduct.__add__

    def __str__(self):
        if(self.coefficient==1):
            return '[' + '*'.join([str(v) for v in self]) + ']'
        return '[('+str(self.coefficient)+')' + '*' + '*'.join([str(v) for v in self]) + ']'

    def __repr__(self):
        if(self.coefficient==1):
            return '[' + '*'.join([repr(v) for v in self]) + ']'
        return '[('+repr(self.coefficient)+')' + '*' + '*'.join([repr(v) for v in self]) + ']'

    def _latex_str_(self):
        if(self.coefficient==1):
            return r'\,\left[' + type(self).LaTeXProductString.join([sympy.latex(v) for v in self]) + r'\right]'
        return r'\left[\left('+sympy.latex(self.coefficient)+r'\right)' + r'\,' \
            + type(self).LaTeXProductString.join([sympy.latex(v) for v in self]) + r'\right]'

    def _repr_latex_(self):
        return '$'+ self._latex_str_() + '$'



class Tensor(object):
    def __init__(self, *tensor_products, **kwargs):
        if(len(tensor_products)==1 and isinstance(tensor_products[0], Tensor)) :
            self.tensor_products = tensor_products[0].tensor_products
        elif(len(tensor_products)==1 and isinstance(tensor_products[0], TensorProduct)) :
            self.tensor_products = [tensor_products[0],]
        else:
            self.tensor_products = list(tensor_products)
        if(len(self.tensor_products)>0):
            rank=self.tensor_products[0].rank()
            for tensor_product in self.tensor_products:
                if(tensor_product.rank() != rank):
                    raise ValueError("Cannot add rank-{0} tensor to rank-{1} tensors.".format(tensor_product.rank(), rank))

    @property
    def _is_tensor(self):
        "Since this is a property, it can be called without parentheses"
        return True

    def rank(self):
        if(len(self.tensor_products)==0):
            return 0
        return self.tensor_products[0].rank()

    def __iter__(self):
        for t_p in self.tensor_products: yield t_p

    def compress(self):
        #display(self)
        removed_elements = []
        for i in range(len(self.tensor_products)):
            if(i in removed_elements):
                continue
            for j in range(i+1,len(self.tensor_products)):
                if(j in removed_elements):
                    continue
                if self.tensor_products[i].has_same_basis_element(self.tensor_products[j]):
                    #print("Removing {0} because {1} is already here".format(j,i))
                    self.tensor_products[i].coefficient = \
                        sympy.simplify( self.tensor_products[i].coefficient + self.tensor_products[j].coefficient )
                    removed_elements.append(j)
            if removed_elements:
                if(self.tensor_products[i].coefficient==0):
                    removed_elements += [i]
                #print("Removing {0} because {1} is already here".format(removed_elements,i))
        self.tensor_products = list(t_p for i,t_p in enumerate(self) if i not in removed_elements)
        if not self.tensor_products:
            return 0
        return self

    def __add__(self, T):
        if(T==0):
            return self
        if(self.rank()==0):
            return T
        if(T.rank()==0):
            return self
        if(T.rank()!=self.rank()):
            raise ValueError("Cannot add rank-{0} tensor to rank-{1} tensor.".format(T.rank(), self.rank()))
        new = Tensor(self)
        if(isinstance(T, Tensor)) :
            new.tensor_products = self.tensor_products + T.tensor_products
        elif(isinstance(T, TensorProduct)) :
            new.tensor_products = self.tensor_products + [T,]
        return new

    def __radd__(self, T):
        return self+T

    def __or__(self, B):
        if(B.rank() != self.rank()):
            raise ValueError("Cannot contract rank-{0} tensor with rank-{1} tensor.".format(self.rank(), B.rank()))
        if(isinstance(B, Tensor)) :
            return sum([T1|T2 for T1 in self for T2 in B])
        elif(isinstance(B, TensorProduct)) :
            return sum([T1|B  for T1 in self])

    def trace(self, j=0, k=1):
        t = sum([T.trace(j,k) for T in self])
        try:
            #print("Trying to compress the Tensor trace")
            return t.compress()
        except:
            #print("Failed to compress the Tensor trace")
            return t

    def __mul__(self, B):
        #print('Tensor.__mul__')
        if(isinstance(B, Tensor)):
            T = Tensor()
            #print('Tensor.__mul__ starting list 1')
            T.tensor_products = list(t_pA*t_pB for t_pA in self for t_pB in B)
            #print('Tensor.__mul__ return 1')
            return T.compress()
        else:
            T = Tensor()
            #print('Tensor.__mul__ starting list 2')
            T.tensor_products = list(t_p*B for t_p in self)
            #print('Tensor.__mul__ return 2')
            return T

    def __rmul__(self, B):
        #print('Tensor.__rmul__')
        if(isinstance(B, Tensor)):
            T = Tensor()
            #print('Tensor.__rmul__ starting list 1')
            T.tensor_products = list(t_pB*t_pA for t_pA in self for t_pB in B)
            #print('Tensor.__rmul__ return 1')
            return T.compress()
        else:
            T = Tensor()
            #print('Tensor.__rmul__ starting list 2')
            T.tensor_products = list(B*t_p for t_p in self)
            #print('Tensor.__rmul__ return 2')
            return T

    def __str__(self):
        return '{'+'\n'.join([str(t_p) for t_p in self])+'}'

    def __repr__(self):
        return '{'+'\n'.join([repr(t_p) for t_p in self])+'}'

    def _latex_str_(self):
        return r'&\left\{' + r' \right. \nonumber \\&\quad \left. + '.join([t_p._latex_str_() for t_p in self]) + r'\right\}'

    def _repr_latex_(self):
        return r'\begin{align}'+ self._latex_str_() + r'\end{align}'


# Since the sum of two `TensorProduct`s is a `Tensor`, we
# had to wait until we got here to define these methods:
TensorProduct.__add__ = lambda self, T: Tensor(self)+T
TensorProduct.__radd__ = TensorProduct.__add__


# In[10]:

class SymmetricTensorProduct(TensorProduct):
    """
    Specialized class for symmetric tensor products

    **Note:**  If you multiply a SymmetricTensorProduct by
    any TensorProduct (even if it's not symmetric), the result
    will be symmetric.  This makes it easy to make STPs, but is
    not how real tensor products work.

    This is a subclass of `TensorProduct` with the necessary
    methods overridden.  Because it is subclassed, and `Tensor`
    isn't very invasive, we can easily create tensors by adding
    symmetric tensor products, and the `Tensor` need not even
    know that it is symmetric.
    """
    LaTeXProductString = r'\otimes_{\mathrm{s}}'

    def __init__(self, *vectors, **kwargs):
        TensorProduct.__init__(self, *vectors, **kwargs)

    def __call__(self, *indices):
        """
        Get the (symmetrized) element at this index.
        """
        if(len(indices)!=self.rank()):
            print(indices)
            raise ValueError("You need {0} indices to index a rank-{0} tensor, not {1}.".format(self.rank(), len(indices)))
        from itertools import permutations
        from math import factorial
        from sympy import prod
        symmetrized_element = 0
        for index in permutations(indices, len(indices)):
            symmetrized_element += prod([v[i] for v,i in zip(self, index)])
            #print(index, element, symmetrized_element)
        return symmetrized_element*self.coefficient*Rational(1,factorial(self.rank()))

    def has_same_basis_element(self, B):
        from collections import Counter
        return Counter(self.vectors) == Counter(B.vectors)

    def ordered_as(self, index_set):
        for i in index_set:
            yield self.vectors[i]

    def __or__(self,B):
        if(B.rank() != self.rank()):
            raise ValueError("Cannot contract rank-{0} tensor with rank-{1} tensor.".format(self.rank(), B.rank()))
        from itertools import permutations
        from sympy import prod
        if(isinstance(B, TensorProduct)):
            # If B is actually a SymmetricTensorProduct, it suffices to just
            # iterate over rearrangements of `self`.
            #return self.coefficient*B.coefficient*prod([sum([v[i]*w[i] for i in range(3)]) for v,w in zip(self, B)])
            coefficient = sympy.simplify(self.coefficient*B.coefficient*frac(1,factorial(self.rank())))
            if(coefficient==0): return 0
            return sympy.simplify( coefficient * sum([prod([v|w for v,w in zip(self.ordered_as(index_set), B)])
                                                      for index_set in permutations(range(self.rank()))]) )
        elif(isinstance(B, Tensor)):
            return sum( [self|t_p for t_p in B] )
        else:
            raise ValueError("Don't know how to contract SymmetricTensorProduct with '{0}'".format(type(B)))

    def trace(self, j=0, k=1):
        """
        Any input elements are ignored, since we will be symmetrizing anyway
        """
        T = Tensor()
        from itertools import permutations
        for j,k in permutations(range(self.rank()), 2):
            coefficient = sympy.simplify( self.coefficient * (self.vectors[j]|self.vectors[k]) )
            if(self.rank()==2):
                return coefficient
            if(coefficient==0):
                continue
            new = SymmetricTensorProduct()
            new.vectors = list(v for i,v in enumerate(self.vectors) if (i!=j and i!=k))
            new.coefficient = coefficient*frac(1,factorial(self.rank()))
            T = T+new
        return T.compress()
