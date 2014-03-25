from __future__ import division
# import sympy
from sympy import simplify, Function, latex, prod, Symbol, flatten, factorial, Derivative
from sympy import Rational as frac
from copy import deepcopy

def DifferentiateString(s, param='t'):
    "Add latex to string indicating differentiation by time."
    if(s.startswith(r'\partial_{0}^{{'.format(param))):
        from re import sub
        def increment_match(match):
            return r'\partial_{0}^{{{1}}}{2}'.format(param, int(match.group(1))+1, match.group(2))
        return sub(r'\\partial_.*?\^\{(.*)\}(.*)', increment_match, s)
    elif(s.startswith(r'\partial_{0}'.format(param))):
        return r'\partial_{0}^{{2}}'.format(param) + s[10:]
    else:
        return r'\partial_{0} '.format(param) + s

def DelimitString(S, latex=True):
    "Surround string by parentheses, brackets, or braces as appropriate"
    if(latex):
        left, right = [r'\left', r'\right']
        DelimiterOpeners, DelimiterClosers = [['(','[','\{'], [')',']','\}']]
    else:
        left, right = ['', '']
        DelimiterOpeners, DelimiterClosers = ['([{', ')]}']
    def FindFirst(S, D):
        for s in S:
            for d in D:
                if(s==d):
                    return d
    FirstDelim = FindFirst(S, DelimiterOpeners)
    if(not FirstDelim):
        return left+'('+S+right+')'
    NewDelimiterIndex = (DelimiterOpeners.index(FindFirst(S, DelimiterOpeners))+1) % len(DelimiterOpeners)
    return r'{0}{1} {2} {3}{4}'.format(left, DelimiterOpeners[NewDelimiterIndex],
                                       S, right, DelimiterClosers[NewDelimiterIndex])


####################################
### First, a few vector thingies ###
####################################
class _VectorFunction(Function):
    """\
    This is just a base class for deriving other vectors from.

    You probably won't need to use this class directly; it is just the
    base class for the class created in VectorFunction factory below.

    This used to be more important because there were other subclasses
    of this, but now it's mostly just for convenient separation of a
    few methods, to clarify what's going on in the factory.

    """
    components = None
    @property
    def name(self):
        return self.__class__.__name__
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
    def __mul__(self, other):
        if(hasattr(other, '_is_tensor_product') or hasattr(other, '_is_tensor')):
            return NotImplemented
    def __rmul__(self, other):
        if(hasattr(other, '_is_tensor_product') or hasattr(other, '_is_tensor')):
            return NotImplemented
    def _eval_derivative(self, s):
        """Return derivative of the function with respect to `s`.

        Note that we must take the chain rule into account in this
        function, but not in `fdiff`.

        """
        return self.fdiff(1) * self.args[0].diff(s)
    def fdiff(self, argindex=1):
        """Returns the first derivative of the function.

        To be precise, this is the first derivative with respect to
        the function's argument, whatever that argument may be.  And
        in particular, we don't need to account for the chain rule.

        """
        from sympy.core.function import ArgumentIndexError
        if (argindex != 1):
            raise ArgumentIndexError(self, argindex)
        return VectorFunction(DifferentiateString(self.name, self.args[0]),
                              [Derivative(c, self.args[0], evaluate=False)
                               for c in self])(self.args[0])


def VectorFunction(Name, ComponentFunctions, DerivativeFunction=None):
    """Create a new vector function

    This function creates a class that is a subclass of
    `_VectorFunction`, which is itself a subclass of
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
    class Vector(_VectorFunction):
        components = ComponentFunctions
    if(DerivativeFunction):
        Vector.fdiff = DerivativeFunction
    Vector.__name__ = Name
    return Vector

def VectorConstant(Name, Components):
    return VectorFunction(Name, Components, lambda self, argindex=1: 0)



################################
### Now, the tensor products ###
################################
class TensorProductFunction(Function):
    @property
    def LaTeXProductString(self):
        if(self.symmetric):
            return r' \otimes_{\mathrm{s}} '
        else:
            return r' \otimes '

    @property
    def _is_tensor_product(self):
         return True

    @property
    def rank(self):
        return len(self.vectors)

    def __iter__(self):
        for v in self.vectors: yield v

    def has_same_basis_element(self, B):
        if(self.symmetric):
            from collections import Counter
            return Counter(self.vectors) == Counter(B.vectors)
        return self.vectors == B.vectors

    def ordered_as(self, index_set):
        for i in index_set:
            yield self.vectors[i]

    def __or__(self,B):
        if(B.rank != self.rank):
            raise ValueError("Cannot contract rank-{0} tensor with rank-{1} tensor.".format(self.rank, B.rank))
        if(isinstance(B, TensorProductFunction)):
            if(self.symmetric):
                from itertools import permutations
                # It suffices to just iterate over rearrangements of `self`.
                coefficient = simplify(self.coefficient*B.coefficient*frac(1,factorial(self.rank)))
                if(coefficient==0): return 0
                return simplify( coefficient * sum([prod([v|w for v,w in zip(self.ordered_as(index_set), B)])
                                                    for index_set in permutations(range(self.rank))]) )
            return (self.coefficient*B.coefficient)*prod([v|w for v,w in zip(self, B)])
        else:
            try:
                return sum( [(self|t_p) for t_p in B] )
            except AttributeError:
                raise ValueError("Don't know how to contract TensorProductFunction with '{0}'".format(type(B)))

    def trace(self, j=-1, k=-1):
        # print("\nTP.trace({0}, {1}) running on rank {2} {3}tensor {4}".format(j,k,self.rank, ('symmetric ' if self.symmetric else ''), self))
        if(not self.symmetric and (j==-1 or k==-1)):
            raise TypeError("trace() takes exactly 3 arguments for non-symmetric tensor products (1 given)")
        coefficient = simplify(self.coefficient * (self.vectors[j]|self.vectors[k]))
        if(self.rank==2):
            # print("Finished TP.trace 1\n")
            return coefficient
        if(self.symmetric):
            from itertools import permutations
            T = 0
            for j,k in permutations(range(self.rank), 2):
                # print(j,k)
                coefficient = simplify( self.coefficient * (self.vectors[j]|self.vectors[k]) )
                if(coefficient==0):
                    continue
                T += TensorProduct(list(v for i,v in enumerate(self.vectors) if (i!=j and i!=k)),
                                   coefficient = coefficient*frac(1,factorial(self.rank)),
                                   symmetric=True)
            # print("Finished TP.trace 2\n")
            try:
                return T.compress()
            except AttributeError:
                return T
        elif(coefficient==0):
            # print("Finished TP.trace 3\n")
            return 0
        # print("Finished TP.trace 4\n")
        return TensorProduct(list(v for i,v in enumerate(self) if (i!=j and i!=k)),
                             coefficient=coefficient, symmetric=False)

    def __mul__(self, B):
        """
        Return the scalar or tensor product
        """
        # print('TensorProductFunction.__mul__')
        # print(type(B), type(self))
        if(simplify(B)==0): return 0
        if(hasattr(B, '_is_tensor') and B._is_tensor):
            # Fall back to Tensor.__rmul__ by doing this:
            return NotImplemented
        if(isinstance(B, TensorProductFunction)):
            # Do tensor product
            # print('TensorProductFunction.__mul__ return 1')
            return TensorProduct(self.vectors + B.vectors,
                                 coefficient=simplify( self.coefficient * B.coefficient ),
                                 symmetric = self.symmetric)
        elif(isinstance(B, _VectorFunction)):
            # print('TensorProductFunction.__mul__ return 2')
            return TensorProduct(self.vectors + [B,],
                                 coefficient=self.coefficient,
                                 symmetric = self.symmetric)
        else:
            # Otherwise, try scalar multiplication
            # print('TensorProductFunction.__mul__ return 3')
            return TensorProduct(self.vectors,
                                 coefficient=self.coefficient*B,
                                 symmetric = self.symmetric)

    def __rmul__(self, B):
        """
        Return the scalar or tensor product
        """
        # print('TensorProductFunction.__rmul__')
        # print(type(B), type(self))
        if(simplify(B)==0): return 0
        if(hasattr(B, '_is_tensor') and B._is_tensor):
            # Fall back to Tensor.__rmul__ by doing this:
            return NotImplemented
        if(isinstance(B, TensorProductFunction)):
            # Do tensor product
            # print('TensorProductFunction.__rmul__ return 1')
            return TensorProduct(B.vectors+self.vectors,
                                 coefficient=simplify( B.coefficient * self.coefficient ),
                                 symmetric = self.symmetric)
        elif(isinstance(B, _VectorFunction)):
            # print('TensorProductFunction.__rmul__ return 2')
            return TensorProduct([B,] + self.vectors,
                                 coefficient=self.coefficient,
                                 symmetric = self.symmetric)
        else:
            # Otherwise, try scalar multiplication
            # print('TensorProductFunction.__rmul__ return 3')
            return TensorProduct(self.vectors,
                                 coefficient=B*self.coefficient,
                                 symmetric = self.symmetric)

    # These two will be defined below, once we have defined Tensor
    # objects:
    # TensorProductFunction.__add__ = lambda self, T: Tensor(self)+T
    # TensorProductFunction.__radd__ = TensorProductFunction.__add__

    def _eval_derivative(self, s):
        """Return derivative of the function with respect to `s`.

        Note that we must take the chain rule into account in this
        function, but not in `fdiff`.

        """
        return self.fdiff(1) * self.args[0].diff(s)

    def fdiff(self, argindex=1):
        """Returns the first derivative of the function.

        To be precise, this is the first derivative with respect to
        the function's argument, whatever that argument may be.  And
        in particular, we don't need to account for the chain rule.

        """
        from sympy.core.function import ArgumentIndexError
        if (argindex != 1):
            raise ArgumentIndexError(self, argindex)
        return VectorFunction(DifferentiateString(self.name, self.args[0]),
                              [Derivative(c, self.args[0], evaluate=False)
                               for c in self])(self.args[0])

    def fdiff(self, argindex=1):
        return TensorProduct(list(self.vectors),
                             coefficient = Derivative(self.coefficient, self.args[0], evaluate=False),
                             symmetric = self.symmetric) \
            + sum( [ TensorProduct([t.fdiff(argindex) if i==j else t for j,t in enumerate(self)],
                                   coefficient = self.coefficient, symmetric=self.symmetric)
                     for i in range(self.rank) ] )

    def __str__(self):
        if(self.coefficient==1):
            return DelimitString('*'.join([str(v) for v in self]), latex=False)
        return DelimitString( DelimitString(str(self.coefficient))
                              + '*' + '*'.join([str(v) for v in self]) )

    def __repr__(self):
        if(self.coefficient==1):
            return DelimitString('*'.join([repr(v) for v in self]), latex=False)
        return DelimitString( DelimitString(str(self.coefficient))
                              + '*' + '*'.join([repr(v) for v in self]) )

    def _latex_str_(self):
        if(self.coefficient==1):
            return DelimitString( self.LaTeXProductString.join([latex(v) for v in self]) )
        return DelimitString( DelimitString(latex(self.coefficient)) + '\, '
                              + self.LaTeXProductString.join([latex(v) for v in self]) )

    def _latex(self, printer):
        "Sympy looks for this when latex printing is on"
        # printer._settings['mode'] = 'equation*'
        return self._latex_str_()

    def _repr_latex_(self):
        return '$'+ self._latex_str_() + '$'

_TensorProduct_count = 0
def TensorProduct(*input_vectors, **kwargs):
    if('coefficient' in kwargs and kwargs['coefficient']==0):
        return 0

    ## First, go through and make the data nice
    if(len(input_vectors)==0):
        # Since TensorProducts are multiplicative, the empty object should be 1
        return kwargs.get('coefficient', 1)
    if(len(input_vectors)==1 and isinstance(input_vectors[0], TensorProductFunction)) :
        vectors = list(input_vectors[0].vectors)
        coefficient = deepcopy(input_vectors[0].coefficient)
        symmetric = bool(input_vectors[0].symmetric)
    else:
        if(len(input_vectors)==1 and isinstance(input_vectors[0], list)):
            input_vectors = input_vectors[0]
        vectors = list(input_vectors)
        coefficient = deepcopy(kwargs.get('coefficient', 1))
        symmetric = bool(kwargs.get('symmetric', False))

    ## Now, create the object and set its data.  Because of sympy's
    ## caching, tensor products with different data need to be created
    ## as classes with different names.  So we just create a
    ## lighweight subclass with a unique name (the number at the end
    ## gets incremented every time we construct a tensor product).
    global _TensorProduct_count
    ThisTensorProductFunction = type('TensorProductFunction_'+str(_TensorProduct_count),
                                     (TensorProductFunction,), {})
    _TensorProduct_count += 1
    TP = ThisTensorProductFunction( *tuple( set( flatten( [v.args for v in vectors] ) ) ) )
    TP.vectors = vectors
    TP.coefficient = coefficient
    TP.symmetric = symmetric
    return TP

def SymmetricTensorProduct(*input_vectors, **kwargs):
    kwargs['symmetric'] = True
    return TensorProduct(*input_vectors, **kwargs)


##############################
### And tensors themselves ###
##############################
class TensorFunction(Function):
    """
    This is just a base class for deriving other tensors from.

    You probably won't need to use this class directly; it is just the
    base class for the class created in Tensor below.
    """

    @property
    def _is_tensor(self):
        """Since this is a property, it can be called without parentheses.
        But I want to make sure it applies to the instance, not the
        class, so it shouldn't just be a variable."""
        return True

    @property
    def rank(self):
        if(len(self.tensor_products)==0):
            return 0
        return self.tensor_products[0].rank

    def __iter__(self):
        if self.tensor_products:
            for t_p in self.tensor_products: yield t_p
        else:
            raise StopIteration()

    def compress(self):
        # print("\ncompressing")
        # print(self)
        removed_elements = []
        for i in range(len(self.tensor_products)):
            if(i in removed_elements):
                continue
            for j in range(i+1,len(self.tensor_products)):
                if(j in removed_elements):
                    continue
                if self.tensor_products[i].has_same_basis_element(self.tensor_products[j]):
                    # print("Removing {0} because {1} is already here".format(j,i))
                    self.tensor_products[i] = TensorProduct(self.tensor_products[i].vectors,
                                                            coefficient=simplify( self.tensor_products[i].coefficient + self.tensor_products[j].coefficient ),
                                                            symmetric=self.tensor_products[i].symmetric)
                    removed_elements.append(j)
            if removed_elements:
                if(self.tensor_products[i].coefficient==0):
                    removed_elements += [i]
                # print("Removing {0} because {1} is already here".format(removed_elements,i))
        self.tensor_products = list(t_p for i,t_p in enumerate(self) if i not in removed_elements)
        if not self.tensor_products:
            # print("compressed to 0")
            return 0
        # print("compressed to {0}\n".format(self))
        return self

    def __add__(self, T):
        if(T==0):
            return self
        if(self.rank==0):
            return T
        if(T.rank==0):
            return self
        if(T.rank!=self.rank):
            raise ValueError("Cannot add rank-{0} tensor to rank-{1} tensor.".format(T.rank, self.rank))
        if(isinstance(T, TensorFunction)) :
            return Tensor(self.tensor_products + T.tensor_products)
        elif(isinstance(T, TensorProductFunction)) :
            return Tensor(self.tensor_products + [T,])

    def __radd__(self, T):
        """Addition is commutative, but python might get here when T is a
        TensorProduct or something"""
        return self+T

    def __or__(self, B):
        if(B.rank != self.rank):
            raise ValueError("Cannot contract rank-{0} tensor with rank-{1} tensor.".format(self.rank, B.rank))
        if(isinstance(B, TensorFunction)) :
            return sum([T1|T2 for T1 in self for T2 in B])
        elif(isinstance(B, TensorProductFunction)) :
            return sum([T1|B  for T1 in self])

    def trace(self, j=0, k=1):
        # print("\nSumming trace")
        # print(self)
        # for T in self:
        #     # print(T)
        #     # print("")
        # print("The sum will run next:")
        t = sum([T.trace(j,k) for T in self])
        # print("Sum:", t)
        try:
            # print("Trying to compress the Tensor trace")
            return t.compress()
        except AttributeError:
            # print("Failed to compress the Tensor trace")
            return t

    def __mul__(self, B):
        # print('Tensor.__mul__')
        if(isinstance(B, TensorFunction)):
            # print('TensorFunction.__mul__ return 1')
            return Tensor(list(t_pA*t_pB for t_pA in self for t_pB in B))
        else:
            # print('TensorFunction.__mul__ return 2')
            return Tensor(list(t_p*B for t_p in self))

    def __rmul__(self, B):
        # print('TensorFunction.__rmul__')
        if(isinstance(B, TensorFunction)):
            # print('TensorFunction.__rmul__ return 1')
            return Tensor(list(t_pB*t_pA for t_pA in self for t_pB in B))
        else:
            # print('TensorFunction.__rmul__ return 2')
            return Tensor(list(B*t_p for t_p in self))

    def fdiff(self, argindex=1):
        return Tensor(list(t_p.fdiff(argindex) for t_p in self))

    def __str__(self):
        return DelimitString( '\n'.join([str(t_p) for t_p in self]), latex=False)

    def __repr__(self):
        return DelimitString( '\n'.join([repr(t_p) for t_p in self]), latex=False)

    def _latex_str_(self):
        return '&' + DelimitString( r' \right. \nonumber \\&\quad \left. + '.join(
            [t_p._latex_str_() for t_p in self]) )

    def _latex(self, printer):
        printer._settings['mode'] = 'align*'
        return self._latex_str_()
            # + DelimitString( ','.join([ str(printer._print(arg)) for arg in self.args ]) )

    def _repr_latex_(self):
        return r'\begin{align}'+ self._latex_str_() + r'\end{align}'

_Tensor_count = 0
def Tensor(*tensor_products):
    """Create a new tensor"""
    ## First, go through and make the data nice
    if(len(tensor_products)==0):
        # Since Tensor objects are additive, the empty object should be zero
        return 0
    if(len(tensor_products)==1 and isinstance(tensor_products[0], TensorFunction)) :
        tensor_products = list(tensor_products[0].tensor_products)
    elif(len(tensor_products)==1 and isinstance(tensor_products[0], TensorProductFunction)) :
        tensor_products = [tensor_products[0],]
    else:
        if(len(tensor_products)==1 and isinstance(tensor_products[0], list)):
            tensor_products = tensor_products[0]
        tensor_products = flatten(list(tensor_products))
        if(len(tensor_products)>0):
            rank=tensor_products[0].rank
            for t_p in tensor_products:
                if(t_p.rank != rank):
                    raise ValueError("Cannot add rank-{0} tensor to rank-{1} tensors.".format(t_p.rank, rank))

    ## Now, create the object and set its data.  Because of sympy's
    ## caching, tensors with different data need to be created as
    ## classes with different names.  So we just create a lighweight
    ## subclass with a unique name (the number at the end gets
    ## incremented every time we construct a tensor).
    global _Tensor_count
    ThisTensorFunction = type('TensorFunction_'+str(_Tensor_count),
                             (TensorFunction,), {})
    _Tensor_count += 1
    T = ThisTensorFunction( *tuple( set( flatten( [t_p.args for t_p in tensor_products] ) ) ) )
    T.tensor_products = tensor_products
    return T


# Since the sum of two `TensorProduct`s is a `Tensor`, we
# had to wait until we got here to define these methods:
TensorProductFunction.__add__ = lambda self, T: Tensor(self)+T
TensorProductFunction.__radd__ = TensorProductFunction.__add__
