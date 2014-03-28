from __future__ import division
# import sympy
from sympy import simplify, Function, latex, prod, Symbol, Mul, Add, flatten, factorial, Derivative, diff, sympify
from sympy import Rational as frac
from copy import deepcopy
from string import maketrans

TranslationTable = maketrans(r'[]{}',r'()()')

def DifferentiateString(s, param='t'):
    "Add latex to string indicating differentiation by time."
    if(s.startswith(r'\partial_{0}^{{'.format(param))):
        from re import sub
        def increment_match(match):
            return r'\partial_{0}^{{{1}}}{2}'.format(param, int(match.group(1))+1, match.group(2))
        return sub(r'\\partial_.*?\^\{(.*?)\}(.*)', increment_match, s)
    elif(s.startswith(r'\partial_{0}'.format(param))):
        return r'\partial_{0}^{{2}}'.format(param) + s[10:]
    else:
        return r'\partial_{0} '.format(param) + s

def DelimitString(S, latex=True):
    "Surround string by parentheses, brackets, or braces as appropriate"
    if(latex):
        left, right = [r'\left', r'\right']
        DelimiterOpeners, DelimiterClosers = [['(','[',r'\{'], [')',']',r'\}']]
    else:
        left, right = ['', '']
        DelimiterOpeners, DelimiterClosers = ['([{', ')]}']
    def FindFirst(S, D):
        for i,s in enumerate(S):
            if(s=='\\' and i+1<len(S)):
                s += S[i+1]
            for d in D:
                if(s==d):
                    return d
    FirstDelim = FindFirst(S, DelimiterOpeners)
    if(not FirstDelim):
        return left+'('+S+right+')'
    NewDelimiterIndex = (DelimiterOpeners.index(FindFirst(S, DelimiterOpeners))+1) % len(DelimiterOpeners)
    return r'{0}{1} {2} {3}{4}'.format(left, DelimiterOpeners[NewDelimiterIndex],
                                       S, right, DelimiterClosers[NewDelimiterIndex])

def LatexSubs(S, subsargs, subskwargs):
    """Construct latex showing that variables have been substituted

    Note that `subs` accepts either one or two arguments.  If two, the
    first is the variable in the expression currently, and the second
    is the replacement value.  If there is only one argument, it may
    be a dictionary of replacement pairs or an iterable of tuples.
    There is some clever logic in the `subs` function itself for
    dealing with a dictionary involving operation counts, etc.; I'll
    will be stupid about it and just iterate through the dictionary.

    """
    if(len(subsargs)>2 or len(subsargs)<1):
        raise ValueError('`subs` got {0} args: "{1}"'.format(len(subsargs), subsargs))
    if(len(subsargs)==2):
        return r'\left.' + S + r'\right|_{{{0}={1}}}'.format(latex(subsargs[0]), latex(subsargs[1]))
    subsargs = subsargs[0]
    subsstring = ','.join('{0}={1}'.format(latex(a),latex(b))
                          for a,b in (subsargs.items() if isinstance(subsargs,dict) else subsargs))
    return r'\left.' + S + r'\right|_{{{0}}}'.format(subsstring)


def ReduceExpr(expr):
    if isinstance(expr, (TensorFunction, TensorProductFunction,)):
        return expr
    if isinstance(expr, Mul):
        # Look for a Tensor here, and multiply everything else by it.
        # First, try to ReduceExpr on everything that will be going
        # into this calculation in hopes of getting some tensors.
        args = list(o if o.is_Atom or isinstance(o, (TensorFunction, TensorProductFunction,))
                    else ReduceExpr(o)
                    for o in expr.args)
        tensors = prod(t for t in args if isinstance(t, (TensorFunction, TensorProductFunction,)))
        others = prod(o for o in args if not isinstance(o, (TensorFunction, TensorProductFunction,)))
        if tensors==1:
            return others
        else:
            return tensors*others
    if isinstance(expr, Add):
        return sum(ReduceExpr(arg) for arg in expr.args)
    return expr

####################################
### First, a few vector thingies ###
####################################
class VectorFunction(Function):
    """\
    This is just a base class for deriving other vectors from.

    You probably won't need to use this class directly; it is just the
    base class for the class created in VectorFunction factory below.

    This used to be more important because there were other subclasses
    of this, but now it's mostly just for convenient separation of a
    few methods, to clarify what's going on in the factory.

    """
    _op_priority = 1000000.0
    components = None
    coefficient = 1
    @property
    def _is_vector(self):
        return True
    def __eq__(self, other):
        try:
            if(len(self.components)!=len(other.components)):
                return False
            for c1,c2 in zip(self, other):
                if c1!=c2:
                    return False
        except:
            return False
        return True
    def __ror__(self,other):
        return self.__or__(other)
    def __or__(self,other):
        """
        In keeping with the notation of various other packages
        (most notably sympy.galgebra), contraction is denoted
        by the bitwise `or` operator.  So the contraction of
        vectors `v` and `w` is just `v|w`.  This notation will
        be used in the tensor classes also.
        """
        return self.coefficient*other.coefficient * sum( s*o for s,o in zip(self, other) )
    def __iter__(self):
        for c in self.__class__.components: yield c
    def __div__(self, other):
        return self.__mul__(sympify(1)/other)
    def __mul__(self, other):
        if(other==1):
            return self
        if(other==0):
            return 0
        if(hasattr(other, '_is_tensor_product') or hasattr(other, '_is_tensor')):
            return NotImplemented
        if(hasattr(other, '_is_vector')):
            return TensorProduct(self, other)
        return Vector(str(other)+'*'+self.name,
                      latex(other)+r'\,'+self.latex_name,
                      [other*c for c in self])
    def __rmul__(self, other):
        if(other==1):
            return self
        if(other==0):
            return 0
        if(hasattr(other, '_is_tensor_product') or hasattr(other, '_is_tensor')):
            return NotImplemented
        if(hasattr(other, '_is_vector')):
            return TensorProduct(other, self)
        return Vector(self.name+'*'+str(other),
                      self.latex_name+r'\,'+latex(other),
                      [c*other for c in self])
    def _eval_derivative(self, s):
        """Return derivative of the function with respect to `s`.

        Note that we must take the chain rule into account in this
        function, but not in `fdiff`.

        Note that `diff` is the general function for evaluating
        derivatives; `Derivative` results in unevaluated derivatives.

        """
        # print(self.fdiff(1))
        # print(self.args[0].diff(s))
        return self.fdiff(1) * self.args[0].diff(s)
    def fdiff(self, argindex=1):
        """Returns the first derivative of the function.

        To be precise, this is the first derivative with respect to
        the function's argument, whatever that argument may be.  And
        in particular, we don't need to account for the chain rule.

        Note that `diff` is the general function for evaluating
        derivatives; `Derivative` results in unevaluated derivatives.

        """
        from sympy.core.function import ArgumentIndexError
        # print("VectorFunction.fdiff")
        if (argindex != 1):
            raise ArgumentIndexError(self, argindex)
        # print("Returning from VectorFunction.fdiff")
        V = Vector('d'+self.name+'d'+str(self.args[0]),
                   DifferentiateString(self.latex_name, self.args[0]),
                   [diff(c, self.args[0]) for c in self])
        if V==0:
            return 0
        return V(self.args[0])
    def subs(self, *args, **kwargs):
        if(len(args)==0):
            return self
        if(len(args)==1 and isinstance(args[0], list) and not args[0]):
            return self
        V = Vector(self.name+'.subs({0}, {1})'.format(args, kwargs),
                   LatexSubs(self.latex_name, args, kwargs),
                   [sympify(c).subs(*args, **kwargs) for c in self])
        if V==0: return 0
        return V(self.args[0])
    def __str__(self):
        return self.name
    def __repr__(self):
        return self.name.translate(TranslationTable)
    def _latex(self, printer):
        return self.latex_name
    def _repr_latex_(self):
        return self.latex_name

_Vector_count = 0
def Vector(Name, LatexName, ComponentFunctions):
    """Create a new vector function

    This function creates a class that is a subclass of
    `VectorFunction`, which is itself a subclass of
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
    if(ComponentFunctions == [0,]*len(ComponentFunctions)):
        return 0
    ## Now, create the object and set its data.  Because of sympy's
    ## caching, vectors with different data need to be created as
    ## classes with different names.  So we just create a lighweight
    ## subclass with a unique name (the number at the end gets
    ## incremented every time we construct a vector).
    # global _Vector_count
    # ThisVectorFunction = type('VectorFunction_'+str(_Vector_count),
    #                           (VectorFunction,), {})
    # _Vector_count += 1
    ThisVectorFunction = type(Name, (VectorFunction,), {})
    ThisVectorFunction.name = Name
    ThisVectorFunction.latex_name = LatexName
    ThisVectorFunction.components = list(ComponentFunctions)
    return ThisVectorFunction


xHat = Vector('xHat', r'\hat{x}', [1,0,0])
yHat = Vector('yHat', r'\hat{y}', [0,1,0])
zHat = Vector('zHat', r'\hat{z}', [0,0,1])


################################
### Now, the tensor products ###
################################
class TensorProductFunction(Function):
    _op_priority = 1000001.0
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

    def __ror__(self,other):
        return self.__or__(other)
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

    def __div__(self, other):
        return self.__mul__(sympify(1)/other)

    def __mul__(self, B):
        """
        Return the scalar or tensor product
        """
        # print('TensorProductFunction.__mul__<{0},{1}>({2},{3})'.format(type(self), type(B), self,B))
        if(hasattr(B, '_is_tensor') and B._is_tensor):
            # Fall back to Tensor.__rmul__ by doing this:
            return NotImplemented
        if(isinstance(B, TensorProductFunction)):
            # Do tensor product
            # print('TensorProductFunction.__mul__ return 1')
            return TensorProduct(self.vectors + B.vectors,
                                 coefficient=simplify( self.coefficient * B.coefficient ),
                                 symmetric = self.symmetric)
        elif(isinstance(B, VectorFunction)):
            # print('TensorProductFunction.__mul__ return 2')
            return TensorProduct(self.vectors + [B,],
                                 coefficient=self.coefficient,
                                 symmetric = self.symmetric)
        else:
            try:
                if(simplify(B)==0): return 0
            except:
                pass
            # Otherwise, try scalar multiplication
            # print('TensorProductFunction.__mul__ return 3')
            return TensorProduct(self.vectors,
                                 coefficient=self.coefficient*B,
                                 symmetric = self.symmetric)

    def __rmul__(self, B):
        """
        Return the scalar or tensor product
        """
        # print('TensorProductFunction.__rmul__<{0},{1}>({2},{3})'.format(type(self), type(B), self, B))
        if(hasattr(B, '_is_tensor') and B._is_tensor):
            # Fall back to Tensor.__rmul__ by doing this:
            return NotImplemented
        if(isinstance(B, TensorProductFunction)):
            # Do tensor product
            # print('TensorProductFunction.__rmul__ return 1')
            return TensorProduct(B.vectors+self.vectors,
                                 coefficient=simplify( B.coefficient * self.coefficient ),
                                 symmetric = self.symmetric)
        elif(isinstance(B, VectorFunction)):
            # print('TensorProductFunction.__rmul__ return 2')
            return TensorProduct([B,] + self.vectors,
                                 coefficient=self.coefficient,
                                 symmetric = self.symmetric)
        else:
            try:
                if(simplify(B)==0): return 0
            except:
                pass
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

        Note that `diff` is the general function for evaluating
        derivatives; `Derivative` results in unevaluated derivatives.

        """
        return self.fdiff(1) * self.args[0].diff(s)

    def fdiff(self, argindex=1):
        """Returns the first derivative of the function.

        To be precise, this is the first derivative with respect to
        the function's argument, whatever that argument may be.  And
        in particular, we don't need to account for the chain rule.

        Note that `diff` is the general function for evaluating
        derivatives; `Derivative` results in unevaluated derivatives.

        """
        from sympy.core.function import ArgumentIndexError
        if (argindex != 1):
            raise ArgumentIndexError(self, argindex)
        TP = TensorProduct(list(self.vectors),
                           coefficient = diff(self.coefficient, self.args[0]),
                           symmetric = self.symmetric) \
            + sum( [ TensorProduct([diff(t, self.args[0]) if i==j else t for j,t in enumerate(self)],
                                   coefficient = self.coefficient, symmetric=self.symmetric)
                     for i in range(self.rank) ] )
        try:
            return TP.compress()
        except AttributeError:
            return TP

    def subs(self, *args, **kwargs):
        TP = TensorProduct([c.subs(*args, **kwargs) for c in self],
                           coefficient = self.coefficient.subs(*args, **kwargs), symmetric=self.symmetric)
        try:
            return TP.compress()
        except:
            try:
                if TP==0: return 0
            except:
                pass
            return TP

    def __str__(self):
        if(self.coefficient==1):
            return DelimitString('*'.join([str(v) for v in self]), latex=False)
        return DelimitString( DelimitString(str(self.coefficient), latex=False)
                              + '*' + '*'.join([str(v) for v in self]), latex=False )

    def __repr__(self):
        if(self.coefficient==1):
            return DelimitString('*'.join([repr(v) for v in self]), latex=False).translate(TranslationTable)
        return DelimitString( DelimitString(str(self.coefficient), latex=False)
                              + '*' + '*'.join([repr(v) for v in self]), latex=False ).translate(TranslationTable)

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
        # Since TensorProducts are multiplicative, the empty object
        # should be 1 (or whatever coefficient was passed, if any)
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

    ## Now, make sure none of the input vectors are zero
    for v in vectors:
        if v==0:
            return 0

    ## Finally, create the object and set its data.  Because of
    ## sympy's caching, tensor products with different data need to be
    ## created as classes with different names.  So we just create a
    ## lighweight subclass with a unique name (the number at the end
    ## gets incremented every time we construct a tensor product).
    global _TensorProduct_count
    ThisTensorProductFunction = type('TensorProductFunction_'+str(_TensorProduct_count),
                                     (TensorProductFunction,), {})
    _TensorProduct_count += 1
    # print('About to construct a tensor with args ',
    #       tuple( set( flatten( [v.args for v in vectors] ) ) ),
    #       input_vectors,
    #       kwargs,
    #       vectors,
    #       [v.args for v in vectors] )
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
    _op_priority = 1000002.0

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
        # def debugcompress(s):
        #     print(s)
        # debugcompress("\ncompressing")
        # debugcompress(self)
        removed_elements = []
        for i in range(len(self.tensor_products)):
            # debugcompress("\t%s"%i)
            if(i in removed_elements):
                # debugcompress("\t\tskipping %s"%i)
                continue
            for j in range(i+1,len(self.tensor_products)):
                # debugcompress("\t\t%s"%j)
                if(j in removed_elements):
                    # debugcompress("\t\t\tskipping %s"%j)
                    continue
                if self.tensor_products[i].has_same_basis_element(self.tensor_products[j]):
                    # debugcompress("Removing {0} because {1} is already here".format(j,i))
                    removed_elements.append(j)
                    NewCoefficient = simplify( self.tensor_products[i].coefficient + self.tensor_products[j].coefficient )
                    if(NewCoefficient==0):
                        # debugcompress("Also removing {0} because {1} cancelled it out".format(i,j))
                        removed_elements.append(i)
                        break
                    self.tensor_products[i] = TensorProduct(self.tensor_products[i].vectors,
                                                            coefficient=NewCoefficient,
                                                            symmetric=self.tensor_products[i].symmetric)
                # debugcompress("\t\t\tfinished %s"%j)
            # debugcompress("\t\tfinished j loop for %s"%i)
            if removed_elements:
                if(self.tensor_products[i].coefficient==0 or self.tensor_products[i]==0):
                    removed_elements += [i]
                # debugcompress("Removing {0} because {1} is already here".format(removed_elements,i))
            # debugcompress("\t\tfinished %s"%i)
        # debugcompress("compression results in {0} elements".format(len(self.tensor_products) - len(removed_elements)))
        self.tensor_products = list(t_p for i,t_p in enumerate(self) if i not in removed_elements and t_p!=0 and t_p.coefficient!=0)
        # debugcompress("compression results in {0} elements:\n{1}\n.\n".format(len(self.tensor_products), self.tensor_products))
        if not self.tensor_products:
            # debugcompress("compressed to 0")
            return 0
        # debugcompress("compressed to {0}\n".format(self))
        return self

    def __add__(self, T):
        if(T==0):
            return self
        if(self.rank==0):
            return T
        if not hasattr(T, 'rank'):
            try:
                T = ReduceExpr(T)
                # T = simplify(T)
                # T = T.doit()
                # T = ReduceExpr(T)
            except Exception as e:
                print('Failed to ReduceExpr({0})\n"{1}"\n\n'.format(T, e))
                pass
            if not hasattr(T, 'rank'):
                print("You probably should never be trying to add something without a rank to a tensor!!!")
                print("Trying with {2}({3}) T={0}; self={1}".format(T,self, T.func, T.args))
                return NotImplemented
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

    def __sub__(self, T):
        if(T==0):
            return self
        if(self.rank==0):
            return -T
        if not hasattr(T, 'rank'):
            try:
                T = ReduceExpr(T)
                # T = simplify(T)
                # T = T.doit()
            except:
                pass
            if not hasattr(T, 'rank'):
                print("This is bad!!!  You probably should never be trying to subtract something without a rank from a tensor!!!\nT={0}; self={1}".format(T,self))
                return NotImplemented
        if(T.rank==0):
            return self
        if(T.rank!=self.rank):
            raise ValueError("Cannot add rank-{0} tensor to rank-{1} tensor.".format(T.rank, self.rank))
        if(isinstance(T, TensorFunction)) :
            return Tensor(self.tensor_products + simplify(-1*T).tensor_products)
        elif(isinstance(T, TensorProductFunction)) :
            return Tensor(self.tensor_products + [simplify(-1*T),])

    def __ror__(self,other):
        return self.__or__(other)
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

    def __div__(self, other):
        print('TensorFunction.__div__<{0},{1}>({2},{3})'.format(type(self), type(B), self, other))
        return Tensor(list(t_p/other for t_p in self))

    def __mul__(self, other):
        # print('TensorFunction.__rmul__<{0},{1}>({2},{3})'.format(type(self), type(other), self, other))
        if(isinstance(other, TensorFunction)):
            # print('TensorFunction.__mul__ return 1')
            return Tensor(list(t_pA*t_pB for t_pA in self for t_pB in other))
        else:
            # print('TensorFunction.__mul__ return 2')
            return Tensor(list(t_p*other for t_p in self))

    def __rmul__(self, other):
        # print('TensorFunction.__rmul__<{0},{1}>({2},{3})'.format(type(self), type(other), self, other))
        if(isinstance(other, TensorFunction)):
            # print('TensorFunction.__rmul__ return 1')
            return Tensor(list(t_pB*t_pA for t_pA in self for t_pB in other))
        else:
            # print('TensorFunction.__rmul__ return 2')
            return Tensor(list(t_p*other for t_p in self))

    def _eval_derivative(self, s):
        """Return derivative of the function with respect to `s`.

        Note that we must take the chain rule into account in this
        function, but not in `fdiff`.

        Note that `diff` is the general function for evaluating
        derivatives; `Derivative` results in unevaluated derivatives.

        """
        return self.fdiff(1) * self.args[0].diff(s)

    def fdiff(self, argindex=1):
        """Returns the first derivative of the function.

        To be precise, this is the first derivative with respect to
        the function's argument, whatever that argument may be.  And
        in particular, we don't need to account for the chain rule.

        Note that `diff` is the general function for evaluating
        derivatives; `Derivative` results in unevaluated derivatives.

        """
        from sympy.core.function import ArgumentIndexError
        if (argindex != 1):
            raise ArgumentIndexError(self, argindex)
        T = sum(diff(t_p, self.args[0]) for t_p in self)
        try:
            return T.compress()
        except AttributeError:
            return T

    def subs(self, *args, **kwargs):
        T = Tensor([c.subs(*args, **kwargs) for c in self])
        try:
            return T.compress()
        except:
            try:
                if T==0: return 0
            except:
                pass
            return T

    def __str__(self):
        return DelimitString( '\n +'.join([str(t_p) for t_p in self]), latex=False)

    def __repr__(self):
        return DelimitString( '\n +'.join([repr(t_p) for t_p in self]), latex=False).translate(TranslationTable)

    def _latex_str_(self):
        return '&' + DelimitString( r' \right. \nonumber \\&\quad \left. + '.join(
            [t_p._latex_str_() for t_p in self]) ) + r'\\'

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
        tensor_products = list(t_p for t_p in tensor_products[0].tensor_products if t_p!=0)
    elif(len(tensor_products)==1 and isinstance(tensor_products[0], TensorProductFunction)) :
        tensor_products = [tensor_products[0],]
    else:
        if(len(tensor_products)==1 and isinstance(tensor_products[0], list)):
            tensor_products = tensor_products[0]
        tensor_products = flatten(list(t_p for t_p in tensor_products if t_p!=0))
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




def contract(Expr):
    pass
