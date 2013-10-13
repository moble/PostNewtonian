from collections import OrderedDict
from sympy import Symbol, Function

class _PNSymbol(Symbol) :
    """
    This is the basic object created by calls to `AddVariable`,
    etc., and is a simple subclass of python.Symbol, as described
    above.

    The method `__new__` is always called first, since this is an
    immutable object, which creates the object, allocating memory
    for it.  Since `__new__` actually returns an object,
    `__init__` is then called with the same arguments as
    `__new__`.  This is why we throw away the three custom
    arguments in `__new__`, and throw away the rest in `__init__`.

    """
    def __new__(cls, name, constant, fundamental, substitution, atoms, **assumptions) :
        from sympy import Symbol
        return Symbol.__new__(cls, name, **assumptions)
    def __init__(self, name, constant, fundamental, substitution, atoms, **kwargs) :
        self.constant = constant
        self.fundamental = fundamental
        self.substitution = substitution
        if atoms:
            self.atoms = atoms
        else:
            try:
                self.atoms = self.substitution.atoms(Symbol)
            except AttributeError:
                self.atoms = None

class PNVariablesCollection(OrderedDict) : # subclass of OrderedDict
    """Subclass of `OrderedDict` to hold PN variables, each of which is a subclasses sympy `Symbol`

    This class has a few extra functions to add variables nicely, and
    include them in the calling scope.  The PN variables are just like
    sympy `Symbol`s, except they also remember three things:

      * `constant` (for things like mass)
      * `fundamental` (for things that aren't defined in terms of other things)
      * `substitutions` (re-express non-fundamental objects in terms of fundamentals)

    """
    def __init__(self, *args):
        OrderedDict.__init__(self, *args)
    def _AddVariable(self, name, constant=False, fundamental=False, substitution=None, atoms=None, **args) :
        from inspect import currentframe
        from sympy import Basic, FunctionClass
        frame = currentframe().f_back.f_back
        try:
            args['constant'] = constant
            args['fundamental'] = fundamental
            args['substitution'] = substitution
            args['atoms'] = atoms
            cls = args.pop('cls', _PNSymbol)
            sym = cls(name, **args)
            if sym is not None:
                if isinstance(sym, Basic):
                    frame.f_globals[sym.name] = sym
                elif isinstance(sym, FunctionClass):
                    frame.f_globals[sym.__name__] = sym
        finally:
            del frame
        self[sym] = name
        return sym
    def AddVariable(self, name, constant=False, fundamental=False, substitution=None, atoms=None, **args) :
        return self._AddVariable(name, constant, fundamental, substitution, atoms, **args)
    def AddBasicConstants(self, names, **args) :
        from re import split
        names = split(',| ', names)
        for name in names :
            if name :
                self._AddVariable(name, constant=True, substitution=None, **args)
    def AddBasicVariables(self, names, **args) :
        from re import split
        names = split(',| ', names)
        for name in names :
            if name :
                self._AddVariable(name, constant=False, substitution=None, **args)
    def AddDerivedConstant(self, name, substitution, atoms=None, **args) :
        self._AddVariable(name, constant=True, fundamental=False, substitution=substitution, atoms=atoms, **args)
    def AddDerivedVariable(self, name, substitution, atoms=None, **args) :
        self._AddVariable(name, constant=False, fundamental=False, substitution=substitution, atoms=atoms, **args)
