from collections import OrderedDict
from sympy import Symbol, Function

class PNSymbol(Symbol) :
    """This is the basic object created by calls to `AddVariable`, etc.,
    and is a simple subclass of python.Symbol, as described above.

    The method `__new__` is always called first, since this is an
    immutable object, which creates the object, allocating memory for
    it.  Since `__new__` actually returns an object, `__init__` is
    then called with the same arguments as `__new__`.  This is why we
    throw away the custom arguments in `__new__`, and throw away the
    rest in `__init__`.

    """
    def __new__(cls, name, constant, fundamental, substitution, substitution_atoms, datatype, **assumptions) :
        from sympy import Symbol
        return Symbol.__new__(cls, name, **assumptions)
    def __init__(self, name, constant, fundamental, substitution, substitution_atoms, datatype, **kwargs) :
        if not fundamental and isinstance(substitution, basestring) and not substitution_atoms:
            raise ValueError('Either `substitution` must be a sympy expression, '
                             +'or `substitution_atoms` must be non-empty for derived variables.')
        self.constant = constant
        self.fundamental = fundamental
        self.substitution = substitution
        if substitution_atoms:
            self.substitution_atoms = substitution_atoms
        else:
            try:
                self.substitution_atoms = self.substitution.atoms(Symbol)
            except AttributeError:
                self.substitution_atoms = None
        self.datatype = datatype
    def ccode(**args):
        from sympy import ccode, horner, N
        if self.fundamental:
            return str(self)
        if isinstance(self.substitution, basestring):
            return self.substitution
        code = self.substitution
        try:
            code=N(code)
        except:
            pass
        try:
            code=horner(code, wrt=args.pop('wrt', None))
        except:
            pass
        return ccode(code, **args)

class PNCollection(OrderedDict) : # subclass of OrderedDict
    """Subclass of `OrderedDict` to hold PN variables, each of which is a
    subclass sympy `Symbol`

    This class has a few extra functions to add variables nicely, and
    include them in the calling scope.  The PN variables are just like
    sympy `Symbol`s, except they also remember these things:

      * `constant` (boolean: doesn't need to be updated once it's defined)
      * `fundamental` (boolean: not defined in terms of other things)
      * `substitution` (string or sympy expression: re-express non-fundamental object in terms of fundamentals)
      * `substitution_atoms` (atomic sympy objects in terms of which the variable is defined)
      * `datatype` (optional name for datatype of variable in code output;
                    if None, the basic real datatype of the output language is assumed)

    """
    def __init__(self, *args):
        OrderedDict.__init__(self, *args)
    def _AddVariable(self, name, **args) :
        from inspect import currentframe
        from sympy import Basic, FunctionClass
        frame = currentframe().f_back.f_back
        try:
            sym = PNSymbol(name, **args)
            if sym is not None:
                if isinstance(sym, Basic):
                    frame.f_globals[sym.name] = sym
                elif isinstance(sym, FunctionClass):
                    frame.f_globals[sym.__name__] = sym
        finally:
            del frame
        self[sym] = name
        return sym
    def AddVariable(self, name, constant=False, fundamental=False, substitution=None, substitution_atoms=None, datatype=None, **args) :
        return self._AddVariable(name, constant=constant, fundamental=fundamental,
                                 substitution=substitution, substitution_atoms=substitution_atoms, datatype=datatype, **args)
    def AddBasicConstants(self, names, datatypes=None, **args) :
        from re import split
        names = split(',| ', names)
        args['constant'] = True
        args['fundamental'] = True
        args['substitution'] = args.pop('substitution', None)
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        for name in names :
            if name :
                self._AddVariable(name, **args)
    def AddBasicVariables(self, names, **args) :
        from re import split
        names = split(',| ', names)
        args['constant'] = False
        args['fundamental'] = True
        args['substitution'] = args.pop('substitution', None)
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        for name in names :
            if name :
                self._AddVariable(name, **args)
    def AddDerivedConstant(self, name, substitution, **args) :
        args['constant'] = True
        args['fundamental'] = False
        args['substitution'] = substitution
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        self._AddVariable(name, **args)
    def AddDerivedVariable(self, name, substitution, **args) :
        args['constant'] = False
        args['fundamental'] = False
        args['substitution'] = substitution
        args['substitution_atoms'] = args.pop('substitution_atoms', None)
        args['datatype'] = args.pop('datatype', None)
        self._AddVariable(name, **args)
