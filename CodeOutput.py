class CodeConstructor:
    """Contains lists of variables and expressions to be written as code.

    `CodeConstructor` objects contain:
      1) An ordered list of atoms for the code to use
      2) A PNCollection of PNSymbol objects
      3) A PNCollection of expressions to be calculated

    Once the `CodeConstructor` is initialized with these objects, it
    can be used to construct various types of code.  For example, the
    `CppDeclarations` method will output a list of declarations of the
    atoms.  Similar methods are available for function input
    arguments, class initializer lists, and the final evaluations
    needed to calculate the input `Expressions`.

    Support for other languages or constructions can be added by
    adding more method functions to this class.

    Note that it is generally necessary to obey a strict ordering for
    defining variables.  The functions of this class assume that the
    ordering in which the variables were defined in python should
    remain the same in the output code.  Because of the structure of
    the `PNCollection` objects, it should be hard to define the
    variables out of order, so this should not require anything from
    the user.  (That is why `PNCollection` is a subclass of the basic
    `OrderedDictionary` object.)  However, if new functions are added
    here, they must obey that ordering.

    """
    def __init__(self, Variables, Expressions):
        AtomSet = set([])
        self.Variables = Variables.copy()
        self.Expressions = Expressions.copy()
        for Expression in self.Expressions:
            AtomSet.update(Expression.substitution_atoms)
        LastAtomsLength = 0
        while(len(AtomSet) != LastAtomsLength):
            LastAtomsLength = len(AtomSet)
            for Atom in list(AtomSet):
                if (Atom.substitution_atoms):
                    AtomSet.update(Atom.substitution_atoms)
        self.Atoms = []
        for sym in self.Variables:
            if sym in AtomSet:
                self.Atoms.append(sym)

    @staticmethod
    def const(e):
        if e.constant:
            return 'const '
        else:
            return ''

    @staticmethod
    def dtype(e):
        if e.datatype:
            return e.datatype
        else:
            return 'double'

    def AddDependencies(self, Expressions):
        AtomSet = set([])
        for Expression in Expressions:
            if (Expression.substitution_atoms):
                AtomSet.update(Expression.substitution_atoms)
        LastAtomsLength = 0
        while(len(AtomSet) != LastAtomsLength):
            LastAtomsLength = len(AtomSet)
            for Atom in list(AtomSet):
                if (Atom.substitution_atoms):
                    AtomSet.update(Atom.substitution_atoms)
        OldAtoms = self.Atoms[:]
        self.Atoms = []
        for sym in self.Variables:
            if sym in AtomSet or sym in OldAtoms:
                self.Atoms.append(sym)

    def CppDeclarations(self, Indent=4):
        """Create declaration statements for C++

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the declaration
        list should be

            const double m1, m2;
            double t, x;

        The code knows which atoms need to be declared at the
        beginning, and which ones should be `const`, for example.  For
        C++, the default datatype is `double`; if the atom was created
        with a different datatype, that will be used appropriately.

        """
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ''
        wrapper.subsequent_indent = ''
        datatype = ''
        Declarations = ''
        Names = []
        for atom in self.Atoms:
            thisdatatype = CodeConstructor.const(atom) + CodeConstructor.dtype(atom) + ' '
            if thisdatatype != datatype:
                if Names:
                    Declarations += wrapper.fill(', '.join(Names)) + ";\n"
                Names = []
                datatype = thisdatatype
                wrapper.initial_indent = ' '*Indent + thisdatatype
                wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
            Names.append(self.Variables[atom])
        if Names:
            Declarations += wrapper.fill(', '.join(Names)) + ";\n"
        return Declarations.rstrip()

    def CppInputArguments(self, Indent=12):
        """Create basic input arguments for C++

        The fundamental variables are listed, along with their data
        types and `const` if the variable is constant.  This would be
        an appropriate string to represent the input arguments for a
        function or class constructor to calculate the `Expressions`
        of this CodeConstructor object.

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the input argument
        list should be

            const double m1_i, const double m2_i, double t_i, double x_i

        """
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ' '*Indent
        wrapper.subsequent_indent = wrapper.initial_indent
        InputArguments = ['const {0} {1}_i'.format(self.dtype(atom), self.Variables[atom])
                          for atom in self.Atoms if atom.fundamental]
        return wrapper.fill(', '.join(InputArguments)).lstrip()

    def CppInitializations(self, Indent=4):
        """Create initialization list for C++

        For example, if the `Variables` object contains atoms m1, m2,
        t, and x referred to in the `Expressions` object, where m1 and
        m2 are constant, and t and x are variables, the initialization
        list should be

            m1(m1_i), m2(m2_i), t(t_i), x(x_i)

        The quantities m1_i, etc., appear in the input-argument list
        output by the method `CppInputArguments`.

        """
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ' '*Indent
        wrapper.subsequent_indent = wrapper.initial_indent
        def Initialization(atom):
            if atom.fundamental:
                return '{0}({0}_i)'.format(self.Variables[atom])
            else:
                return '{0}({1})'.format(self.Variables[atom], atom.ccode())
        Initializations  = [Initialization(atom) for atom in self.Atoms]
        return wrapper.fill(', '.join(Initializations))

    def CppEvaluations(self, Indent=4):
        """Evaluate all derived variables in C++

        This function uses the `substitution` expressions for the
        derived variables.  This output is appropriate for updating
        the values of the variables at each step of an integration,
        for example.

        """
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ' '*Indent
        wrapper.subsequent_indent = wrapper.initial_indent + '  '
        return '\n'.join([wrapper.fill('{0} = {1};'.format(self.Variables[atom], atom.ccode()))
                          for atom in self.Atoms if not atom.fundamental and not atom.constant])

    def CppEvaluateExpressions(self, Indent=4, Expressions=None):
        """Declare and define the `Expressions` for C++

        The output of this function declares are defines the
        `Expressions` as individual variables.  An optional dictionary
        of expressions allows just a subset of this object's
        expressions to be output; if this argument is not present, all
        will be output.

        """
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ' '*Indent
        wrapper.subsequent_indent = wrapper.initial_indent+'  '
        Evaluations = []
        if not Expressions:
            Expressions=self.Expressions
        for Expression in Expressions:
            Evaluations.append(wrapper.fill('{0}{1} {2} = {3};'.format(self.const(Expression), self.dtype(Expression),
                                                                       Expressions[Expression], Expression.ccode())))
        return '\n'.join(Evaluations)

    def CppExpressionsAsFunctions(self, Indent=4, Expressions=None):
        """Define functions to calculate the `Expressions` in C++

        The output of this function gives C++ functions to calculate
        the `Expressions`, assuming the functions are member methods
        in a class, and so can access the atoms of the expression
        without explicit arguments.  An optional dictionary of
        expressions allows just a subset of this object's expressions
        to be output; if this argument is not present, all will be
        output.

        """
        def dtype(e):
            if e.datatype:
                return e.datatype
            else:
                return 'double'
        from textwrap import TextWrapper
        wrapper = TextWrapper(width=120)
        wrapper.initial_indent = ' '*Indent + '  return '
        wrapper.subsequent_indent = ' '*Indent + '    '
        Evaluations = []
        if not Expressions:
            Expressions=self.Expressions
        for Expression in Expressions:
            Evaluations.append(
                ' '*Indent + dtype(Expression) + ' ' + Expressions[Expression] + '() {\n'
                + wrapper.fill(Expression.ccode())
                + ';\n' + ' '*Indent + '}'
            )
        return '\n'.join(Evaluations)
