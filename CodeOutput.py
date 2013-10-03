def Expression(Quantity, ns) :
    if (isinstance(Quantity, list) and len(Quantity)==2) :
	if (isinstance(Quantity[1], basestring)) :
	    try :
		# test if it names something that exists
		eval(Quantity[1], ns)
		return Quantity[1]
	    except NameError :
		return None
	elif ('sympy' in type(Quantity[1]).__name__) :
	    return str(Quantity[1])
    elif (isinstance(Quantity, basestring)) :
	if (' ' not in Quantity) :
	    try :
		# test if it names a variable that exists
		eval(Quantity, ns)
		return Quantity
	    except NameError :
		return None
	else :
	    return None
    else :
	return None


def FindAtoms(Expressions, PNVariables, ns) :
    from sympy import Symbol
    Atoms = set([])
    for Expression in Expressions :
	Atoms.update(eval(Expression, ns).atoms(Symbol))
    LastAtomsLength = 0
    while(len(Atoms) != LastAtomsLength) :
	LastAtomsLength = len(Atoms)
	for Atom in list(Atoms) :
	    if (Atom.substitution) :
		Atoms.update(Atom.substitution.atoms(Symbol))
    BasicConstantAtoms = []
    BasicVariableAtoms = []
    DerivedConstantAtoms = []
    DerivedVariableAtoms = []
    for key in PNVariables.iterkeys() :
	if (key in Atoms) :
	    if (key.constant) :
		if (key.substitution) :
		    DerivedConstantAtoms += [key]
		else :
		    BasicConstantAtoms += [key]
	    else :
		if (key.substitution) :
		    DerivedVariableAtoms += [key]
		else :
		    BasicVariableAtoms += [key]
    return BasicConstantAtoms, BasicVariableAtoms, DerivedConstantAtoms, DerivedVariableAtoms

def CCodeOutput(Quantities, PNVariables, ns, Utilities=[], Indent=[2,4,4,4]) :
    """
    Return four strings for C/C++ code compilations.

    The four strings are:

      * `Declarations`: Declarations of the variables
      * `Initializations`: Initializer list of the variables
      * `Evaluations`: Re-evaluation of all non-constant variables
      * `Computations`: Computation of the requested quantities themselves
    """
    from sympy import ccode, horner, N
    from textwrap import TextWrapper
    wrapper = TextWrapper(width=120)
    def codify(s, **args) :
        try :
            s = horner(s)
        except PolynomialError :
            pass
        return ccode(N(s), **args)

    # Accept a single string as input
    if (not isinstance(Quantities, list)) :
	Quantities = [Quantities]

    # Get a list of names of sympy expressions in the input
    Expressions = [E for Quantity in Quantities+Utilities for E in [Expression(Quantity, ns)] if E]
    BasicConstantAtoms, BasicVariableAtoms, DerivedConstantAtoms, DerivedVariableAtoms = FindAtoms(Expressions, PNVariables, ns)

    # basic const declarations
    names = ['{0}'.format(PNVariables[atom]) for atom in BasicConstantAtoms]
    wrapper.initial_indent = ' '*Indent[0] + "const double "
    wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
    Declarations = wrapper.fill(', '.join(names)) + ";\n"
    # basic non-const declarations
    names = ['{0}'.format(PNVariables[atom]) for atom in BasicVariableAtoms]
    wrapper.initial_indent = ' '*Indent[0] + "double "
    wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
    Declarations += wrapper.fill(', '.join(names)) + ";\n"
    # variable const declarations
    names = ['{0}'.format(PNVariables[atom]) for atom in DerivedConstantAtoms]
    wrapper.initial_indent = ' '*Indent[0] + "const double "
    wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
    Declarations += wrapper.fill(', '.join(names)) + ";\n"
    # variable non-const declarations
    names = ['{0}'.format(PNVariables[atom]) for atom in DerivedVariableAtoms]
    wrapper.initial_indent = ' '*Indent[0] + "double "
    wrapper.subsequent_indent = ' '*len(wrapper.initial_indent)
    Declarations += wrapper.fill(', '.join(names)) + ";"


    # Initializations
    names = ['{0}({0}_in)'.format(PNVariables[atom]) for atom in BasicConstantAtoms]
    names += ['{0}({0}_0)'.format(PNVariables[atom]) for atom in BasicVariableAtoms]
    names += ['{0}({1})'.format(PNVariables[atom],
				codify(atom.substitution)) for atom in DerivedConstantAtoms]
    names += ['{0}({1})'.format(PNVariables[atom],
				codify(atom.substitution)) for atom in DerivedVariableAtoms]
    wrapper.initial_indent = ' '*Indent[1]
    wrapper.subsequent_indent = wrapper.initial_indent
    Initializations = wrapper.fill(', '.join(names))


    # Evaluations
    wrapper.initial_indent = ' '*Indent[2]
    wrapper.subsequent_indent = wrapper.initial_indent+' '*4
    names = [wrapper.fill('{0} = {1};'.format(PNVariables[atom],
					      codify(atom.substitution))) for atom in DerivedVariableAtoms]
    Evaluations = '\n'.join(names)


    # Computations
    wrapper.initial_indent = ' '*Indent[3]
    wrapper.subsequent_indent = wrapper.initial_indent+' '*4
    names = []
    for Quantity in Quantities :
	if (isinstance(Quantity, list) and len(Quantity)==2) :
	    if (isinstance(Quantity[1], basestring)) :
		names += [wrapper.fill(str(Quantity[0])+codify(eval(Quantity[1], ns)))]
	    elif ('sympy' in type(Quantity[1]).__name__) :
		names += [wrapper.fill(str(Quantity[0])+codify(Quantity[1]))]
	elif (isinstance(Quantity, basestring)) :
	    if (' ' not in Quantity) :
		names += [wrapper.fill(codify(eval(Quantity, ns), assign_to='const double '+Quantity))]
	    else :
		names += [wrapper.fill(str(Quantity))]
	else :
	    names += [wrapper.fill(str(Quantity))]
    Computations = '\n'.join(names)

    return Declarations, Initializations, Evaluations, Computations
