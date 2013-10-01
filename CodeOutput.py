# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def CCodeOutput(Quantities, Indent=4) :
    """
    Print C/C++ code to evaluate the given sympy quantities.

    The input argument should be a list of strings or a list of pairs
    of strings (or just one single string).  For each item in
    Quantities, the following output is given:

    * If the item is a pair of strings, the first is taken to be the
      constructing statement (e.g., `const std::vector<double> x = `
      equation), while the second is evaluated and placed on the right
      of the first string, followed by a semicolon.

    * If the item is a single string containing no space, it is
      defined as a new `const double` with the given name, and set
      equal to that name evaluated.

    * Everything else is simply converted to a string and added to the
      output.

    Note that all atoms (basic sympy symbols) are re-calculated if
    they have entries in the dictionary `BasicSubstitutions`.  Some of
    these may be unnecessary, and can be removed from the output before
    compilation.

    >>> from sympy.abc import x,y
    >>> x = y**2
    >>> print CCodeOutput(['x', 'say hello'])

    const double x = pow(y, 2);
    say hello
    """
    from sympy import Symbol, ccode, horner
    ReturnList = ['']
    Atoms = set([])
    if (not isinstance(Quantities, list)) :
        Quantities = [Quantities]
    for Quantity in Quantities :
        if (isinstance(Quantity, list) and len(Quantity)==2) :
            if (isinstance(Quantity[1], basestring)) :
                Atoms.update(eval(Quantity[1], globals(), locals()).atoms(Symbol))
            elif ('sympy' in type(Quantity[1]).__name__) :
                Atoms.update(Quantity[1].atoms(Symbol))
        elif (isinstance(Quantity, basestring)) :
            if (' ' not in Quantity) :
                Atoms.update(eval(Quantity, globals(), locals()).atoms(Symbol))
            else :
                pass # Quantity will just get printed as is
        else :
            pass # Quantity will just get printed as is
    for Atom in Atoms :
        if (Atom in BasicSubstitutions and Atom not in VariableConstants) :
            ReturnList += [ccode(N(BasicSubstitutions[Atom]), assign_to=str(Atom))]
    for Quantity in Quantities :
        if (isinstance(Quantity, list) and len(Quantity)==2) :
            if (isinstance(Quantity[1], basestring)) :
                ReturnList += [str(Quantity[0])+ccode(N(horner(eval(Quantity[1], globals(), locals()))))]
            elif ('sympy' in type(Quantity[1]).__name__) :
                ReturnList += [str(Quantity[0])+ccode(N(horner(Quantity[1])))]
        elif (isinstance(Quantity, basestring)) :
            if (' ' not in Quantity) :
                ReturnList += [ccode(N(horner(eval(Quantity, globals(), locals()))),
                                     assign_to='const double '+Quantity)]
            else :
                ReturnList += [str(Quantity)]
        else :
            ReturnList += [str(Quantity)]
    return ('\n'+' '*Indent).join(ReturnList)

