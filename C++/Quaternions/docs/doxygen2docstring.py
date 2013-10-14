#!/usr/bin/env python

"""
Convert doxygen XML to docstring.

The class defined in this file may be used externally.  However, when
called as a script, the docstrings are output in SWIG format.

Usage:

  doxygen2docstring.py ../xml/index.xml GWFrames_Doc.i

index.xml is your doxygen-generated XML file which comes from running
doxygen on the C++ code, and GWFrames_Doc.i is where the output will
be written (the file will be clobbered).

Also, doxygen seems to be failing for my particular versions by
occasionally including the return argument in the docstrings for
certain template types.  The following (regex) replacement in
GWFrames_Doc.i may also be necessary:

  %feature("docstring").*>\**  ->  %feature("docstring")

There are a few specializations to the GWFrames module in this file.
For example, deletion of 'GWFrames::' in certain strings.  Otherwise,
this should be fairly general.

"""


# This code is adapted from doxy2swig.py, which had the following
# license information:
## This code is implemented using Mark Pilgrim's code as a guideline:
##   http://www.faqs.org/docs/diveintopython/kgp_divein.html
##
## Author: Prabhu Ramachandran
## License: BSD style



# Construct dictionary containing function names and docstrings and
# print out docstring features for SWIG.
#
# Each individual version of the function should have the following
# format, so that overloaded functions may be fully documented.
# Different versions of the function should be separated by a couple
# newlines.  If parameters are not documented, that's okay; they
# should still be listed.  Empty parameter lists should just have a
# parenthetical statement like "(no parameters)".  If no description
# is given, that section should be omitted.
#
## Short description as imperative statement.
## ==========================================
##   Parameters
##   ----------
##     const Thingie& Param1
##       Text describing this parameter in some level of detail, which
##       I won't go into here, but should be sufficient to understand
##       this parameter and its function.
##     const Mabob Param2 = Doohicky
##       More info, about Param2 that goes on and on, giving us some
##       useful information about it.
##   Returns
##   -------
##     Doodad& ReturnValue
##   Description
##   -----------
##     Some lengthy description of the workings of this function that
##     provides significant detail.
##     
##     Note that multiple paragraphs are possible -- and likely.  The
##     indentation should be here because this is the only section
##     with line breaks, so the indentation should help the reader
##     scan for function overloads.
# 
# Finally, when all the docstrings are generated for
#   * classes
#   * namespaces
#   * class member functions, free functions in files
# they should be output as `%feature("docstring")` statements



def word_wrap(strings, width=80, ind1=0, ind2=0, prefix=''):
    """ word wrapping function.
        string: the string to wrap
        width: the column number to wrap at
        prefix: prefix each line with this string (goes before any indentation)
        ind1: number of characters to indent the first line
        ind2: number of characters to indent the rest of the lines
        
    """
    stringlist = re.split('\n\s*\n', strings)
    returnstring = ''
    for i,string in enumerate(stringlist):
        string = prefix + ind1 * " " + string.strip()
        newstring = ""
        while len(string) > width:
            # find position of nearest whitespace char to the left of "width"
            marker = width - 1
            while not string[marker].isspace():
                marker = marker - 1
            
            # remove line from original string and add it to the new string
            newline = string[0:marker] + "\n"
            newstring = newstring + newline
            string = prefix + ind2 * " " + string[marker + 1:]
        
        returnstring += newstring + string
        if i<len(stringlist)-1: returnstring += '\n' + prefix + ind1*' ' + '\n'
    
    return returnstring


from xml.dom import minidom
import re
import textwrap
import sys
import types
import os.path
import warnings


class doxy2docstring :
    """
    Convert doxygen XML output to docstring.
    
    Input is the index.xml file doxygen creates.  The file is then
    read, recursively, and the various namespaces, classes, and
    functions are documented.  The result is held in a dictionary,
    where each key is the name of the object and the value is its
    docstring.
    
    This can be output to SWIG format using the class method
    `output_swig("OutputFileName.i")`.
    
    """
    
    def __init__(self, src) :
        """
        Initialize an instance given a source object (file or
        filename).
        
        """
        f = self.open_read(src)
        self.my_dir = os.path.dirname(f.name)
        self.xmldoc = minidom.parse(f).documentElement
        f.close()
        self.documented_objects = {}
        self.tags_re = re.compile(r'<[^>]*>')
        
        # i = 0
        # for memberdef in self.xmldoc.getElementsByTagName('name'):
        #     #print getText(memberdef.childNodes)
        #     self.documented_objects[self.getText(memberdef.childNodes)] = i
        #     i = i+1
        comps = self.xmldoc.getElementsByTagName('compound')
        for c in comps :
            refid = c.attributes['refid'].value
            fname = refid + '.xml'
            if not os.path.exists(fname) :
                fname = os.path.join(self.my_dir,  fname)
            #print("parsing file: {0}".format(fname))
            ## Read file into variable cdoc
            f = self.open_read(fname)
            cdoc = minidom.parse(f).documentElement
            f.close()
            ## Look for classes and namespaces (<compounddef>[kind=="class" || kind=="namespace"])
            for compounddef in cdoc.getElementsByTagName('compounddef') :
                kind = compounddef.attributes['kind'].value
                if kind == 'class' :
                    documented_object = self.human_readable(compounddef.getElementsByTagName('compoundname')[0])
                    docstring = 'class {0}\n{1}\n'.format(documented_object, "="*(6+len(documented_object)))
                    #print(docstring)
                    ## Document class description
                    for child in compounddef.childNodes : # only loop over direct children
                        if ((child.nodeName == 'briefdescription') or (child.nodeName == 'detaileddescription')) :
                            string = self.human_readable(child)
                            if len(string)>0 :
                                docstring += word_wrap(string, ind1=2, ind2=2) + '\n  \n'
                    ## Document variables (children with <memberdef>[kind=="variable"])
                    MemberVariables = False
                    memberdefs = compounddef.getElementsByTagName('memberdef')
                    for memberdef in memberdefs :
                        if (memberdef.attributes['kind'].value=='variable') :
                            if (MemberVariables != True) :
                                docstring += '  Member variables\n  ' + '-'*16 + '\n'
                                MemberVariables = True
                            declaration = self.human_readable(memberdef.getElementsByTagName('type')[0]).replace('GWFrames::','').replace('std::','') + ' ' \
                                + self.human_readable(memberdef.getElementsByTagName('name')[0]) + '\n'
                            docstring += '    ' + declaration
                            ## Add description
                            for child in memberdef.childNodes :
                                if ((child.nodeName == 'briefdescription') or (child.nodeName == 'detaileddescription')) :
                                    string = self.human_readable(child)
                                    if len(string)>0 :
                                        docstring += '\n      ' + word_wrap(string, ind1=6, ind2=6)
                    if(MemberVariables == True) : docstring += '  \n'
                    ## Document non-public functions (children with <memberdef>[kind=="function" && prot!="public"])
                    MemberFunctions = False
                    memberdefs = compounddef.getElementsByTagName('memberdef')
                    for memberdef in memberdefs :
                        if ((memberdef.attributes['kind'].value=='function') and (memberdef.attributes['prot'].value!='public')) :
                            if (MemberFunctions != True) :
                                docstring += '  Non-public member functions\n  ' + '-'*27 + '\n'
                                MemberFunctions = True
                            declaration = self.human_readable(memberdef.getElementsByTagName('type')[0]).replace('GWFrames::','').replace('std::','') + ' ' \
                                + self.human_readable(memberdef.getElementsByTagName('name')[0]) + '\n'
                            docstring += '    ' + declaration
                            for child in memberdef.childNodes :
                                if ((child.nodeName == 'briefdescription') or (child.nodeName == 'detaileddescription')) :
                                    string = self.human_readable(child)
                                    if len(string)>0 :
                                        docstring += word_wrap(string, ind1=6, ind2=6)
                    if(MemberFunctions == True) : docstring += '  \n'
                    ## Add to self.documented_objects
                    self.documented_objects[documented_object] = docstring
                elif kind == 'namespace' :
                    documented_object = self.human_readable(compounddef.getElementsByTagName('compoundname')[0])
                    docstring = 'namespace {0}\n{1}\n'.format(documented_object, "="*(10+len(documented_object)))
                    #print(docstring)
                    ## Document namespace description
                    for child in compounddef.childNodes : # only loop over direct children
                        if ((child.nodeName == 'briefdescription') or (child.nodeName == 'detaileddescription')) :
                            string = self.human_readable(child)
                            if len(string)>0 :
                                docstring += word_wrap(string, ind1=2, ind2=2) + '\n  \n'
                    ## Add to self.documented_objects
                    self.documented_objects[documented_object] = docstring
                    
            ## Look for <memberdef>[kind="function" prot="public"]
            for memberdef in cdoc.getElementsByTagName('memberdef') :
                if ((memberdef.attributes['kind'].value=='function') and (memberdef.attributes['prot'].value=='public')) :
                    definition = self.human_readable(memberdef.getElementsByTagName('definition')[0])
                    typestring = self.human_readable(memberdef.getElementsByTagName('type')[0])
                    documented_object = definition.replace(typestring+' ', '')
                    ## Make header with brief description
                    brief = self.human_readable(memberdef.getElementsByTagName('briefdescription')[0])
                    #print(documented_object)
                    docstring = ''
                    if documented_object in self.documented_objects :
                        docstring = self.documented_objects[documented_object] + '\n'
                    ## Document parameters
                    docstring += brief + '\n' + len(brief)*'=' + '\n  Parameters\n  ----------\n'
                    params = memberdef.getElementsByTagName('param')
                    if len(params)==0 :
                        docstring += '    (none)\n  \n'
                    else :
                        for param in params :
                            type = self.human_readable(param.getElementsByTagName('type')[0])
                            try:
                                declname = self.human_readable(param.getElementsByTagName('declname')[0])
                            except IndexError:
                                declname = ''
                            docstring += '    ' + (type + ' ' + declname).replace('GWFrames::','').replace('std::','')
                            default = param.getElementsByTagName('defval')
                            if len(default)>0 : docstring += ' = ' + self.human_readable(default[0]).replace('GWFrames::','').replace('std::','')
                            docstring += '\n'
                            for parametername in memberdef.getElementsByTagName('parametername') :
                                if self.getText(parametername.childNodes)==declname :
                                    #print(parametername.parentNode.nextSibling.nextSibling.toprettyxml())
                                    parameterdescription = self.human_readable(parametername.parentNode.nextSibling.nextSibling).strip()
                                    if len(parameterdescription)>0 : docstring += word_wrap(parameterdescription, ind1=6, ind2=6) + '\n'
                        docstring += '  \n'
                    ## Document return quantity
                    docstring += '  Returns\n  -------\n    '
                    returnType = self.human_readable(memberdef.getElementsByTagName('type')[0]).replace('GWFrames::','').replace('std::','')
                    if len(returnType.strip())==0 :
                        returnType = self.human_readable(memberdef.getElementsByTagName('name')[0])
                    docstring += returnType + '\n  \n'
                    ## Document other desciptions
                    Descriptions=False
                    for child in memberdef.childNodes :
                        if ((child.nodeName == 'detaileddescription') or (child.nodeName == 'inbodydescription')) :
                            parameterlist = child.getElementsByTagName('parameterlist')
                            if len(parameterlist)>0 : parameterlist[0].parentNode.removeChild(parameterlist[0])
                            string = self.human_readable(child)
                            if len(string)>0 :
                                if (Descriptions == False) :
                                    docstring += '  Description\n  -----------\n'
                                    Descriptions = True
                                docstring += word_wrap(string, ind1=4, ind2=4) + '\n  \n'
                    #if(Descriptions == True) : docstring += '  \n'
                    ##   Add to self.documented_objects
                    self.documented_objects[documented_object] = docstring
    
    def open_read(self, source):
        if hasattr(source, "read"):
            return source
        else:
            return open(source)
    
    def open_write(self, dest):
        if hasattr(dest, "write"):
            return dest
        else:
            return open(dest, 'w')
    
    def human_readable(self, node) :
        return re.sub(self.tags_re, '',
                      node.toxml()
                      .replace('</para>','\n\n')) \
                      .replace('&lt; ','<') \
                      .replace('&lt;','<') \
                      .replace(' &gt;','>') \
                      .replace('&gt;','>') \
                      .replace(' &amp;','&') \
                      .replace('&amp;','&') \
                      .replace('&quot;',"'") \
                      .replace('&apos;',"'") \
                      .replace('\\','\\\\') \
                      .strip()
    
    def getText(self, nodelist):
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(node.data)
        return ''.join(rc)
    
    def output_swig(self, OutputFileName):
        o = self.open_write(OutputFileName)
        for key, value in self.documented_objects.iteritems():
            o.write('%feature("docstring") {0} """\n{1}"""\n\n'.format(key, value))
        o.close()




def main(input, output):
    p = doxy2docstring(input)
    p.output_swig(output)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print __doc__
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
