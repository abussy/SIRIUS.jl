### SOME DOCS

import re
import sys

def is_array(argattr):
    if "dimension" in argattr:
        return True
    else:
        return False

def is_input(argattr):
    if "in" in argattr:
        return True
    else:
        return False

def is_output(argattr):
    if "out" in argattr:
        return True
    else:
        return False

def is_req(argattr):
    if "required" in argattr:
        return True
    else:
        return False

def is_opt(argattr):
    if "optional" in argattr:
        return True
    else:
        return False

def get_type(argtype):
    types_lookup = {"int": "Cint", "double": "Cdouble", "bool": "Bool",
                    "string": "Cchar"}

    for ftype in types_lookup.keys():
        if ftype in argtype:
            return types_lookup[ftype]

def get_function_by_name(fname, content):

    pattern = r'"""\n\s*{}\(.*?end\b'.format(fname)
    matches = re.findall(pattern, content, re.DOTALL)

    if len(matches) != 1:
        sys.exit("Funciton {} not found or not unique.".format(fname))
    
    return matches[0].strip()

def get_function_arguments(fstring):
    arguments = {}

    pattern = r'^(\w+:\ntype:.*?\nattr:.*?$)'
    regex = re.compile(pattern, re.MULTILINE)
    matches = regex.finditer(fstring)
    
    for match in matches:
        blk = match.group(1)
        argname = blk.splitlines()[0].strip(":")
        argtype = blk.splitlines()[1].split()[1]
        argattr = blk.splitlines()[2].strip("attr: ")

        if "error_code" in argname:
            continue

        arguments[argname] = {"type": argtype, "attr": argattr}

    return arguments

def get_function_call(fname, args, prefix):

    call = "\n"+prefix+"error_code__ = Ref{Cint}(0)\n"

    call += prefix+"LibSirius."+fname+"("
    for argname, argspecs in args.items():
        if "handler" in argspecs["type"]:
            call += argspecs["type"]+".handler_ptr, "
        elif is_array(argspecs["attr"]):
            call += argname+", "
        else:
            call += argname+"__, "

    call += "error_code__)\n"

    call += prefix+'if error_code__[] != 0\n'
    call += prefix+'   error("SIRIUS.{} failed with error code", error_code__[])\n'.format(fname.replace("sirius_", ""))
    call += prefix+'end\n\n' 

    return call

def get_function_output(args, prefix):

    out = ""

    for argname, argspecs in args.items():
        if is_output(argspecs["attr"]) and not is_array(argspecs["attr"]):
            out += "{}{} = {}__[]\n".format(prefix, argname, argname)

    return out

def get_function_signature(fname, args):

    sign = "function {}(".format(fname.replace("sirius_", ""))

    is_mutable = False
    required_so_far = True
    for argname, argspecs in args.items():

        if is_array(argspecs["attr"]) and "get" in fname:
            is_mutable = True

        if is_output(argspecs["attr"]):
            is_mutable = True

        delim = ", "
        if required_so_far and is_opt(argspecs["attr"]):
            delim = "; "
            required_so_far = False

        default = ""
        if is_opt(argspecs["attr"]):
            default = "=nothing"

        name = argname
        if "handler" in argspecs["type"]:
            name = argspecs["type"]

        sign += "{}{}{}".format(delim, name, default)
        
    if is_mutable:
        sign = sign.replace("(, ", "!(")
    else:
        sign = sign.replace("(, ", "(")

    return sign+")\n\n"


def prep_input_arg(argname, argspecs, prefix):

    arg_str = ""

    if is_output(argspecs["attr"]):
        return arg_str

    if is_array(argspecs["attr"]):
        return arg_str

    if "handler" in argspecs["type"]:
        return arg_str

    if is_req(argspecs["attr"]):
        arg_str += prefix+ argname+"__ = Ref{"+get_type(argspecs["type"])+"}("+argname+")\n"

    if is_opt(argspecs["attr"]):
    
        arg_str += prefix + "if isnothing({})\n".format(argname)
        arg_str += prefix + "   "+argname+"__ = Ptr{"+get_type(argspecs["type"])+"}(C_NULL)\n"
        arg_str += prefix + "else\n"
        arg_str += prefix + "   "+argname+"__ = Ref{"+get_type(argspecs["type"])+"}("+argname+")\n"
        arg_str += prefix + "end\n\n"

    return arg_str

def prep_output_arg(argname, argspecs, prefix):

    arg_str = ""

    if is_input(argspecs["attr"]):
        return arg_str

    if is_array(argspecs["attr"]):
        return arg_str

    arg_str += prefix+ argname+"__ = Ref{"+get_type(argspecs["type"])+"}(0)\n"

    return arg_str

def parse_to_julia(fname, args):

    idt = "   "

    # function signature
    func = get_function_signature(fname, args)
    

    # analysis of arguments
    func += idt+"#input arguments (non-array)\n"
    for argname, argspecs in args.items():
        func += prep_input_arg(argname, argspecs, idt)

    func += "\n"+idt+"#output arguments (non-array)\n"
    for argname, argspecs in args.items():
        func += prep_output_arg(argname, argspecs, idt)

    # function call and error handling
    func += get_function_call(fname, args, idt)

    # function output (simple types, no alloc)
    func += get_function_output(args, idt)

    func += "end\n\n"

    return func

with open("functions_to_parse.txt", "r") as myfile:
    tmp = myfile.read()
    functions_to_parse = tmp.splitlines()

with open("LibSirius.jl", "r") as myfile:
    content = myfile.read()

formatted_api = ""
for fname in functions_to_parse:
    body = get_function_by_name(fname, content)
    args = get_function_arguments(body)

    parsed_function = parse_to_julia(fname, args)
    formatted_api += parsed_function

with open("module_skeleton.jl", "r") as myfile:
    skeleton = myfile.read()

with open("SIRIU.jl", "w") as myfile:
    myfile.write(skeleton.replace("__insert_generated_code_here__", formatted_api))

