#!/usr/bin/env python

### This python script parses the raw julia wrappers around the SIRIUS library, and spits out
### formatted Julia wrappers that are easier to use
### This deals with the Fortran style of the SIRIUS API, in which most variables are passed as pointers:
### scalar input pointers are initialized as Ref{Ctype}(input), while arrays are directly fed to the C code
### The script also adds error handling for the API calls
###
### The script expects 3 files to exist:
### 1) LibSirius.jl: the raw julia wrappers generated from the C headers by Clang.jl
### 2) functions_to_parse.txt: the list of functions of interest to be parsed by this script.
###                            The SIRIUS API is full of QE specific routines not needed here.
###                            Moreover, some functions/structs are wrapped by hand, as we add
###                            additional information on the julia side (e.g. handlers)
### 3) module_skeleton.jl: the skeleton code for the SIRIUS.jl module, which includes the hand
###                        written wrappers for the handlers and a place holder for the code generated
###                        by this script

import re
import sys

### Collection of functions to check the properties of a given argument based on the API documentation
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

def is_handle(argtype):
    if "handle" in argtype:
        return True
    else:
        return False

def is_string(argtype):
    if "string" in argtype:
        return True
    else:
        return False

def get_type(argtype):
    types_lookup = {"int": "Cint", "double": "Cdouble", "bool": "Bool",
                    "string": "Cchar"}

    for ftype in types_lookup.keys():
        if ftype in argtype:
            return types_lookup[ftype]

### Uses regex to extract the code block of a given function in a file
def get_function_by_name(fname, content):

    pattern = r'"""\n\s*{}\(.*?end\b'.format(fname)
    matches = re.findall(pattern, content, re.DOTALL)

    if len(matches) != 1:
        sys.exit("Funciton {} not found or not unique.".format(fname))
    
    return matches[0].strip()

### Uses regex to extract the argument documentation from a function code block
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

### Generates a function call based on a name and argument properties
def gen_function_call(fname, args, prefix):

    call = "\n"+prefix+"error_code__ = Ref{Cint}(0)\n"

    call += prefix+"redirect_stdio(;stdout=get_outstream()) do\n"

    call += 2*prefix+"LibSirius."+fname+"("
    for argname, argspecs in args.items():
        if is_handle(argspecs["type"]):
            call += argspecs["type"]+".handler_ptr, "
        elif is_array(argspecs["attr"]):
            call += argname+", "
        elif is_string(argspecs["type"]):
            call += argname+", "
        else:
            call += argname+"__, "

    call += "error_code__)\n"

    call += 2*prefix+"Base.Libc.flush_cstdio()\n"
    call += prefix+"end\n"

    call += prefix+'if error_code__[] != 0\n'
    call += prefix+'   error("SIRIUS.{} failed with error code", error_code__[])\n'.format(fname.replace("sirius_", ""))
    call += prefix+'end\n\n' 

    return call

### Generates the return statement of a function based on the argument documentation, for scalar variables
def gen_function_output(args, prefix):

    out = ""

    for argname, argspecs in args.items():
        if is_output(argspecs["attr"]) and not is_array(argspecs["attr"]):
            out += argname+"__[], "

    if len(out) > 0:
        return prefix + "return " + out[:-2] + "\n"
    else:
        return ""

### Generates a function signature based on a function name and argument documentation
### Pass Fortran optional arguments as keword arguments defaulting to nothing
def gen_function_signature(fname, args):

    sign = "function {}(".format(fname.replace("sirius_", ""))

    is_mutable = False
    required_so_far = True
    for argname, argspecs in args.items():

        if is_output(argspecs["attr"]) and not is_array(argspecs["attr"]):
            continue

        if is_array(argspecs["attr"]) and "get" in fname:
            is_mutable = True

        delim = ", "
        if required_so_far and is_opt(argspecs["attr"]):
            delim = "; "
            required_so_far = False

        default = ""
        if is_opt(argspecs["attr"]):
            if is_array(argspecs["attr"]):
                default = "=C_NULL"
            else:
                default = "=nothing"

        name = argname
        if is_handle(argspecs["type"]):
            name = argspecs["type"]

        sign += "{}{}{}".format(delim, name, default)
        
    if is_mutable:
        sign = sign.replace("(, ", "!(")
    else:
        sign = sign.replace("(, ", "(")

    return sign+")\n\n"

# Prepares input scalar arguments for API call: initializes Ref{Ctype}(input_arg) and deals
# with Fortran optional arguments by passing C_NULL pointers by default
def prep_input_arg(argname, argspecs, prefix):

    arg_str = ""

    if is_output(argspecs["attr"]):
        return arg_str

    if is_array(argspecs["attr"]):
        return arg_str

    if is_handle(argspecs["type"]):
        return arg_str

    if is_string(argspecs["type"]):
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

# Prepares output scalar arguments for API call: create a Ref{Ctype} to be converted later
def prep_output_arg(argname, argspecs, prefix):

    arg_str = ""

    if is_input(argspecs["attr"]):
        return arg_str

    if is_array(argspecs["attr"]):
        return arg_str

    arg_str += prefix+ argname+"__ = Ref{"+get_type(argspecs["type"])+"}(0)\n"

    return arg_str

# Creates a Julia code block for the given function and arguments
def parse_to_julia(fname, args):

    idt = "   "

    # function signature
    func = gen_function_signature(fname, args)
    

    # analysis of arguments
    func += idt+"#input arguments (non-array)\n"
    for argname, argspecs in args.items():
        func += prep_input_arg(argname, argspecs, idt)

    func += "\n"+idt+"#output arguments (non-array)\n"
    for argname, argspecs in args.items():
        func += prep_output_arg(argname, argspecs, idt)

    # function call and error handling
    func += gen_function_call(fname, args, idt)

    # function output (simple types, no alloc)
    func += gen_function_output(args, idt)

    func += "end\n\n"

    return func

def main():
    ### Parsing all required functions
    print("Generating second layer of wrappers with Python...")
    
    with open("functions_to_parse.txt", "r") as myfile:
        tmp = myfile.read()
        functions_to_parse = tmp.splitlines()
    
    with open("../../src/LibSirius.jl", "r") as myfile:
        content = myfile.read()
    
    formatted_api = ""
    for fname in functions_to_parse:
        body = get_function_by_name(fname, content)
        args = get_function_arguments(body)
    
        parsed_function = parse_to_julia(fname, args)
        formatted_api += parsed_function
    
    with open("module_skeleton.jl", "r") as myfile:
        skeleton = myfile.read()
    
    with open("../../src/SIRIUS.jl", "w") as myfile:
        myfile.write(skeleton.replace("__insert_generated_code_here__", formatted_api))

    print("Done!")

if __name__ == "__main__":
    main()
