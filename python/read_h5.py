#!/usr/bin/python
# -*- coding: utf-8 -*-
import h5py
import os
import sys
from colored_io import prtinfo, prtwarn, read_var, die
if (sys.version[0:3] != "2.7"):
    raw_input = input
    not_py27 = True
else:
    not_py27 = False

# For use with PIERNIK.
# Script recognizes integers, floats, logical and string variables.
# String parameter values cannot contain spaces.
# On input it receives a list of names of variables to be read, while it returns a list of parameter values.
# Having array of names and array of values one can simply `exec ("%s=%s" %(name,value))

# searches for variable name, if found splits and appends the value ----------


def append_split_var(line, variable_name, var_array_to_append, param_found):
    if (not_py27):
        words = (str(line).replace(" ", "").replace("\t", "").replace("b'", "")).split('=')  # remove leading b' (byte type) in python3
    else:
        words = str(line.replace(" ", "").replace("\t", "")).split('=')
    if words[0] == variable_name:
        var_value = words[1].split("!")[0]
        if (not_py27):
            var_value = var_value.strip("'")  # remove remaining b' (byte type) in python3
        var_value = determine_type_append(var_value)
        var_array_to_append.append(var_value)
        param_found = True

    return var_array_to_append, param_found
# reads problem.par file embedded in PIERNIK h5 file --------------------------


def read_par(hdf5_filename, var_nam, default_values):  # , var_array):
    if len(hdf5_filename) <= 1:
        sys.exit("Exiting: no filename provided")
    if len(str(var_nam[0])) <= 1 and len(var_nam) <= 1:
        sys.exit("Exiting: no variables to read provided")

    found_parameter = [False] * len(var_nam)
    var_array = []
    value = 0.
    for i in range(len(var_nam)):
        h5File = h5py.File(hdf5_filename, 'r')
        parfile = h5File['problem.par']
        for line in parfile:
            if found_parameter[i] == False:
                value, found_parameter[i] = append_split_var(line, var_nam[i], var_array, found_parameter[i])
    for i in range(len(var_nam)):
        if found_parameter[i] == False:
            prtwarn("Warning: some parameters were not included in problem.par, i.e: %s (default value: %s) . Please provide it:" % (var_nam[i], str(default_values[i])))
            value = input_names_array()
            if (len(value) > 1):
                var_array.append(value)
            else:
                var_array.append(default_values[i])
    return var_array
# for given string value it determines the type of value and returns it -------
# Value types supported: integer, float, boolean and string.
# No lists or arrays so far.


def determine_type_append(var):
    try:
        if ("." in var and ("e" in var or "E" in var)) or ("." in var):  # "." and "e" or just "." means that var is float or fortran boolean
            raise ValueError
        else:
            int(var)
            var = int(var)
            return var
    except ValueError:
        try:
            float(var)
            var = float(var)
            return var
        except ValueError:
            try:
                var_tmp = var.split(' ')
                var_tmp = var_tmp[1]
                if var_tmp == ".true." or var_tmp == ".TRUE.":
                    var = True
                    return var
                elif var_tmp == ".false." or var_tmp == ".FALSE.":
                    var = False
                    return var
                else:
                    var = var  # assume string variable
                    return var
            except:
                pass
        except:
            pass
    except:
        print("Type for provided variable %s not recognized - ")
        return var
# read names if nothing provided


def input_names_array():
    var_names = read_var("Give variables to be read from problem.par file (separated by space): ")
    var_names = var_names.split(' ')
    return var_names
