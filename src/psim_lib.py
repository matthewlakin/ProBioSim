
##########################################################################################
# 
# psim_lib.py
# 
#  - General utility functions for ProBioSim
# 
##########################################################################################
# 
# ProBioSim
# Copyright (C) 2021 Matthew Lakin
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
##########################################################################################

import math
import inspect
import sys
import subprocess

##########################################################################################

# Which version of ProBioSim is this?
__version__ = '0.1'

##########################################################################################

# Display a regular message
def message(x):
    print('... ' + x)

# Display a warning message
def warning(x):
    print('!!! WARNING: ' + x)

# A custom exception class for ProBioSim errors
class ProBioSimError(Exception):
    def __init__(self, message):
        super().__init__(message)
    
# Display an error message and quit
def error(x):
    print('*** ERROR: ' + x)
    raise ProBioSimError(x)

# Signal an internal error and quit
def internalError():
    file = inspect.getfile(inspect.currentframe().f_back)
    lineno = inspect.currentframe().f_back.f_lineno
    print('*** INTERNAL ERROR in file ' + str(file) + ', line ' + str(lineno))
    sys.exit(1)

# Try to execute a system call
def trySystemCall(cmd):
    if subprocess.call(cmd, shell=True) != 0:
        error('Failed to execute command "' + cmd + '"')
def trySystemCall_Quiet(cmd):
    if subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL) != 0:
        error('Failed to execute command "' + cmd + '"')    
def trySystemCall_Silent(cmd):
    if subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        error('Failed to execute command "' + cmd + '"')

##########################################################################################

# Print a float to specified number of decimal places
def fstr(f, dp=1):
    if not (type(dp) is int):
        error('Decimal places (dp) kward to fstr must be a positive integer')
    if not (dp >= 1):
        error('Decimal places (dp) kward to fstr must be a positive integer')
    formatStr = '{0:.'+str(dp)+'f}'
    return formatStr.format(f)

# Quote a string
def quote(x):
    return '"'+x+'"'

##########################################################################################

# See whether token is a float
def isFloat(x):
    try:
        a = float(x)
        return True
    except ValueError:
        return False

# See whether token is an int (or a float that represents an actual integer)
def isInt(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

# See whether a token is an int OR a float
def isNumber(token):
    return isFloat(token) or isInt(token)

# Get a float from a token, if possible
def getFloat(token):
    if not isFloat(token):
        error('Argument to getFloat is not a valid representation of a float ('+str(token)+')')
    return float(token)

# Does a token represent None?
def isNone(token):
    return token == 'None'

# Is a token either float or None?
def isFloatOrNone(token):
    return isFloat(token) or isNone(token)

# Get a float from a token, or None
def getFloatOrNone(token):
    if not (isFloatOrNone(token)):
        error('Argument to getFloatOrNone is not a valid representation of a float or None ('+str(token)+')')
    if isNone(token):
        return None
    else:
        return float(token)

# Get an int from a token that represents an int
def getInt(token):
    if not (isInt(token)):
        error('Argument to getInt is not a valid representation of an int ('+str(token)+')')
    return int(token)

# Is a token either int or None?
def isIntOrNone(token):
    return isInt(token) or isNone(token)

# Get an int from a token, or None
def getIntOrNone(token):
    if not (isIntOrNone(token)):
        error('Argument to getIntOrNone is not a valid representation of an int or None ('+str(token)+')')
    if isNone(token):
        return None
    else:
        return int(token)

# Check whether a token is a boolean
def isBool(token):
    return token in ['true','false','True','False']

# Get a boolean from a token, if possible
def getBool(token):
    if not (isBool(token)):
        error('Argument to getBool is not a valid representation of a bool ('+str(token)+')')
    if token in ['true','True']:
        return True
    elif token in ['false','False']:
        return False
    else:
        internalError()

##########################################################################################

# Return list of pairs of pairwise neighbors in a list
def pairwise(xs):
    return [(xs[idx],xs[idx+1]) for idx in range(len(xs)-1)]

# Return list of pairs of pairwise neighbor INDEXES in a list
def pairwise_indexes(xs):
    return [(idx,idx+1) for idx in range(len(xs)-1)]

# Are the list elements all distinct?
def distinct(xs):
    temp = []
    for x in xs:
        if x in temp:
            return False
        else:
            temp.append(x)
    return True

# Is one float less than or equal to the other, allowing for floating-point imprecision?
def float_leq(f1, f2):
    return (f1 < f2) or math.isclose(f1, f2)

##########################################################################################
