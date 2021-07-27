
##########################################################################################
# 
# psim_parse.py
# 
#  - Utility functions for parsing text-based inputs to ProBioSim.
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

import os
import re
import sys
from psim_lib import error, internalError, quote, isFloat, getFloat, isBool, getBool, isNumber, isInt, getInt, isIntOrNone, getIntOrNone

##########################################################################################

#
# Functions that define valid names for species and rates
#

# Is given token a valid species name?
def isValidSpeciesName(x):
    if not isinstance(x, str):
        error('Argument to '+quote('isValidSpecies')+' should be a string - found: ' + os.linesep + quote(str(x)))
    speciesNameRegex = r'\A[A-Za-z][A-Za-z0-9\_]*\Z'        
    return re.match(speciesNameRegex, x) is not None

# Is given token a valid rate name?
def isValidRateName(x):
    if not isinstance(x, str):
        error('Argument to '+quote('isValidRateName')+' should be a string - found: ' + os.linesep + quote(str(x)))
    rateNameRegex = r'\A[A-Za-z][A-Za-z0-9\_]*\Z'
    return re.match(rateNameRegex, x) is not None

##########################################################################################

#
# Some helper functions for model parsing
#

# Parse a list (one or two) chemical reactions from a string
def getReactionsListFromString(x):
    if not isinstance(x, str):
        error('Argument to '+quote('getReactionsFromString')+' should be a string - found: ' + os.linesep + quote(str(x)))
    irrevArrowRegex1 = r'[^<]->\{([^\}]+)\}' ## Won't match if no reactants (i.e., irrev arrow is the very first thing in the string - no character before the - to match)
    irrevArrowRegex2 = r'^->\{([^\}]+)\}'    ## This covers the special case where the irrev arrow is the very first thing in the string - no character before the - to match
    revArrowRegex = r'\{([^\}]+)\}<->\{([^\}]+)\}'
    numIrrevArrows1 = len(re.findall(irrevArrowRegex1, x))
    numIrrevArrows2 = len(re.findall(irrevArrowRegex2, x))
    numRevArrows = len(re.findall(revArrowRegex, x))
    if not ((numIrrevArrows1 == 1 and numIrrevArrows2 == 0 and numRevArrows == 0) or
            (numIrrevArrows1 == 0 and numIrrevArrows2 == 1 and numRevArrows == 0) or
            (numIrrevArrows1 == 0 and numIrrevArrows2 == 0 and numRevArrows == 1)):
        totalNumArrows = numIrrevArrows1 + numIrrevArrows + numRevArrows
        error ('In '+quote('getReactionsFromString')+', should be precisely one reaction arrow (reversible or irreversible), but found ' +
               str(totalNumArrows) + ': ' + os.linesep + quote(x))
    if numRevArrows == 0:
        if not (((numIrrevArrows1 == 1 and numIrrevArrows2 == 0) or
                 (numIrrevArrows1 == 0 and numIrrevArrows2 == 1))):
            internalError()
        parts = re.split(irrevArrowRegex1, x) if numIrrevArrows1 == 1 else re.split(irrevArrowRegex2, x)
        if len(parts) != 3:
            error('In '+quote('getReactionsFromString')+', there was a syntax error in this string: ' + os.linesep + quote(x))
        lhs = parts[0].lstrip().rstrip()
        kStr = parts[1].lstrip().rstrip()
        rhs = parts[2].lstrip().rstrip()
        if lhs == '':
            lhsParts = []
        else:
            lhsParts = [p.lstrip().rstrip() for p in lhs.split('+')]
        for l in lhsParts:
            if not isValidSpeciesName(l):
                error('In '+quote('getReactionsFromString')+', expected a reactant species but found: ' + os.linesep + quote(str(l)))
        if rhs == '':
            rhsParts = []
        else:
            rhsParts = [p.lstrip().rstrip() for p in rhs.split('+')]
        for r in rhsParts:
            if not isValidSpeciesName(r):
                error('In '+quote('getReactionsFromString')+', expected a product species but found: ' + os.linesep + quote(str(r)))
        k = getFloat(kStr) if isFloat(kStr) else kStr
        return [{'reactants':lhsParts, 'rate':k, 'products':rhsParts}]
    elif numRevArrows == 1:
        if not (numIrrevArrows1 == 0 and numIrrevArrows2 == 0):
            internalError()
        parts = re.split(revArrowRegex, x)
        if len(parts) != 4:
            error('In '+quote('getReactionsFromString')+', there was a syntax error in this string: ' + os.linesep + quote(x))
        lhs = parts[0].lstrip().rstrip()
        kBwdStr = parts[1].lstrip().rstrip()
        kFwdStr = parts[2].lstrip().rstrip()
        rhs = parts[3].lstrip().rstrip()
        lhsParts = [p.lstrip().rstrip() for p in lhs.split('+')]
        if not (0 <= len(lhsParts) <= 2):
            error('In '+quote('getReactionsFromString')+', should be 0, 1, or 2 reactants for a reversible reaction - found: ' + os.linesep + quote(str(lhsParts)))
        for l in lhsParts:
            if not isValidSpeciesName(l):
                error('In '+quote('getReactionsFromString')+', expected a reactant species but found: ' + os.linesep + quote(str(l)))
        rhsParts = [p.lstrip().rstrip() for p in rhs.split('+')]
        # # Omitting this check: it shouldn't be required and unnecessarily limits expressiveness
        # #if not (0 <= len(rhsParts) <= 2):
        # #   error('In '+quote('getReactionsFromString')+', should be 0, 1, or 2 products for a reversible reaction - found: ' + os.linesep + quote(str(lhsParts)))
        for r in rhsParts:
            if not isValidSpeciesName(r):
                error('In '+quote('getReactionsFromString')+', expected a product species but found: ' + os.linesep + quote(str(r)))
        kFwd = getFloat(kFwdStr) if isFloat(kFwdStr) else kFwdStr
        kBwd = getFloat(kBwdStr) if isFloat(kBwdStr) else kBwdStr
        if len(lhsParts) == 0 and len(rhsParts) == 0:
            error('In '+quote('getReactionsFromString')+', cannot have both reactants and products empty, but found: ' + os.linesep + quote(x))
        return [{'reactants':lhsParts, 'rate':kFwd, 'products':rhsParts},
                {'reactants':rhsParts, 'rate':kBwd, 'products':lhsParts}]
    else:
        internalError()

# Get perturbation from string
def getPerturbationInfoFromString(x):
    tokens = x.split()
    if (len(tokens) == 5 and isValidSpeciesName(tokens[0]) and
        (tokens[1] in ['=', '+=', '-=']) and isFloat(tokens[2]) and
        tokens[3] == '@' and isFloat(tokens[4])):
        speciesName = tokens[0]
        operator = tokens[1]
        amount = getFloat(tokens[2])
        time = getFloat(tokens[4])
        if operator == '-=':
            operator = '+='
            amount *= -1.0
        relative = operator != '='
        action = {'speciesName':speciesName, 'amount':amount, 'relative':relative}
        return (time, action)
    else:
        error('In '+quote('getPerturbationInfoFromString')+', there was a syntax error in this string: ' + os.linesep + quote(x))

##########################################################################################

#
# Modify model based on parsing text input
#
    
# Update rate definition from string
def updateRateDefinitionsFromString(model, x, strict=False):
    tokens = x.split()
    if (len(tokens) == 3 and isValidRateName(tokens[0]) and
        tokens[1] == '=' and isFloat(tokens[2])):
        rateName = tokens[0]
        rateVal = getFloat(tokens[2])
        model.setRateConstant(rateName, rateVal, strict=strict)
    else:
        error('In '+quote('updateRateDefinitionsFromString')+', there was a syntax error in this string: ' + os.linesep + quote(x))

# Update species initialization value from string
def updateSpeciesInitsFromString(model, x, strict=False):
    tokens = x.split()
    if (len(tokens) == 3 and isValidSpeciesName(tokens[0]) and
        (tokens[1] == '=') and isFloat(tokens[2])):
        speciesName = tokens[0]
        initVal = getFloat(tokens[2])
        model.setSpeciesInit(speciesName, initVal, strict=strict)
    else:
        error('In '+quote('updateSpeciesInitsFromString')+', there was a syntax error in this string: ' + os.linesep + quote(x))

# Update reactions list from string
def updateReactionsFromString(model, x):
    model.addReactions(getReactionsListFromString(x))

# Update simulation settings from string
def updateSimulationSettingsFromString(model, x):
    tokens = x.split()
    if (len(tokens) == 2 and tokens[0] == 'length' and isNumber(tokens[1])):
        model.setSimulationLength(getFloat(tokens[1]))
    elif (len(tokens) == 2 and tokens[0] == 'points' and isInt(tokens[1])):
        model.setSimulationPoints(getInt(tokens[1]))
    elif (len(tokens) == 2 and tokens[0] == 'stiff' and isBool(tokens[1])):
        model.setSimulationStiff(getBool(tokens[1]))
    elif (len(tokens) == 2 and tokens[0] == 'seed' and isIntOrNone(tokens[1])):
        model.setSimulationSeed(getIntOrNone(tokens[1]))
    elif (len(tokens) == 2 and tokens[0] == 'rtol' and isFloat(tokens[1])):
        model.setSimulationRtol(getFloat(tokens[1]))
    elif (len(tokens) == 2 and tokens[0] == 'atol' and isFloat(tokens[1])):
        model.setSimulationAtol(getFloat(tokens[1]))
    elif (len(tokens) == 2 and tokens[0] == 'recordAllIfStochastic' and isBool(tokens[1])):
        model.setSimulationRecordAllIfStochastic(getBool(tokens[1]))
    else:
        error('In '+quote('updateSimulationSettingsFromString')+', there was a syntax error in this string: ' + os.linesep + quote(x))
        
# Update perturbations list from string
def updatePerturbationsFromString(model, x):
    (time, action) = getPerturbationInfoFromString(x)
    model.addPerturbation(time, action)

# Update model from a list of strings
def updateFromStrings(model, xs, strict=False):
    for x in xs:
        commentIdx = x.find('#')
        if commentIdx != -1:
            x = x[:commentIdx]
        x = x.strip()
        if x == '':
            pass
        elif x.startswith('simulation '):
            this_line = x[len('simulation '):]
            updateSimulationSettingsFromString(model, this_line)
        elif x.startswith('init '):
            this_line = x[len('init '):]
            updateSpeciesInitsFromString(model, this_line, strict=strict)
        elif x.startswith('rate '):
            this_line = x[len('rate '):]
            updateRateDefinitionsFromString(model, this_line, strict=strict)
        elif x.startswith('reaction '):
            this_line = x[len('reaction '):]
            updateReactionsFromString(model, this_line)
        elif x.startswith('perturbation '):
            this_line = x[len('perturbation '):]
            updatePerturbationsFromString(model, this_line)
        else:
            error('In '+quote('parseString')+', there was a syntax error in this string: ' + os.linesep + quote(x))

# Use a single string to update a model data structure
def updateFromString(model, x, strict=False):
    updateFromStrings(model, x.splitlines(), strict=strict)

##########################################################################################
