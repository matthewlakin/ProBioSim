
##########################################################################################
#
# psim_model.py
#
#  - Class for creating and manipulating ProBioSim models, along with related functions.
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
from psim_lib import error, internalError, warning, message, trySystemCall, quote
import psim_parse

##########################################################################################

#
# Produce positive and negative species names for dual rail.
# This function codifies a notational convenion whereby,
# for a dual rail signal X, the positive and negative species
# are Xp and Xm respectively (for X-plus and X-minus).
# By using this, we ensure we always stick to a consistent naming
# convention for dual-rail species.
#
def positiveDualRailSpecies(x):
    return x+'p'
def negativeDualRailSpecies(x):
    return x+'m'

##########################################################################################

#
# A class for models to be simulated
#

class Model(object):

    ######################################################################################

    #
    # Basic setup and wellformedness checking
    #
    
    # When creating a model, it is empty by default and has default simulation settings
    def __init__(self):
        defaultSimSettings = {'length':1.0, 'points':101, 'stiff':True, 'seed':None, 'rtol':1e-8, 'atol':1e-8, 'recordAllIfStochastic':False}
        self.reactions = []
        self.rateDefinitions = {}
        self.speciesInits = {}
        self.perturbations = []
        self.simSettings = defaultSimSettings
        self.allSpeciesNames = []
        self.checkWellFormedness()

    # Return string representation of model
    def __str__(self):
        return self.getASCII()

    # Do some basic checking of well-formedness
    def checkWellFormedness(self):
        if not (self.simSettings['length'] > 0):
            error('Model not well-formed because simulation length is <= 0')
        if not (self.simSettings['points'] > 0):
            error('Model not well-formed because number of points is <= 0')
        if not (type(self.simSettings['points']) is int):
            error('Model not well-formed because number of points must be an integer')
        if not ((self.simSettings['seed'] is None) or (type(self.simSettings['seed']) is int)):
            error('Model not well-formed because "seed" parameter value must either be None or an integer')
        if not (type(self.simSettings['stiff']) is bool):
            error('Model not well-formed because "stiff" parameter value must be a bool')
        if not (type(self.simSettings['rtol']) is float and self.simSettings['rtol'] > 0):
            error('Model not well-formed because "rtol" parameter value must be a float')
        if not (type(self.simSettings['atol']) is float and self.simSettings['atol'] > 0):
            error('Model not well-formed because "atol" parameter value must be a float')
        if not (type(self.simSettings['recordAllIfStochastic']) is bool):
            error('Model not well-formed because "recordAllIfStochastic" parameter value must be a bool')
        for rateName in self.rateDefinitions:
            if not (self.rateDefinitions[rateName] >= 0):
                error('Model not well-formed because value for rate parameter '+quote(str(rateName))+' is <0')
        for speciesName in self.speciesInits:
            if not (self.speciesInits[speciesName] >= 0):
                error('Model not well-formed because initial value for species '+quote(str(speciesName))+' is <0')
        for speciesName in self.allSpeciesNames:
            if not (psim_parse.isValidSpeciesName(speciesName)):
                error('Model not well-formed because species name '+quote(str(speciesName))+' is not syntactically valid (must be a letter followed by 0+ letters, numbers, or underscores)')
            if not (speciesName in self.speciesInits):
                error('Model not well-formed because species name '+quote(str(speciesName))+' from initial species values is not represented in the list of all species names')
        # Check for well-formedness of individual reactions
        for reaction in self.reactions:
            if not (len(reaction['reactants']) > 0 or len(reaction['products']) > 0):
                error('Model not well-formed because a reaction has no reactants or products')
            if not (sorted(reaction['reactants']) != sorted(reaction['products'])):
                error('Model not well-formed because a reaction has identical reactants and products')
            if type(reaction['rate']) is float:
                if not (reaction['rate'] >= 0):
                    error('Model not well-formed because a reaction has a negative float as its rate constant value: '+str(reaction['rate']))
            elif type(reaction['rate']) is str:
                if not psim_parse.isValidRateName(reaction['rate']):
                    error('Model not well-formed because a reaction has an invalid name for its rate constant (must be a letter followed by 0+ letters, numbers, or underscores)')
                if not (reaction['rate'] in self.rateDefinitions):
                    error('Model not well-formed because a reaction has a rate constant name '+quote(str(reaction['rate']))+' not featured in the rate definitions')
            else:
                error('Model not well-formed because a reaction has a rate entry ('+str(reaction['rate'])+') that is neither a float nor a string')
            if 'annotation' in reaction:
                if not (reaction['annotation'] is None or type(reaction['annotation']) is str):
                    error('Model not well-formed because a reaction has an annotation that is neither None nor a string')
        # Check that there are no duplicated reactions (same reactions and products)
        for (i,ri) in enumerate(self.reactions):
            for (j,rj) in enumerate(self.reactions):
                if i != j:
                    if ((sorted(ri['reactants']) == sorted(rj['reactants'])) and
                        (sorted(ri['products']) == sorted(rj['products']))):
                        error('Model not well-formed because it contains a duplicate reaction (reactants: '+str(sorted(ri['reactants']))+', products: '+str(sorted(ri['products'])))
        # Check for well-formedness of perturbations
        lastTime = 0.0
        for perturbation in self.perturbations:
            for key in ['time', 'actions', 'functions']:
                if key not in perturbation:
                    error('Model not well-formed because a perturbation was missing the key: '+str(quote(key)))
            if not (perturbation['time'] > lastTime):
                error('Model not well-formed because perturbation times are not in ascending order.')
            lastTime = perturbation['time']
            for action in perturbation['actions']:
                if not (type(action) is dict):
                    error('Model not well-formed because a perturbation has an "action" entry that is not a dict')
                for key in ['speciesName', 'amount', 'relative']:
                    if key not in action:
                        error('Model not well-formed because a perturbation action was missing the key: '+str(quote(key)))
                if not (action['speciesName'] in self.allSpeciesNames):
                    error('Model not well-formed because species name '+quote(str(speciesName))+' in perturbation action is not mentioned in the list of all species in the model')
                if not ((type(action['amount']) is float) or (type(action['amount']) is int)):
                    error('Model not well-formed because a perturbation has an action entry whose "amount" is not a float or an int')
                if not (type(action['relative']) is bool):
                    error('Model not well-formed because a perturbation has an action entry whose "relative" value is not a bool')
            for f in perturbation['functions']:
                if not callable(f):
                    error('Model not well-formed because a perturbation has a function entry that is not callable')

    # Run post-update fixes on model to keep it well-formed
    def runPostUpdateFixes(self):
        # Update the list of all species mentioned
        for reaction in self.reactions:
            for x in reaction['reactants'] + reaction['products']:
                if x not in self.allSpeciesNames:
                    self.allSpeciesNames += [x]
        for x in self.speciesInits:
            if x not in self.allSpeciesNames:
                self.allSpeciesNames += [x]
        for p in self.perturbations:
            for a in p['actions']:
                x = a['speciesName']
                if x not in self.allSpeciesNames:
                    self.allSpeciesNames += [x]
        # Set 0 as default rate constant for any named rate constants not explicitly set
        for reaction in self.reactions:
            if type(reaction['rate']) is str:
                if reaction['rate'] not in self.rateDefinitions:
                    self.rateDefinitions[reaction['rate']] = 0.0
        # Fix the species inits by adding 0 as the default initial value for any species not explicitly initialized
        for speciesName in self.allSpeciesNames:
            if speciesName not in self.speciesInits:
                self.speciesInits[speciesName] = 0.0
        # Sort all perturbations in ascending time order
        self.perturbations.sort(key=lambda p:p['time'])
        # Check well-formedness
        self.checkWellFormedness()

    ######################################################################################

    #
    # Modify model programmatically - all species names
    #

    # Get all species names referenced in the model
    def getAllSpeciesNames(self):
        return self.allSpeciesNames
    
    ######################################################################################

    #
    # Modify model programmatically - rate constants
    #

    # Get all rate constant values
    def getRateConstants(self):
        return self.rateDefinitions

    # Get a particular rate constant value
    def getRateConstant(self, rateName):
        if rateName in self.rateDefinitions:
            return self.rateDefinitions[rateName]
        else:
            error('No rate constant value defined for '+quote(str(rateName)))

    # Set rate constant
    def setRateConstant(self, rateName, value, strict=False):
        if not (value >= 0):
            error('Value provided for rate parameter '+quote(str(rateName))+' is <0 '+quote(str(value)))
        if strict and (rateName in self.rateDefinitions):
            error('In "setRateConstant", rate '+rateName+' already defined')
        self.rateDefinitions[rateName] = value
        self.runPostUpdateFixes()

    ######################################################################################

    #
    # Modify model programmatically - species initial values
    #

    # Get all species init values
    def getSpeciesInits(self):
        return self.speciesInits

    # Get a particular species init value
    def getSpeciesInit(self, species):
        if species in self.speciesInits:
            return self.speciesInits[species]
        else:
            error('No initial value defined for '+quote(str(species)))

    # Set species init value
    def setSpeciesInit(self, species, value, strict=False):
        if not (value >= 0):
            error('Initial value provided for species '+quote(str(species))+' is <0 '+quote(str(value)))
        if strict and (species in self.speciesInits):
            error('In "setSpeciesInit", species '+species+' already initialized')
        self.speciesInits[species] = float(value)
        self.runPostUpdateFixes()

    # Get a particular species init value for a dual rail species
    def getSpeciesInitDual(self, signal):
        positiveSpeciesName = positiveDualRailSpecies(signal)
        negativeSpeciesName = negativeDualRailSpecies(signal)
        return getSpeciesInit(positiveSpeciesName) - getSpeciesInit(negativeSpeciesName)

    # Set species init values for a dual rail species
    def setSpeciesInitDual(self, signal, value, strict=False):
        positiveSpeciesName = positiveDualRailSpecies(signal)
        negativeSpeciesName = negativeDualRailSpecies(signal)
        if value < 0.0:
            self.setSpeciesInit(positiveSpeciesName, 0.0, strict=strict)
            self.setSpeciesInit(negativeSpeciesName, -value, strict=strict)
        elif value > 0.0:
            self.setSpeciesInit(positiveSpeciesName, value, strict=strict)
            self.setSpeciesInit(negativeSpeciesName, 0.0, strict=strict)
        else:
            self.setSpeciesInit(positiveSpeciesName, 0.0, strict=strict)
            self.setSpeciesInit(negativeSpeciesName, 0.0, strict=strict)

    ######################################################################################

    #
    # Modify model programmatically - reactions
    #

    # Get reactions
    def getReactions(self):
        return self.reactions
    
    # Add a reaction
    def addReaction(self, reaction):
        self.reactions += [reaction]
        self.runPostUpdateFixes()

    # Add multiple reactions, in a single call
    def addReactions(self, reactions):
        for reaction in reactions:
            self.addReaction(reaction)

    # Add a reaction (pieces given separately as individual arguments, rather than as a dict)
    def addReactionSeparately(self, reactants, rate, products, annotation=None):
        r = {'reactants':reactants, 'rate':rate, 'products':products}
        if annotation is not None:
            r['annotation'] = annotation
        self.addReaction(r)

    ######################################################################################

    #
    # Modify model programmatically - simulation settings
    #

    # Get simulation length setting
    def getSimulationLength(self):
        return self.simSettings['length']

    # Update simulation length setting
    def setSimulationLength(self, length):
        self.simSettings['length'] = length
        self.runPostUpdateFixes()

    # Get simulation points setting
    def getSimulationPoints(self):
        return self.simSettings['points']

    # Update simulation points setting
    def setSimulationPoints(self, points):
        self.simSettings['points'] = points
        self.runPostUpdateFixes()

    # Get (deterministic) simulation stiffness setting
    def getSimulationStiff(self):
        return self.simSettings['stiff']

    # Update (deterministic) simulation stiffness setting
    def setSimulationStiff(self, stiff):
        self.simSettings['stiff'] = stiff
        self.runPostUpdateFixes()

    # Get (stochastic) simulation seed setting
    def getSimulationSeed(self):
        return self.simSettings['seed']

    # Update (stochastic) simulation seed setting
    def setSimulationSeed(self, seed):
        self.simSettings['seed'] = seed
        self.runPostUpdateFixes()

    # Get (deterministic) simulation relative tolerance setting
    def getSimulationRtol(self):
        return self.simSettings['rtol']

    # Update (deterministic) simulation relative tolerance setting
    def setSimulationRtol(self, rtol):
        self.simSettings['rtol'] = rtol
        self.runPostUpdateFixes()
    
    # Get (deterministic) simulation absolute tolerance setting
    def getSimulationAtol(self):
        return self.simSettings['atol']

    # Update (deterministic) simulation absolute tolerance setting
    def setSimulationAtol(self, atol):
        self.simSettings['atol'] = atol
        self.runPostUpdateFixes()

    # Get (stochastic) simulation setting to sample after every reaction
    def getSimulationRecordAllIfStochastic(self):
        return self.simSettings['recordAllIfStochastic']

    # Update (stochastic) simulation setting to sample after every reaction
    def setSimulationRecordAllIfStochastic(self, b):
        self.simSettings['recordAllIfStochastic'] = b
        self.runPostUpdateFixes()

    ######################################################################################

    #
    # Modify model programmatically - perturbations
    #

    # Get perturbations
    def getPerturbations(self):
        return self.perturbations
    
    # Add a perturbation
    def addPerturbation(self, time, action):
        # Get a list of perturbation times (should be sorted)
        currentTimes = []
        for p in self.perturbations:
            if p['time'] not in currentTimes:
                currentTimes += [p['time']]
        # Add new time if necessary, otherwise, add action in appropriate place
        if time not in currentTimes:
            self.perturbations += [{'time':time, 'actions':[action], 'functions':[]}]
        else:
            for p in self.perturbations:
                if p['time'] == time:
                    p['actions'] += [action]
        self.runPostUpdateFixes()

    # Add a perturbation (pieces given separately as individual arguments, rather than as a dict)
    def addPerturbationSeparately(self, time, speciesName, amount, relative):
        self.addPerturbation(time, {'speciesName':speciesName, 'amount':float(amount), 'relative':relative})

    # Add a perturbation to model, with multiple actions
    def addPerturbationMultipleActions(self, time, actions):
        for action in actions:
            self.addPerturbation(time, action)

    # Add a perturbation for a dual rail signal
    def addPerturbationDual(self, time, dual_action):
        positiveSpeciesName = positiveDualRailSpecies(dual_action['signalName'])
        negativeSpeciesName = negativeDualRailSpecies(dual_action['signalName'])
        if dual_action['relative']:
            if dual_action['amount'] < 0.0:
                action = {'speciesName':negativeSpeciesName, 'amount':-1.0*float(dual_action['amount']), 'relative':True}
                self.addPerturbation(time, action)
            elif dual_action['amount'] > 0.0:
                action = {'speciesName':positiveSpeciesName, 'amount':float(dual_action['amount']), 'relative':True}
                self.addPerturbation(time, action)
            else:
                pass
        else:
            if dual_action['amount'] < 0.0:
                pos_action = {'speciesName':positiveSpeciesName, 'amount':0.0,                               'relative':False}
                ned_action = {'speciesName':negativeSpeciesName, 'amount':-1.0*float(dual_action['amount']), 'relative':False}
                self.addPerturbationMultipleActions(time, [pos_action, neg_action])
            elif dual_action['amount'] > 0.0:
                pos_action = {'speciesName':positiveSpeciesName, 'amount':float(dual_action['amount']), 'relative':False}
                ned_action = {'speciesName':negativeSpeciesName, 'amount':0.0,                          'relative':False}
                self.addPerturbationMultipleActions(time, [pos_action, neg_actiona])
            else:
                pos_action = {'speciesName':positiveSpeciesName, 'amount':0.0, 'relative':False}
                ned_action = {'speciesName':negativeSpeciesName, 'amount':0.0, 'relative':False}
                self.addPerturbationMultipleActions(time, [pos_action, neg_actiona])

    # Add a perturbation for a dual rail signal (pieces given separately as individual arguments, rather than as a dict)
    def addPerturbationSeparatelyDual(self, time, signalName, amount, relative):
        self.addPerturbationDual(time, {'signalName':signalName, 'amount':float(amount), 'relative':relative})

    # Add a perturbation for a dual-rail signal, with multiple actions
    def addPerturbationMultipleActionsDual(self, time, dual_actions):
        for dual_action in dual_actions:
            self.addPerturbationDual(time, dual_action)

    # Add a perturbation to the model, containing arbitrary code to run!
    def addPerturbationAsCode(self, time, f):
        # Get a list of perturbation times (should be sorted)
        currentTimes = []
        for p in self.perturbations:
            if p['time'] not in currentTimes:
                currentTimes += [p['time']]
        # Add new time if necessary, otherwise, add action in appropriate place
        if time not in currentTimes:
            self.perturbations += [{'time':time, 'actions':[], 'functions':[f]}]
        else:
            for p in self.perturbations:
                if p['time'] == time:
                    p['functions'] += [f]
        self.runPostUpdateFixes()
    
    # Add a perturbation to the model, containing arbitrary code to run, with multiple functions to fire at one time!
    def addPerturbationAsCodeMultipleFunctions(self, time, fs):
        for f in fs:
            self.addPerturbationAsCode(time, f)
    
    ######################################################################################
            
    # Allow incremental updating of model via parsed string inputs
    def updateFromString(self, x):
        return psim_parse.updateFromString(self, x)
    def updateFromStrings(self, xs):
        return psim_parse.updateFromStrings(self, xs)
    
    ######################################################################################

    # Infer dual rail signal species
    def inferDualRailSignals(self):
        possible_positives = []
        possible_negatives = []
        for species in self.allSpeciesNames:
            if species.endswith('p'):
                possible_positives += [species[:-1]]
            elif species.endswith('m'):
                possible_negatives += [species[:-1]]
        signals = []
        for s in possible_positives:
            if s in possible_negatives:
                signals += [s]
        return signals

    ######################################################################################

    #
    # Helpful functions for model output below
    #

    # Get all reaction annotations used in the model
    def getAllReactionAnnotations(self):
        res = []
        for r in self.reactions:
            if 'annotation' in r and r['annotation'] is not None and r['annotation'] not in res:
                res += [r['annotation']]
        return res
    
    # Sort reactions into classes based on their annotations, and return
    def sortReactionsByAnnotations(self, reaction_annotation_sections):
        if reaction_annotation_sections is None:
            error('Cannot call "sortReactionsByAnnotations" with reaction_annotation_sections=None')
        otherTag = '__Other__'
        odict = {}
        for s in reaction_annotation_sections:
            if s in odict:
                error('Duplicated name '+quote(str(s))+' in reaction_annotation_sections argument to "sortReactionsByAnnotations"')
            odict[s] = []
        if otherTag in odict:
            error('The reserved name for unannotated reactions ('+otherTag+') was used as a reaction annotation!')
        odict[otherTag] = []
        num_reactions_added = 0
        for r in self.reactions:
            if r['annotation'] is None:
                odict[otherTag] += [r]
                num_reactions_added += 1
            else:
                if r['annotation'] in odict:
                    odict[r['annotation']] += [r]
                    num_reactions_added += 1
                else:
                    odict[otherTag] += [r]
                    num_reactions_added += 1
        if not (num_reactions_added == len(self.reactions)):
            internalError()
        res = []
        for s in reaction_annotation_sections:
            res += [(s,odict[s])]
        if otherTag in odict and len(odict[otherTag]) > 0:
            res += [(otherTag, odict[otherTag])]
        return res

    ######################################################################################

    #
    # Parseable ASCII representation
    #

    # Produce an ASCII text representation of the model, that could be re-parsed
    def getASCII(self, title=None, includeZeroInits=False, reaction_annotation_sections=None):
        strings = []
        if title is not None:
            strings += ['#', '#', '# '+title, '#', '#', '']
        strings += ['#', '# Number of reactions: '+str(len(self.reactions)), '# Number of species: '+str(len(self.speciesInits.items())), '#', '']
        if self.reactions != []:
            strings += ['#', '# Reactions', '#', '']
            if reaction_annotation_sections is None:
                for reaction in self.reactions:
                    thisReactantsString = ' + '.join(reaction['reactants'])
                    thisRateString = ' ->{'+str(reaction['rate'])+'} '
                    thisProductsString = ' + '.join(reaction['products'])
                    if 'annotation' in reaction:
                        thisAnnotationString = ' # '+reaction['annotation'] if reaction['annotation'] is not None else ''
                    else:
                        thisAnnotationString = ''
                    strings += ['reaction ' + thisReactantsString + thisRateString + thisProductsString + thisAnnotationString]
            else:
                sortedByAnnotation = self.sortReactionsByAnnotations(reaction_annotation_sections)
                numSections = len(sortedByAnnotation)
                for (idx,(s, sreactions)) in enumerate(sortedByAnnotation):
                    numReactions = len(sreactions)
                    if s == '__Other__':
                        s = 'Other reactions (not covered by specified annotations)'
                        if numReactions == 0:
                            continue
                    strings += ['# '+s, '']
                    if numReactions == 0:
                        strings += ['# NB: no reactions found for this annotation!']
                    else:
                        for reaction in sreactions:
                            thisReactantsString = ' + '.join(reaction['reactants'])
                            thisRateString = ' ->{'+str(reaction['rate'])+'} '
                            thisProductsString = ' + '.join(reaction['products'])
                            strings += ['reaction ' + thisReactantsString + thisRateString + thisProductsString]
                    if idx < (numSections - 1):
                        strings += ['']
            strings += ['']
        if self.rateDefinitions != {}:
            strings += ['#', '# Rate constant definitions', '#', '']
            for (name,value) in self.rateDefinitions.items():
                strings += ['rate ' + name + ' = ' + str(value)]
            strings += ['']
        if self.speciesInits != {}:
            strings += ['#']
            strings += ['# Species initial concentrations'] if includeZeroInits else ['# Species initial concentrations (all others zero)']
            strings += ['#', '']
            for (species,value) in self.speciesInits.items():
                if includeZeroInits or value != 0:
                    strings += ['init ' + species + ' = ' + str(value)]
            strings += ['']
        if self.perturbations != []:
            strings += ['#', '# Perturbations', '#', '']
            for p in self.perturbations:
                for a in p['actions']:
                    if a['relative']:
                        if a['amount'] < 0:
                            opStr = '-='
                            valStr = str(-1 * a['amount'])
                        else:
                            opStr = '+='
                            valStr = str(a['amount'])
                    else:
                        opStr = '='
                        valStr = str(a['amount'])
                    strings += ['perturbation ' + a['speciesName'] + ' ' + opStr + ' ' + valStr + ' @ ' + str(p['time'])]
                for f in p['functions']:
                    warningMsg = 'AN ARBITRARY CODE PERMUTATION WAS OMITTED FROM ASCII OUTPUT (TIME='+str(p['time'])+')'
                    strings += ['# !!!!!! WARNING - '+warningMsg+' !!!!!!']
                    warning(warningMsg)
            strings += ['']
        if self.simSettings != {}:
            strings += ['#', '# Simulation settings', '#', '']
            strings += ['simulation length ' + str(self.simSettings['length'])]
            strings += ['simulation points ' + str(self.simSettings['points'])]
            strings += ['simulation stiff ' + ('true' if self.simSettings['stiff'] else 'false')]
            strings += ['simulation seed ' + str(self.simSettings['seed'])]
            strings += ['simulation rtol ' + str(self.simSettings['rtol'])]
            strings += ['simulation atol ' + str(self.simSettings['atol'])]
            strings += ['simulation recordAllIfStochastic ' + str(self.simSettings['recordAllIfStochastic'])]
        return os.linesep.join(strings)

    # Write an ASCII text representation of the model to a file
    def writeASCIIToFile(self, outfile, title=None, includeZeroInits=False, reaction_annotation_sections=None):
        with open(outfile, 'w') as f:
            f.write(self.getASCII(title=title, includeZeroInits=includeZeroInits, reaction_annotation_sections=reaction_annotation_sections))
        message('Wrote parseable description of model to '+outfile)
        
    ######################################################################################

    #
    # LaTeX representation of model
    #

    # Produce a LaTeX representation of the model, that could be compiled
    def getLaTeX(self, title=None, includeZeroInits=False, reaction_annotation_sections=None):
        strings = [r'\documentclass[11pt]{article}',
                   r'\usepackage{amssymb}',
                   r'\usepackage{amsmath}',
                   r'\usepackage{mathpazo}',
                   r'\usepackage[hmargin=1in,vmargin=1in]{geometry}',
                   r'\newcommand{\reaction}[1]{\ensuremath{\mathbin{\,\xrightarrow{\ensuremath{\,#1\,}}\,}}}',
                   r'\begin{document}']
        if title is not None:
            strings += [r'\begin{center}',
                        r'{\bf\LARGE '+title+'}',
                        r'\end{center}',
                        '']
        strings += [r'\subsection*{Model size}',
                    r'\begin{itemize}',
                    r'\item Number of reactions: '+str(len(self.reactions)),
                    r'\item Number of species: '+str(len(self.speciesInits.items())),
                    r'\end{itemize}',
                    '']
        if self.reactions != []:
            strings += [r'\subsection*{Reactions}']
            if reaction_annotation_sections is None:
                strings += [r'\begin{itemize}']
                for reaction in self.reactions:
                    thisReactantsString = ' + '.join(reaction['reactants'])
                    thisRateString = r' \reaction{'+str(reaction['rate'])+'} '
                    thisProductsString = ' + '.join(reaction['products'])
                    if 'annotation' in reaction:
                        thisAnnotationString = ' --- '+reaction['annotation'] if reaction['annotation'] is not None else ''
                    else:
                        thisAnnotationString = ''
                    strings += [r'\item $' + thisReactantsString + thisRateString + thisProductsString + '$' + thisAnnotationString]
                strings += [r'\end{itemize}', '']
            else:
                sortedByAnnotation = self.sortReactionsByAnnotations(reaction_annotation_sections)
                for (s,sreactions) in sortedByAnnotation:
                    numReactions = len(sreactions)
                    if s == '__Other__':
                        s = 'Other reactions (not covered by specified annotations)'
                        if numReactions == 0:
                            continue
                    strings += [r'\subsubsection*{'+s+'}']
                    strings += [r'\begin{itemize}']
                    if numReactions == 0:
                        strings += [r'\item \textbf{NB: no reactions found for this annotation!}']
                    else:
                        for reaction in sreactions:
                            thisReactantsString = ' + '.join(reaction['reactants'])
                            thisRateString = r' \reaction{'+str(reaction['rate'])+'} '
                            thisProductsString = ' + '.join(reaction['products'])
                            strings += [r'\item $' + thisReactantsString + thisRateString + thisProductsString + '$']
                    strings += [r'\end{itemize}', '']
        if self.rateDefinitions != {}:
            strings += [r'\subsection*{Rate constant definitions}',
                        r'\begin{itemize}']
            for (name,value) in self.rateDefinitions.items():
                strings += [r'\item $' + name + ' = ' + str(value) + '$']
            strings += [r'\end{itemize}', '']
        if self.speciesInits != {}:
            strings += [r'\subsection*{Initial concentrations}'] if includeZeroInits else [r'\subsection*{Initial concentrations (all others zero)}']
            strings += [r'\begin{itemize}']
            for (species,value) in self.speciesInits.items():
                if includeZeroInits or value != 0:
                    strings += [r'\item $['+str(species)+']_0 = ' + str(value) + '$']
            strings += [r'\end{itemize}', '']
        if self.perturbations != []:
            strings += [r'\subsection*{Perturbations}',
                        r'\begin{itemize}']
            for p in self.perturbations:
                for a in p['actions']:
                    if a['relative']:
                        if a['amount'] < 0:
                            opStr = '-='
                            valStr = str(-1 * a['amount'])
                        else:
                            opStr = '+='
                            valStr = str(a['amount'])
                    else:
                        opStr = '='
                        valStr = str(a['amount'])
                    strings += [r'\item Perturbation: ' + a['speciesName'] + ' ' + opStr + ' ' + valStr + ' at time ' + str(p['time'])]
                for f in p['functions']:
                    warningMsg = 'AN ARBITRARY CODE PERMUTATION WAS OMITTED FROM LATEX OUTPUT (TIME='+str(p['time'])+')'
                    strings += [r'\item\textbf{!!!!!! WARNING - '+warningMsg+' !!!!!!}']
                    warning(warningMsg)
            strings += [r'\end{itemize}', '']
        if self.simSettings != {}:
            strings += [r'\subsection*{Simulation settings}',
                        r'\begin{itemize}']
            strings += [r'\item Simulation length = ' + str(self.simSettings['length'])]
            strings += [r'\item Simulation points = ' + str(self.simSettings['points'])]
            strings += [r'\item Simulation stiff? ' + ('true' if self.simSettings['stiff'] else 'false')]
            strings += [r'\item Stochastic simulation seed = ' + str(self.simSettings['seed'])]
            strings += [r'\item Deterministic simulation relative tolerance = ' + str(self.simSettings['rtol'])]
            strings += [r'\item Deterministic simulation absolute tolerance = ' + str(self.simSettings['atol'])]
            strings += [r'\item Record results after every reaction if simulation is stochastic? ' + str(self.simSettings['recordAllIfStochastic'])]
            strings += [r'\end{itemize}', '']
        strings += [r'\end{document}']
        return os.linesep.join(strings)

    # Write LaTeX representation of the model to a file
    def writeLaTeXToFile(self, outfile, title=None, includeZeroInits=False, reaction_annotation_sections=None):
        with open(outfile, 'w') as f:
            f.write(self.getLaTeX(title=title, includeZeroInits=includeZeroInits, reaction_annotation_sections=reaction_annotation_sections))
        message('Wrote LaTeX description of model to '+outfile)

    # Write LaTeX representation of the model to a file and make a system call to typeset it
    def writeLaTeXToFileAndTypeset(self, outfile, latexProgram='pdflatex', title=None, includeZeroInits=False, reaction_annotation_sections=None):
        self.writeLaTeXToFile(outfile, title=title, includeZeroInits=False, reaction_annotation_sections=reaction_annotation_sections)
        trySystemCall(latexProgram+' '+outfile)
        message('Finished compiling LaTeX description of model from '+outfile+' using '+latexProgram)

##########################################################################################

#
# Model creation functions
#

# Function to get an empty model
def getEmptyModel():
    return Model()

# Parse a list of strings and create a model data structure
def parseStrings(xs, strict=False):
    m = getEmptyModel()
    psim_parse.updateFromStrings(m, xs, strict=strict)
    return m

# Parse a single string and create a model data structure
def parseString(x, strict=False):
    return parseStrings(x.splitlines(), strict=strict)

# Parse a file and create a model data structure
def parseFile(infile, strict=False, quiet=False):
    with open(infile, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
    m = parseStrings(lines, strict=strict)
    if not quiet:
        message('Loaded model definition from '+infile)
    return m

##########################################################################################
