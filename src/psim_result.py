
##########################################################################################
# 
# psim_result.py
#
#  - Class for representing the results from ProBioSim simulations, including methods to
#    access and export the data.
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
import numpy
import os
from psim_lib import message, error, fstr, pairwise, pairwise_indexes, distinct, float_leq
from psim_model import positiveDualRailSpecies, negativeDualRailSpecies
try:
    import pandas
    PandasAvailable = True
except ImportError:
    PandasAvailable = False
try:
    from openpyxl import Workbook
    ExcelAvailable = True
except ImportError:
    ExcelAvailable = False

##########################################################################################

#
# A class to wrap the dictionary of simulation results produced by the simulator from psim_simulate.py.
# This enables some cleaner and easier lookups, especially for trying to find out the values of species at a particular time.
#

class SimulationResult(object):
    # Initialize the SimulationResult object
    def __init__(self, times, concsDict):
        # Sanity check the supplied times
        num_times = len(times)
        if not (num_times > 0):
            error('Cannot create SimulationResult object because an empty sequence of times was supplied')
        for (t1,t2) in pairwise(times):
            if not (float_leq(0, t1) and float_leq(t1, t2)):
                error('Cannot create SimulationResult object because times are not non-decreasing')
        # Check concentrations map to correct number of concentrations
        # NB: Don't check for >= 0 as sometimes numerical error could take values slightly below zero...
        if not (len(concsDict.keys()) > 0):
            error('Cannot create SimulationResult object because no timecourse data was supplied')
        for x in concsDict.keys():
            if not (len(concsDict[x]) == num_times):
                error('Cannot create SimulationResult object because different number of times ('+str(num_times)+') than entries for species '+str(x)+' ('+str(len(concsDict[x]))+')')
        # Assign inputs to object parameters
        self._times_ = numpy.array(times)
        self._concsDict_ = concsDict
    # Get index associated with specified times. Returns an error if any times are outside the simulated range.
    # Also raises an error if multiple times in the list are identical.
    # NB: returns the "after" index if the requested time falls precisely on a perturbation time.
    # NB: this could be made more efficient using a binary search, though may not be a major problem in practice.    
    def _getTimeIndexes_(self, ts):
        if not (len(ts) > 0):
            error('Cannot _getTimeIndexes_ because empty sequence of times supplied')
        if not distinct(ts):
            error('Cannot _getTimeIndexes_ because sequence of times supplied contains duplicate elements')
        final_idx = len(self.times()) - 1
        initial_time = self.timeInitial()
        final_time = self.timeFinal()
        ts_remaining = list(ts.copy())
        for t in ts_remaining:
            if not (math.isclose(t, initial_time) or math.isclose(t, final_time) or (initial_time < t < final_time)):
                error('Cannot _getTimeIndexes_ because supplied time '+str(t)+' is outside of the available range')
        time_index_pairs = pairwise_indexes(self.times())
        idxs = []
        t = ts_remaining.pop(0)
        for (idx1,idx2) in time_index_pairs:
            t1 = self.times()[idx1]
            t2 = self.times()[idx2]
            if not float_leq(t1, t2):
                error('Cannot _getTimeIndexes_ because supplied times are not non-decreasing ('+str(t1)+', '+str(t2)+')')
            new_idxs = []
            if math.isclose(t1,t2):
                while True:
                    if math.isclose(t1,t):
                        new_idxs += [idx2] # Want the second of the two perturbation pair values.
                        if ts_remaining != []:
                            t = ts_remaining.pop(0)
                        else:
                            return idxs + new_idxs
                    else:
                        break
            else:
                if not (t1 < t2):
                    internalError()
                while True:
                    if math.isclose(t1,t):
                        new_idxs += [idx1]
                        if ts_remaining != []:
                            t = ts_remaining.pop(0)
                        else:
                            return idxs + new_idxs
                    elif math.isclose(t,t2): # Only keep the second one if it is the final index (otherwise, go around again - it will get found next time)
                        if idx2 == final_idx:
                            new_idxs += [idx2]
                            if ts_remaining != []:
                                t = ts_remaining.pop(0)
                            else:
                                return idxs + new_idxs
                        else:
                            break
                    elif t1 < t < t2:
                        new_idxs += [idx1] # If t is between these indexes, return the first one
                        if ts_remaining != []:
                            t = ts_remaining.pop(0)
                        else:
                            return idxs + new_idxs
                    else:
                        break
            idxs += new_idxs
        error('requested times '+str(ts_remaining)+' outside simulated time range ('+str(self.timeInitial())+' to '+str(self.timeFinal())+')')
    # Get index associated with specified time. Returns an error if outside simulated range, using the _getTimeIndexes_ method above.
    def _getTimeIndex_(self, t):
        idxs = self._getTimeIndexes_([t])
        if not (len(idxs) == 1):
            internalError()
        return idxs[0]
    # Get times
    def times(self):
        return self._times_
    def timeInitial(self):
        return self.times()[0]
    def timeFinal(self):
        return self.times()[-1]
    # Get raw species concentrations
    def get(self, x):
        if x in self._concsDict_:
            return numpy.array(self._concsDict_[x]) # Ensure this returns a numpy array.
        else:
            error('species '+x+' not found')
    def getInitial(self, x):
        return self.get(x)[0]
    def getFinal(self, x):
        return self.get(x)[-1]
    def _getIndex_(self, x, idx):
        return self.get(x)[idx]
    def _getIndexes_(self, x, idxs):
        return numpy.array([self._getIndex_(x, idx) for idx in idxs])
    def getAt(self, x, t):
        return self._getIndex_(x, self._getTimeIndex_(t))
    def getMultiple(self, xs, t):
        idx = self._getTimeIndex_(t)
        return [self._getIndex_(x, idx) for x in xs]
    def getAtMultiple(self, x, ts):
        return self._getIndexes_(x, self._getTimeIndexes_(ts))
    def getMultipleAtMultiple(self, xs, ts):
        idxs = self._getTimeIndexes_(ts)
        return [self._getIndexes_(x, idxs) for x in xs]
    # Get dual-rail species concentrations
    def getDual(self, x):
        return self.get(positiveDualRailSpecies(x)) - self.get(negativeDualRailSpecies(x)) # These work as they are numpy arrays.
    def getInitialDual(self, x):
        return self.getInitial(positiveDualRailSpecies(x)) - self.getInitial(negativeDualRailSpecies(x))
    def getFinalDual(self, x):
        return self.getFinal(positiveDualRailSpecies(x)) - self.getFinal(negativeDualRailSpecies(x))
    def _getIndexDual_(self, x, idx):
        return self._getIndex_(positiveDualRailSpecies(x), idx) - self._getIndex_(negativeDualRailSpecies(x), idx)
    def _getIndexesDual_(self, x, idxs):
        return numpy.array([self._getIndexDual_(x, idx) for idx in idxs])
    def getAtDual(self, x, t):
        return self._getIndexDual_(x, self._getTimeIndex_(t))
    def getMultipleDual(self, xs, t):
        idx = self._getTimeIndex_(t)
        return [self._getIndexDual_(x, idx) for x in xs]    
    def getAtMultipleDual(self, x, ts):
        return self._getIndexesDual_(x, self._getTimeIndexes_(ts))
    def getMultipleAtMultipleDual(self, xs, ts):
        idxs = self._getTimeIndexes_(ts)
        return [self._getIndexesDual_(x, idxs) for x in xs]
    # Get all species names (sorted)
    def allSpecies(self):
        return sorted(self._concsDict_.keys())
    # Get all species concentrations at certain points in time, as a dict
    def getInitialAll(self):
        d = {}
        for s in self.allSpecies():
            if s in d:
                internalError()
            d[s] = self.getInitial(s)
        return d
    def getFinalAll(self):
        d = {}
        for s in self.allSpecies():
            if s in d:
                internalError()
            d[s] = self.getFinal(s)
        return d
    def _getIndexAll_(self, idx):
        d = {}
        for s in self.allSpecies():
            if s in d:
                internalError()
            d[s] = self._getIndex_(s, idx)
        return d        
    def getAtAll(self, t):
        idx = self._getTimeIndex_(t)
        return self._getIndexAll_(idx)
    # Get data in columns, to save
    def _getColumnsToSave_(self, savetimes, savespecies):
        savespecies = savespecies if savespecies is not None else self.allSpecies()
        if len(savespecies) == 0:
            error('Cannot get data to save, as supplied list of species to save is empty')
        if not distinct(savespecies):
            error('Cannot get data to save, as supplied list of species to save contains duplicate elements')
        if savetimes is None:
            columns = [self.times()] + [self.get(s) for s in savespecies]
        else:
            if len(savetimes) == 0:
                error('Cannot get data to save, as suppied list of times at which to save is empty')
            columns = [savetimes] + self.getMultipleAtMultiple(savespecies, savetimes)
        return columns
    # Save to a text file
    def toFile(self, fname, savetimes=None, savespecies=None, separator='\t', decimalPlaces=None):
        savespecies = savespecies if savespecies is not None else self.allSpecies()
        columns = self._getColumnsToSave_(savetimes, savespecies)
        rows = zip(*columns)
        formatter = (lambda x: str(x)) if decimalPlaces is None else (lambda x: fstr(x, dp=decimalPlaces))
        with open(fname, 'w') as f:
            f.write(separator.join(['Times']+savespecies) + os.linesep)
            for r in rows:
                f.write(separator.join([formatter(x) for x in r]) + os.linesep)
        message('Wrote simulation results to '+fname)
    # Save to a Pandas dataFrame, if possible
    def toPandas(self, savetimes=None, savespecies=None):
        if not PandasAvailable:
            error('Pandas library not loaded, cannot save to DataFrame')
        savespecies = savespecies if savespecies is not None else self.allSpecies()
        columns = self._getColumnsToSave_(savetimes, savespecies)
        df = pandas.DataFrame(index=columns[0])
        for (speciesname, data) in zip(savespecies, columns[1:]):
            if speciesname in df:
                internalError()
            df[speciesname] = data
        return df
    # Save to an Excel spreadsheet, if possible
    def toExcel(self, fname, savetimes=None, savespecies=None):
        if not ExcelAvailable:
            error('Openpyxl library not loaded, cannot save to DataFrame')
        savespecies = savespecies if savespecies is not None else self.allSpecies()
        columns = self._getColumnsToSave_(savetimes, savespecies)
        rows = zip(*columns)
        wb = Workbook()
        ws = wb.active
        ws.append(['Times']+savespecies)
        for r in rows:
            ws.append(r)
        wb.save(fname)
        message('Wrote simulation results to '+fname)

##########################################################################################

# Load results from a text file
def fromFile(fname, separator='\t'):
    with open(fname, 'r') as f:
        tokenized = [line.rstrip().split(separator) for line in f.readlines()]
    if not (len(tokenized) > 0):
        error('Cannot load results from file '+quote(str(fname))+' because the file contains no lines')
    first_line = tokenized[0]
    if not (first_line[0] == 'Times'):
        error('Cannot load results from file '+quote(str(fname))+' because the first line\'s first heading is not '+quote('Times'))
    orderedSpecies = first_line[1:]
    if not distinct(orderedSpecies):
        error('Cannot load results from file '+quote(str(fname))+' because the header line contains duplicate elements')
    theTimes = []
    theDict = {}
    for species in orderedSpecies:
        if species in theDict:
            internalError()
        theDict[species] = []
    for (ldx,line) in enumerate(tokenized[1:], 2):
        if not (len(line) == len(first_line)):
            error('Cannot load results from file '+quote(str(fname))+' because line '+str(ldx)+' does not contain the same number of entries as the header line')
        theTimes += [numpy.float64(line[0])]
        for (species,value) in zip(orderedSpecies, line[1:]):
            theDict[species] += [numpy.float64(value)]
    theTimes = numpy.array(theTimes)
    for species in orderedSpecies:
        theDict[species] = numpy.array(theDict[species])
    message('Loaded simulation results from '+fname)
    return SimulationResult(theTimes, theDict)

##########################################################################################
