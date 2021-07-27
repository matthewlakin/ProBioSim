
##########################################################################################
# 
# psim_examples.py
# 
#  - Some examples to illustrate how perturbations work in probiosim
#    (defined both explicitly and via the execution of arbitrary Python code)
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

import random
import numpy as np
import psim_result
from psim_lib import message, error, warning, internalError
from psim_model import getEmptyModel, parseFile
from psim_simulate import run
try:
    import matplotlib.pyplot as plt
    MatplotlibAvailable = True
except ImportError:
    MatplotlibAvailable = False

##########################################################################################

def getYLabel(simType):
    if simType == 'deterministic':
        return 'Concentration'
    elif simType == 'stochastic':
        return 'Count'
    else:
        error('Unrecognized simulation type: '+str(simType))

def example001_perturbations_explicit():
    model = getEmptyModel()
    model.setSimulationLength(80)
    model.setSimulationPoints(801)
    model.updateFromString('reaction A ->{ka}')
    model.updateFromString('reaction B ->{kb}')
    model.updateFromString('reaction C ->{kc}')
    model.setRateConstant('ka', 1.0)
    model.setRateConstant('kb', 0.5)
    model.setRateConstant('kc', 0.25)
    model.setSpeciesInit('A', 100.0)
    model.setSpeciesInit('B', 50.0)
    model.setSpeciesInit('C', 25.0)
    for t in [10.0, 20.0, 30.0]:
        model.addPerturbationSeparately(t, 'A', 100.0, True)
        model.addPerturbationSeparately(t, 'B', 50.0,  True)
        model.addPerturbationSeparately(t, 'C', 25.0,  True)
    for label in ['A', 'B']:
        for simType in ['deterministic', 'stochastic']:
            for seed in [1,100,1000,None]:
                model.setSimulationSeed(seed)
                res = run(model, simType=simType)
                outfilebase = 'example001_perturbations_explicit_'+label+'_'+simType+'_seed='+str(seed)
                res.toFile(outfilebase+'.SIMRESULTS.txt')
                if MatplotlibAvailable:
                    plt.figure()
                    plt.plot(res.times(), res.get('A'), 'r', label='A')
                    plt.plot(res.times(), res.get('B'), 'b', label='B')
                    plt.plot(res.times(), res.get('C'), 'g', label='C')
                    plt.xlabel('Time')
                    plt.ylabel(getYLabel(simType))
                    plt.legend()
                    plt.savefig(outfilebase+'.pdf', bbox_inches='tight')
                    plt.close()
                    message('Saved plot to '+outfilebase+'.pdf')
                else:
                    warning('Skipped plotting results; matplotlib not available')
                model.writeASCIIToFile(outfilebase+'.txt')
                #model.writeLaTeXToFileAndTypeset(outfilebase+'_LaTeX.tex')
    del model
    del outfilebase
    for label in ['A', 'B']:
        for simType in ['deterministic', 'stochastic']:
            for seed in [1,100,1000,None]:
                infilebase =  'example001_perturbations_explicit_'+label+'_'+simType+'_seed='+str(seed)
                message('Loading model from: '+infilebase+'.txt...')
                model = parseFile(infilebase+'.txt')
                res = run(model, simType=simType)
                outfilebase = 'example001_perturbations_explicit_RELOADED_'+label+'_'+simType+'_seed='+str(seed)
                res.toFile(outfilebase+'.SIMRESULTS.txt')
                if MatplotlibAvailable:
                    plt.figure()
                    plt.plot(res.times(), res.get('A'), 'r', label='A')
                    plt.plot(res.times(), res.get('B'), 'b', label='B')
                    plt.plot(res.times(), res.get('C'), 'g', label='C')
                    plt.xlabel('Time')
                    plt.ylabel(getYLabel(simType))
                    plt.legend()
                    plt.savefig(outfilebase+'.pdf', bbox_inches='tight')
                    plt.close()
                    message('Saved plot to '+outfilebase+'.pdf')
                else:
                    warning('Skipped plotting results; matplotlib not available')
                model.writeASCIIToFile(outfilebase+'.txt')
                del res
                res2 = psim_result.fromFile(outfilebase+'.SIMRESULTS.txt')
                if MatplotlibAvailable:
                    plt.figure()
                    plt.plot(res2.times(), res2.get('A'), 'r', label='A')
                    plt.plot(res2.times(), res2.get('B'), 'b', label='B')
                    plt.plot(res2.times(), res2.get('C'), 'g', label='C')
                    plt.xlabel('Time')
                    plt.ylabel(getYLabel(simType))
                    plt.legend()
                    plt.savefig(outfilebase+'_RELOADEDRESULTS.pdf', bbox_inches='tight')
                    plt.close()
                else:
                    warning('Skipped plotting results; matplotlib not available')

def example002_perturbations_arbitraryCode():
    for perturbationType in ['explicit', 'arbitrarycode']:
        for simseed in [42, 99]:
            model = getEmptyModel()
            model.setSimulationLength(80)
            model.setSimulationPoints(801)
            model.updateFromString('reaction A ->{ka}')
            model.updateFromString('reaction B ->{kb}')
            model.updateFromString('reaction C ->{kc}')
            model.setRateConstant('ka', 1.0)
            model.setRateConstant('kb', 0.5)
            model.setRateConstant('kc', 0.25)
            model.setSpeciesInit('A', 100.0)
            model.setSpeciesInit('B', 50.0)
            model.setSpeciesInit('C', 25.0)
            model.setSimulationSeed(simseed)
            for t in [10.0, 20.0, 30.0]:
                if perturbationType == 'explicit':
                    model.addPerturbationSeparately(t, 'A', 100.0, True)
                    model.addPerturbationSeparately(t, 'B', 50.0,  True)
                    model.addPerturbationSeparately(t, 'C', 25.0,  True)
                elif perturbationType == 'arbitrarycode':
                    fA = lambda t, x0, g, s, adjust: adjust(x0, 'A', 100.0)
                    fB = lambda t, x0, g, s, adjust: adjust(x0, 'B', 50.0)
                    fC = lambda t, x0, g, s, adjust: adjust(x0, 'C', 25.0)
                    model.addPerturbationAsCodeMultipleFunctions(t, [fA, fB, fC])
                else:
                    internalError()
            simType = 'stochastic'
            res = run(model, simType=simType)
            outfilebase = 'example002_perturbations_arbitraryCode_'+perturbationType+'_simseed='+str(simseed)
            res.toFile(outfilebase+'.SIMRESULTS.txt')
            if MatplotlibAvailable:
                plt.figure()
                plt.plot(res.times(), res.get('A'), 'r', label='A')
                plt.plot(res.times(), res.get('B'), 'b', label='B')
                plt.plot(res.times(), res.get('C'), 'g', label='C')
                plt.xlabel('Time')
                plt.ylabel(getYLabel(simType))
                plt.legend()
                plt.savefig(outfilebase+'.pdf', bbox_inches='tight')
                plt.close()
                message('Saved plot to '+outfilebase+'.pdf')
            else:
                warning('Skipped plotting results; matplotlib not available')
            model.writeASCIIToFile(outfilebase+'.txt')    

def example003_perturbations_arbitraryCode_randomized():
    for simseed in [42, 99]:
        for perturbationseed in [8675309, 909693, 1000000, 4815162342]:
            model = getEmptyModel()
            model.setSimulationLength(100)
            model.setSimulationPoints(1001)
            model.updateFromString('reaction A ->{ka}')
            model.updateFromString('reaction B ->{kb}')
            model.updateFromString('reaction C ->{kc}')
            model.setRateConstant('ka', 1.0)
            model.setRateConstant('kb', 0.5)
            model.setRateConstant('kc', 0.25)
            model.setSpeciesInit('A', 100.0)
            model.setSpeciesInit('B', 50.0)
            model.setSpeciesInit('C', 25.0)
            model.setSimulationSeed(simseed)
            perturbationRNG = random.Random(perturbationseed)
            def rpFun(t, x0, get, set, adjust):
                (sp,amt) = perturbationRNG.choice([('A',100.0),('B',50.0),('C',25.0)])
                adjust(x0, sp, amt)
            for t in range(10, 91, 10):
                model.addPerturbationAsCode(t, rpFun)
            simType = 'stochastic'
            res = run(model, simType=simType)
            outfilebase = 'example003_perturbations_arbitraryCode_randomized_simseed='+str(simseed)+'_perturbationseed='+str(perturbationseed)
            res.toFile(outfilebase+'.SIMRESULTS.txt')
            if MatplotlibAvailable:
                plt.figure()
                plt.plot(res.times(), res.get('A'), 'r', label='A')
                plt.plot(res.times(), res.get('B'), 'b', label='B')
                plt.plot(res.times(), res.get('C'), 'g', label='C')
                plt.xlabel('Time')
                plt.ylabel(getYLabel(simType))
                plt.legend()
                plt.savefig(outfilebase+'.pdf', bbox_inches='tight')
                plt.close()
                message('Saved plot to '+outfilebase+'.pdf')
            else:
                warning('Skipped plotting results; matplotlib not available')
            model.writeASCIIToFile(outfilebase+'.txt')    

def example004_perturbations_arbitraryCode_randomized_timeDependent():
    for simseed in [42, 99]:
        for perturbationseed in [8675309, 909693, 1000000, 4815162342]:
            model = getEmptyModel()
            model.setSimulationLength(200)
            model.setSimulationPoints(2001)
            model.updateFromString('reaction A ->{ka}')
            model.updateFromString('reaction B ->{kb}')
            model.updateFromString('reaction C ->{kc}')
            model.setRateConstant('ka', 1.0)
            model.setRateConstant('kb', 0.5)
            model.setRateConstant('kc', 0.25)
            model.setSpeciesInit('A', 100.0)
            model.setSpeciesInit('B', 50.0)
            model.setSpeciesInit('C', 25.0)
            model.setSimulationSeed(simseed)
            perturbationRNG = random.Random(perturbationseed)
            def rpFun(t, x0, get, set, adjust):
                if t <= 100:
                    these_weights = [0.8, 0.1, 0.1]
                else:
                    these_weights = [0.1, 0.1, 0.8]
                (sp, amt) = perturbationRNG.choices(population=[('A',100.0),('B',50.0),('C',25.0)], weights=these_weights, k=1)[0]
                adjust(x0, sp, amt)
            for t in range(10, 191, 10):
                model.addPerturbationAsCode(t, rpFun)
            simType = 'stochastic'
            res = run(model, simType=simType)
            outfilebase = 'example004_perturbations_arbitraryCode_randomized_timeDependent_simseed='+str(simseed)+'_perturbationseed='+str(perturbationseed)
            res.toFile(outfilebase+'.SIMRESULTS.txt')
            if MatplotlibAvailable:
                plt.figure()
                plt.plot(res.times(), res.get('A'), 'r', label='A')
                plt.plot(res.times(), res.get('B'), 'b', label='B')
                plt.plot(res.times(), res.get('C'), 'g', label='C')
                plt.xlabel('Time')
                plt.ylabel(getYLabel(simType))
                plt.legend()
                plt.savefig(outfilebase+'.pdf', bbox_inches='tight')
                plt.close()
                message('Saved plot to '+outfilebase+'.pdf')
            else:
                warning('Skipped plotting results; matplotlib not available')
            model.writeASCIIToFile(outfilebase+'.txt')

def example005_perturbations_arbitraryCode_stateInspection():
    for simseed in [42, 99, 7234, 89234]:
        model = getEmptyModel()
        model.setSimulationLength(200)
        model.setSimulationPoints(2001)
        model.updateFromString('reaction A ->{ka}')
        model.updateFromString('reaction B ->{kb}')
        model.updateFromString('reaction C ->{kc}')
        model.setRateConstant('ka', 1.0)
        model.setRateConstant('kb', 0.5)
        model.setRateConstant('kc', 0.25)
        model.setSpeciesInit('A', 100.0)
        model.setSpeciesInit('B', 50.0)
        model.setSpeciesInit('C', 25.0)
        model.setSimulationSeed(simseed)
        def pFun(t, x0, get, set, adjust):
            message('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            message('Time = '+str(t))
            message('Initial state:')
            for sp in ['A', 'B', 'C']:
                message('#'+sp+' = '+str(get(x0, sp)))
            for (sp, amt) in [('A',100.0),('B',50.0),('C',25.0)]:
                adjust(x0, sp, amt)
            message('Final state:')
            for sp in ['A', 'B', 'C']:
                message('#'+sp+' = '+str(get(x0, sp)))
            message('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        for t in range(10, 191, 10):
            model.addPerturbationAsCode(t, pFun)
        simType = 'stochastic'
        res = run(model, simType=simType)
        outfilebase = 'example005_perturbations_arbitraryCode_stateInspection_simseed='+str(simseed)
        res.toFile(outfilebase+'.SIMRESULTS.txt')
        if MatplotlibAvailable:
            plt.figure()
            plt.plot(res.times(), res.get('A'), 'r', label='A')
            plt.plot(res.times(), res.get('B'), 'b', label='B')
            plt.plot(res.times(), res.get('C'), 'g', label='C')
            plt.xlabel('Time')
            plt.ylabel(getYLabel(simType))
            plt.legend()
            plt.savefig(outfilebase+'.pdf', bbox_inches='tight')
            plt.close()
            message('Saved plot to '+outfilebase+'.pdf')
        else:
            warning('Skipped plotting results; matplotlib not available')
        model.writeASCIIToFile(outfilebase+'.txt')

##########################################################################################

if __name__ == '__main__':
    example001_perturbations_explicit()
    example002_perturbations_arbitraryCode()
    example003_perturbations_arbitraryCode_randomized()
    example004_perturbations_arbitraryCode_randomized_timeDependent()
    example005_perturbations_arbitraryCode_stateInspection()

##########################################################################################
