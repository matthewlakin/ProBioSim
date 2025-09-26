
##########################################################################################
# 
# psim_simulate.py
#
#  - Functions implementing the core simulation engine.
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
import math
from scipy.integrate import solve_ivp
from psim_lib import float_leq, error, warning, message, quote
from psim_result import SimulationResult

##########################################################################################

# General-purpose Gillespie simulator.
# Samples at times if times is not None.
# If times is None, samples after every reaction.
# Not intended to be called directly: pass your model object to the "run" function below!
def gillespie(N, V, y0, tspan, times, rng, verbose=False):
    (t, tend) = tspan
    if not (t < tend):
        internalError()
    if times is not None:
        if not (len(times) > 0):
            internalError()
        for (t1,t2) in zip(times, times[1:]):
            if not (float_leq(t, t1)):
                internalError()
            if not (t1 < t2):
                internalError()
            if not (float_leq(t2, tend)):
                internalError()
        next_sample_time_idx = 0
    y = np.copy(y0)
    ts = []
    res = []
    if times is None:
        ts.append(t)
        res.append(list(y))
    while True:        
        props = V(y)
        total_prop = np.sum(props)
        if total_prop == 0.0:
            break
        if total_prop < 0.0:
            error('TOTAL PROPENSITY IS NEGATIVE: '+str(total_prop))
        dt = rng.expovariate(total_prop)
        if times is not None:
            if next_sample_time_idx is not None:
                while True:
                    if (t + dt) > times[next_sample_time_idx]:
                        if verbose:
                            message('Sampling stochastic simulation at time = '+
                                    ('{:.3f}'.format(times[next_sample_time_idx]))+
                                    ' --- simulation ends at time = '+('{:.3f}'.format(tend)))
                        ts.append(times[next_sample_time_idx])
                        res.append(list(y))
                        next_sample_time_idx += 1
                        if next_sample_time_idx >= len(times):
                            next_sample_time_idx = None
                            break
                    else:
                        break
        if (t + dt) > tend:
            break
        next_reaction_idx = rng.choices(population=range(len(props)), weights=props, k=1)
        change_to_apply = N[:,next_reaction_idx]
        change_to_apply.shape = (len(change_to_apply),)
        t += dt
        y += change_to_apply
        if times is None:
            ts.append(t)
            res.append(list(y))
    if times is not None:
        if next_sample_time_idx is not None:
            while True:
                if verbose:
                    message('Sampling stochastic simulation at time = '+
                            ('{:.3f}'.format(times[next_sample_time_idx]))+
                            ' --- simulation ends at time = '+('{:.3f}'.format(tend)))
                ts.append(times[next_sample_time_idx])
                res.append(list(y))
                next_sample_time_idx += 1
                if next_sample_time_idx >= len(times):
                    next_sample_time_idx = None
                    break
    else:
        ts.append(tend)
        res.append(list(y))
    return (np.array(ts), np.transpose(np.array(res)))

##########################################################################################

# Perturbation format: {'time':float, 'actions':[action], 'functions':[function]}
# Action format: {'species':int, 'amount':float, 'relative':bool}
def simulateWithPerturbations(raw_times_list, x0, species_list, reactions_list, rate_names_values_list, perturbations=[],
                              simType='deterministic', verbose=False, stiff=True, rngSeed=None, rtol=1e-8, atol=1e-8,
                              recordAllIfStochastic=False):
    # Check that simType is valid
    if simType not in ['deterministic', 'stochastic']:
        error('In simulateWithPerturbations, the simulation type '+quote(str(simType))+' is not valid (must be "deterministic" or "stochastic")')
    # If simulation is stochastic, do some checks on the species initial values and perturbation values
    if simType == 'stochastic':
        for initval in x0:
            if not valueRepresentsAnInteger(initval):
                error('Non-integer species initial value in stochastic simulation: '+str(initval))
        for p in perturbations:
            for action in p['actions']:
                amt = action['amount']
                if not valueRepresentsAnInteger(amt):
                    error('Non-integer species perturbation amount in stochastic simulation: '+str(amt))
    # Figure out simulation time spans, factoring in perturbations etc
    start_time = raw_times_list[0]
    end_time = raw_times_list[-1]
    # Ensure that results are always output at the time a perturbation occurs, even if off the sampling grid.
    # This is required to ensure that perturbations and their results are sampled correctly within the simulator.
    for p in perturbations:
        if p['time'] not in raw_times_list and float_leq(p['time'], end_time):
            raw_times_list.append(p['time'])
    times = np.array(sorted(raw_times_list))
    # Filter any perturbations that occur at or after the end time, and sort the remaining perturbations
    for p in perturbations:
        if float_leq(end_time, p['time']):
            warning('Perturbation omitted at time '+str(p['time'])+ ' (simulation end time = '+str(end_time)+')')
    perturbations = [p for p in perturbations if p['time'] < end_time]
    perturbations = sorted(perturbations, key=lambda d: d['time'])
    # Figure out the correct time spans to run simulations over
    if perturbations == []:
        t_spans = [(start_time, end_time)]
        timess = [times]
        perturbations = [None]
    else:
        t_spans = []
        timess = []
        this_span_start = start_time
        for p in perturbations:
            if not (p['time'] < end_time):
                internalError()
            this_span_end = p['time']
            if not (float_leq(start_time, this_span_start)):
                internalError()
            if not (this_span_start < this_span_end):
                internalError()
            if not (float_leq(this_span_end, end_time)):
                internalError()
            t_spans += [(this_span_start, this_span_end)]
            timess += [[t for t in times if float_leq(this_span_start, t) and float_leq(t, this_span_end)]]
            this_span_start = this_span_end
        if this_span_start < end_time:
            t_spans += [(this_span_start, end_time)]
            timess += [[t for t in times if float_leq(this_span_start, t) and float_leq(t, end_time)]]
            perturbations += [None]
        elif math.isclose(this_span_start, end_time):
            pass
        else:
            internalError()
    # Create data structures to describe the mass action model
    N = createStoichiometryMatrix(species_list, reactions_list)
    V = createRateFunction(species_list, reactions_list, rate_names_values_list)
    params = [rate_value for (rate_name, rate_value) in rate_names_values_list]
    # These helper functions simplify the implementation of perturbations using arbitrary code.
    # The difficulty is that the simulator uses indices to refer to species internally.
    # To shield the user from this implementation detail, the simulator passes these functions,
    # which are preloaded with species_list which contains the correct mapping.
    # These can then be used to get, set, and adjust species values.
    # NB: these could perhaps be abstracted into a single object?
    def getSpeciesValue(currState, speciesName):
        thisSpeciesIdx = species_list.index(speciesName)
        return currState[thisSpeciesIdx]
    def setSpeciesValue(currState, speciesName, amount):
        thisSpeciesIdx = species_list.index(speciesName)
        currState[thisSpeciesIdx] = amount
    def adjustSpeciesValue(currState, speciesName, amount):
        thisSpeciesIdx = species_list.index(speciesName)
        currState[thisSpeciesIdx] += amount
    # Set up the solver functions / random number generators required to run the simulations
    if simType == 'deterministic':
        # ODE solver parameters
        method = 'LSODA' if stiff else 'RK45'
        # Make derivatives function to pass into ODE solver.
        # Two separate definitions here, so we don't check the value of "verbose" every single time!
        if verbose:
            def odeFun(t, y):
                message('Deterministic simulation at time = '+('{:.3f}'.format(t))+' --- simulation ends at time = '+('{:.3f}'.format(end_time)))
                return np.array((np.matmul(N, V(y))).T)[0,:]
        else:
            def odeFun(t, y):
                return np.array((np.matmul(N, V(y))).T)[0,:]
    elif simType == 'stochastic':
        # Set up random number generator using the supplied seed (for stochastic simulations only)
        rng = random.Random(rngSeed)
    else:
        internalError()
    # Actually run the simulation
    all_t = None
    all_y = None
    if not (len(t_spans) == len(timess)):
        internalError()
    if not (len(t_spans) == len(perturbations)):
        internalError()
    for (times, t_span, perturbation) in zip(timess, t_spans, perturbations):
        if simType == 'deterministic':
            this_res = solve_ivp(odeFun, t_span, x0, t_eval=times, rtol=rtol, atol=atol, method=method)
            this_res_t = this_res.t
            this_res_y = this_res.y
        elif simType == 'stochastic':
            if recordAllIfStochastic:
                times = None # This tells the "gillespie" function to record the state after every single reaction
            this_res = gillespie(N, V, x0, t_span, times, rng, verbose=verbose)
            this_res_t = this_res[0]
            this_res_y = this_res[1]
        else:
            internalError()
        all_t = this_res_t if all_t is None else np.append(all_t, this_res_t, axis=0)
        all_y = this_res_y if all_y is None else np.append(all_y, this_res_y, axis=1)
        if perturbation is not None:
            ptime = t_span[1]
            if verbose:
                message('Applying perturbation(s) at time '+str(ptime))
            x0 = np.copy(all_y[:,-1])
            for action in perturbation['actions']:
                if action['relative']:
                    x0[action['species']] += action['amount']
                else:
                    x0[action['species']] = action['amount']
            for f in perturbation['functions']:
                f(ptime, x0, getSpeciesValue, setSpeciesValue, adjustSpeciesValue)
            for idx in range(len(x0)):
                x0[idx] = max(x0[idx], 0) # We now do this for ALL species, as with arbitrary perturbations we can't know which species were modified...
            if simType == 'stochastic':
                for idx in range(len(x0)):
                    x0[idx] = int(x0[idx]) # During a stochastic simulation, recast all state values to ints (this will round _down_)...
    # Return all times and all solution values (these are numpy arrays)
    return (all_t, all_y)

##########################################################################################

# Compute the stoichiometry of a species for a given reaction
def stoichiometry(species, reaction):
    return reaction['products'].count(species) - reaction['reactants'].count(species)

# Create a stoichiometry matrix for given species and reactions
def createStoichiometryMatrix(species_list, reactions_list):
    num_species = len(species_list)
    num_reactions = len(reactions_list)
    N = np.zeros((num_species, num_reactions), dtype=np.int64)
    for (row_idx, species) in enumerate(species_list):
        for (col_idx, reaction) in enumerate(reactions_list):
            N[row_idx][col_idx] = stoichiometry(species, reaction)
    return N

# Get index of species from species list
def getSpeciesIndex(species, species_list):
    return species_list.index(species)

# Helper functions for creation of rate function
def getRateConstantValue(reaction, rate_names_values_list):
    if type(reaction['rate']) in [int, float]:
        return float(reaction['rate'])
    else:
        for (rate_name, rate_value) in rate_names_values_list:
            if rate_name == reaction['rate']:
                return rate_value
    error('Rate constant type not recognized: '+str(reaction['rate']))
def getReactantIndexes(reaction, species_list):
    return [getSpeciesIndex(reactant, species_list) for reactant in reaction['reactants']]
def lookUpIndexesAndSum(y, indexes):
    ans = 1.0
    for idx in indexes:
        ans *= y[idx]
    return ans

# Create an efficient rate function for given species and reactions.
# Efficiency is enhanced by looking up the indexes corresponding to the reactants
# in each reaction and storing these. Rate constants are pre-looked up and stored.
# Then, in the inner loop, just have to look up these values and multiply them
# together to make the multiplied concentrations vector, and multiply that by the
# rate constants vector to get V(y).
def createRateFunction(species_list, reactions_list, rate_names_values_list):
    rate_constants_vector = np.array([[getRateConstantValue(reaction, rate_names_values_list)] for reaction in reactions_list])
    reactant_indexes_vectors = [getReactantIndexes(reaction, species_list) for reaction in reactions_list]
    def V(y):
        multiplied_concs_vector = np.array([[lookUpIndexesAndSum(y,indexes)] for indexes in reactant_indexes_vectors])
        return rate_constants_vector * multiplied_concs_vector
    return V

##########################################################################################

# Check if value represents an integer
def valueRepresentsAnInteger(n):
    if isinstance(n, (int, np.int32, np.int64, np.integer)):
        return True
    elif isinstance(n, (float, np.float32, np.float64)):
        return n.is_integer()
    else:
        return False

# Fix perturbation to use indexes rather than names
def changePerturbationNamesToIndexes(p, species_ordering):
    new_actions = [{'species':getSpeciesIndex(a['speciesName'], species_ordering),
                    'amount':a['amount'],
                    'relative':a['relative']}
                   for a in p['actions']]
    return {'time':p['time'], 'actions':new_actions, 'functions':p['functions']}

# Prepare data structures and run simulation.
# This version takes a kwarg to specify simulation type ('deterministic' or 'stochastic').
# This function does some preprocessing to convert species names into indexes.
# Everything else is handled with the "simulateWithPerturbations" function.
def run(model, simType='deterministic', verbose=False):
    if simType not in ['deterministic', 'stochastic']:
        error('In psim_simulate.run, the simulation type '+quote(str(simType))+' is not valid (must be "deterministic" or "stochastic")')
    species_list = model.allSpeciesNames
    reactions_list = model.reactions
    rate_names_values_list = list(model.rateDefinitions.items())
    x0 = np.array([model.speciesInits[x] for x in species_list])
    perturbations = [changePerturbationNamesToIndexes(p, species_list) for p in model.perturbations]
    raw_times_list = list(np.linspace(0, model.simSettings['length'], model.simSettings['points']))
    (t, y) = simulateWithPerturbations(raw_times_list, x0, species_list, reactions_list, rate_names_values_list,
                                       perturbations=perturbations, simType=simType, verbose=verbose,
                                       stiff=model.simSettings['stiff'], rngSeed=model.simSettings['seed'],
                                       rtol=model.simSettings['rtol'], atol=model.simSettings['atol'],
                                       recordAllIfStochastic=model.simSettings['recordAllIfStochastic'])
    concs = {}
    for (idx,x) in enumerate(species_list):
        if x in concs:
            internalError()
        concs[x] = y[idx]
    return SimulationResult(t, concs)

# Prepare data structures and run a DETERMINISTIC (ODE) simulation.
def runDeterministic(model, verbose=False):
    return run(model, simType='deterministic', verbose=verbose)

# Prepare data structures and run a STOCHASTIC (GILLESPIE) simulation.
def runStochastic(model, verbose=False):
    return run(model, simType='stochastic', verbose=verbose)

##########################################################################################
