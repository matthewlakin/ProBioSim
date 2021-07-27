
##########################################################################################
# 
# psim_seeds.py
#
#  - Utility functions for handling seed values for stochastic simulations.
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
import random
import sys
from psim_lib import message, error, distinct

##########################################################################################

def saveSeeds(fname, n, comment=None, overwrite=False, failOnDuplicate=True):
    if not overwrite:
        if os.path.exists(fname):
            error('Seed file '+str(fname)+' already exists, and overwriting was disabled')
    if type(n) is not int:
        error('Number of seeds must be a positive integer')
    if not (n > 0):
        error('Number of seeds must be a positive integer')
    random.seed()
    the_seeds = []
    for x in range(n):
        new_seed = random.randrange(sys.maxsize)
        if failOnDuplicate:
            if new_seed in the_seeds:
                error('Accidentally generated a duplicate seed - you are extremely unlucky!')
        the_seeds += [new_seed]
    with open(fname, 'w') as f:
        if comment is not None:
            f.write('# '+comment+os.linesep)
        for s in the_seeds:
            f.write(str(s) + os.linesep)
    message('Wrote '+str(n)+' randomly chosen seeds to '+fname)
    message('They were '+('' if failOnDuplicate else 'NOT ')+'checked for uniqueness!')

def loadSeeds(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    seeds = [int(l.strip()) for l in lines if not l.startswith('#')]
    return seeds

def appendSeeds(fname, n, comment=None, failOnDuplicate=True):
    if not os.path.isfile(fname):
        error('Seed file '+str(fname)+' does not exist (or is not a file), so cannot append to it')
    if type(n) is not int:
        error('Number of seeds must be a positive integer')
    if not (n > 0):
        error('Number of seeds must be a positive integer')
    random.seed()
    previous_seeds = loadSeeds(fname)
    extra_seeds = []
    for x in range(n):
        new_seed = random.randrange(sys.maxsize)
        if failOnDuplicate:
            if (new_seed in previous_seeds) or (new_seed in extra_seeds):
                error('accidentally generated a duplicate seed - you are extremely unlucky!')
        extra_seeds += [new_seed]
    with open(fname, 'a') as f:
        if comment is not None:
            f.write('# '+comment+os.linesep)
        for s in extra_seeds:
            f.write(str(s) + os.linesep)
    message('Appended '+str(n)+' randomly chosen seeds to '+fname)
    message('They were '+('' if failOnDuplicate else 'NOT ')+'checked for uniqueness!')

def seedsAreUnique(seeds):
    return distinct(seeds)

##########################################################################################
