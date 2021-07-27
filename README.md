# ProBioSim

ProBioSim is a simulator for mass-action chemical reaction networks.

## Description

- ProBioSim includes both deterministic (ODE-based) and stochastic simulation engines.
- Most notably, ProBioSim provides powerful features for modifying the state of the simulation at specified timepoints, by specifying these modifications as arbitrary Python functions.
- This enables state- and time-dependent modifications to the simulation, which can be used to model external interventions in the evolution of the system by a responsive internal environment.

## Getting Started

### Dependencies

- ProBioSim should work on any platform where Python3 runs.
- The following libraries are required: numpy, scipy.
- The following libraries are optional: matplotlib, openpyxl, pandas.

### Installing

- To install, simply extract the source files into a folder of your choice, then add that folder to your $PYTHONPATH environment variable.

## Help

- See the doc/probiosim-doc.pdf file for further details on how to use ProBioSim, including an example and a discussion of important aspects of the codebase.
- Executable examples are provided in the src/psim\_examples.py source file.

## Authors

- Matthew Lakin
  - Email: mlakin@cs.unm.edu
  - WWW: https://cs.unm.edu/~mlakin

## Version History

* 0.1
    * Initial release

## Acknowledgments

This material is based upon work supported by the National Science Foundation under Grants 1518861, 1525553, and 1935087.
