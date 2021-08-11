# tailored-RRob
This repository contains a MATLAB implementation to reproduce the results of [1]. This includes code to load 10 channel-wise invivo 3D B1+ datasets of the human body acquired during deep breathing at 7T as described in [2] and to compute and evaluate tailored non-selective respiration specific and respiration robust kT-points pTx pulses in the human heart [1].

##### Authors:
- Christoph S. Aigner  (<christoph.aigner@ptb.de>)
- Sebastian Dietrich   (<sebastian.dietrich@ptb.de>)
- Tobias Schäffter     (<tobias.schaeffter@ptb.de>)
- Sebastian Schmitter  (<sebastian.schmitter@ptb.de>)

Usage
--------

Run script main.m: This script takes the user through two pulse designs and evaluations as described in [1]. This script shows three cases: 1) default shim setting: phase and equal magnitude set by the coil manufacturer to provide sufficient B1+ throughout the heart and the aorta, 2) tailored-RSpec: XXX and 3) tailored-RRob. The functions rely on B1 maps and ROIs available at TBA.


Contents
--------

##### Test scripts (run these):
    main.m          test script to compute and evaluate tailored-RSpec and tailored-RRob at 7T
