# tailored-RRob
This repository contains a MATLAB implementation to reproduce the results of [1]. This includes code to load 10 channel-wise invivo 3D B1+ datasets of the human body acquired during deep breathing at 7T as described in [2] and to compute and evaluate tailored non-selective respiration specific and respiration robust kT-points pTx pulses in the human heart [1].

##### Authors:
- Christoph S. Aigner  (<christoph.aigner@ptb.de>)
- Sebastian Dietrich   (<sebastian.dietrich@ptb.de>)
- Sebastian Schmitter  (<sebastian.schmitter@ptb.de>)

Usage
--------

Run script main.m: This script takes the user through two pulse designs and evaluations as described in [1]. This script shows three cases: 1) default shim setting: phase and equal magnitude set by the coil manufacturer to provide sufficient B1+ throughout the heart and the aorta, 2) tailored-RSpec: tailored for one out of three subject-specific respiration states (inhale, intermediate and exhale) and 3) tailored-RRob: tailored for all three subject-specific respiration states. The functions rely on B1 maps and ROIs available at TBA.


Contents
--------

##### Test scripts (run these):
    main.m          TBA

Dependencies
------------
These routines were tested under MATLAB R2019a under Windows, but should also run under older versions.

The 10 channel-wise invivo B1+ datasets of the human body at 7T are available at: https://doi.org/10.6084/m9.figshare.15172899.v1 and were computed as described in [2].

The optimization of the kT-points is performed using code by Will Grissom and Zhipeng Cao ([3,4] and https://bitbucket.org/wgrissom/acptx/) who have given permission for inclusion within this package. 

Please cite appropriately.

License
-------

This software is published under GNU GPLv3. 
In particular, all source code is provided "as is" without warranty of any kind, either expressed or implied. 
For details, see the attached LICENSE.

Reference
---------

[1] Aigner, CS, Dietrich, S, Schaeffter, T, and Schmitter, S, Respiration induced B1+ changes and their impact on universal and tailored 3D kT point pulses for 7T cardiac imaging, submitted to Magn. Reson. Med. 2021

[2] Dietrich, S, Aigner, CS, Kolbitsch, C, et al. 3D Free-breathing multichannel absolute B1+ Mapping in the human body at 7T. Magn Reson Med. 2021; 85: 2552– 2567. https://doi.org/10.1002/mrm.28602

[3] Grissom, W.A., Khalighi, M.-M., Sacolick, L.I., Rutt, B.K. and Vogel, M.W. (2012), Small-tip-angle spokes pulse design using interleaved greedy and local optimization methods. Magn Reson Med, 68: 1553-1562. https://doi.org/10.1002/mrm.24165

[4] Cao, Z., Yan, X. and Grissom, W.A. (2016), Array-compressed parallel transmit pulse design. Magn. Reson. Med., 76: 1158-1169. https://doi.org/10.1002/mrm.26020

Created by Christoph S. Aigner, PTB, August 2021.
Email: christoph.aigner@ptb.de
