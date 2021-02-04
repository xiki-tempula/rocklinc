RocklinC
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.com/xiki-tempula/RocklinC.svg?branch=master)](https://travis-ci.com/xiki-tempula/RocklinC)
[![codecov](https://codecov.io/gh/xiki-tempula/RocklinC/branch/master/graph/badge.svg)](https://codecov.io/gh/xiki-tempula/RocklinC/branch/master)


A python module for performing Rocklin Correction.

Installing from source
----------------------

from source. Clone the source from GitHub with::

    git clone https://github.com/xiki-tempula/rocklinc.git

then do::

    cd rocklinc
    pip install .

Uasge
-----
To compute Rocklin correction. The following data is required.

 - u: MDAnalysis.Universe object that contains `positions`, `charges`, `radii` and `dimensions`.
 - box: The dimension (Å) of the simulation box in the form of a list (e.g. [100, 100, 100]).
 - lig_netq: The total charge of the ligand.
 - protein_netq: The total charge of the rest of the system excluding ligand.
 - temp: The temperature of the simulation (K).
 - water: The water model being used. Only rocklinc.waters.TIP3P and rocklinc.waters.TIP4P are supported for now.
 - lig_selection: The MDAnalysis selection string for the ligand.
 - apbs_exe: The executable path of the APBS software.

A full automatic calculation of Rocklin correction could be performed with ::

    import rocklinc
    correction = rocklinc.RocklinCorrection(box, lig_netq, protein_netq, temp)
    correction.make_APBS_input(u, lig_selection)
    correction.run_APBS(apbs_exe=apbs_exe)
    correction.read_APBS()
    result = correction.compute()
    correction.write('correction.txt')
 
`result` is the correction energy, while the details are written in `'correction.txt'`.

Please see the full documentation from RTD. https://rocklinc.readthedocs.io/

### Copyright

Copyright (c) 2020, Zhiyi Wu


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
